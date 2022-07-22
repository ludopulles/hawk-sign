/*
 * Hawk signature generation.
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2017-2019  Falcon Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ===========================(LICENSE END)=============================
 *
 * @author   Thomas Pornin <thomas.pornin@nccgroup.com>
 * @author   Ludo Pulles <ludo.pulles@cwi.nl>
 */

#include "inner.h"

// =============================================================================

/*
 * Discrete Gaussian Sampler
 *
 * The table below contains the cumulative probability table for two discrete
 * gaussian distribution. The first is a discrete gaussian on 2Z with standard
 * deviation 2 sigma, and the other is a discrete gaussian on 2Z + 1 with the
 * same standard deviation.  For n = 512, we have sigma = 1.278.
 *
 * The elements in the even indices of the table (starting from zero) contain
 * the probabilities P(|X| >= 2), P(|X| >= 4), etc. when X is sampled from
 * D_{2Z, 2sigma}, and the number is scaled by a factor of 2^{78}, and then the
 * highest 15 bits are stored in gauss_hi while the lowest 63 bits are stored
 * in gauss_lo. Similarly the odd indices in the table contain P(|X| >= 3),
 * P(|X| >= 5), etc.
 *
 * The generation is contant-time so the whole probability is read fully such
 * that the time consumption does not depend on the coset that is sampled from
 * nor on the outcome of the sampling.
 *
 * To generate the values, run `sage code/renyi.sage`.
 */

/*
 * Precomputed CDT with 78 bits of precision.
 * RD_{513}(sign_sampler, mu=0/2) = 1 + 1.334079E-24 < 1 + 2^-79
 * RD_{513}(sign_sampler, mu=1/2) = 1 + 1.449739E-24 < 1 + 2^-79
 */
static const uint16_t gauss_hi_512[10] = {
	0x580B, 0x35F9,
	0x1D34, 0x0DD7,
	0x05B7, 0x020C,
	0x00A2, 0x002B,
	0x000A, 0x0001,
};
static const uint64_t gauss_lo_512[26] = {
	0x0C27920A04F8F267, 0x3C689D9213449DC9,
	0x1C4FF17C204AA058, 0x7B908C81FCE3524F,
	0x5E63263BE0098FFD, 0x4EBEFD8FF4F07378,
	0x56AEDFB0876A3BD8, 0x4628BC6B23887196,
	0x061E21D588CC61CC, 0x7F769211F07B326F,
	0x2BA568D92EEC18E7, 0x0668F461693DFF8F,
	0x00CF0F8687D3B009, 0x001670DB65964485,
	0x000216A0C344EB45, 0x00002AB6E11C2552,
	0x000002EDF0B98A84, 0x0000002C253C7E81,
	0x000000023AF3B2E7, 0x0000000018C14ABF,
	0x0000000000EBCC6A, 0x000000000007876E,
	0x00000000000034CF, 0x000000000000013D,
	0x0000000000000006, 0x0000000000000000,
};

/*
 * Precomputed CDT with 78 bits of precision.
 * RD_{513}(sign_sampler, mu=0/2) = 1 + 3.446518E-25 < 1 + 2^-81
 * RD_{513}(sign_sampler, mu=1/2) = 1 + 2.235739E-24 < 1 + 2^-78
 */
static const uint16_t gauss_hi_1024[10] = {
	0x58B0, 0x36FE,
	0x1E3A, 0x0EA0,
	0x0632, 0x024A,
	0x00BC, 0x0034,
	0x000C, 0x0002,
};
static const uint64_t gauss_lo_1024[26] = {
	0x3AAA2EB76504E560, 0x01AE2B17728DF2DE,
	0x70E1C03E49BB683E, 0x6A00B82C69624C93,
	0x55CDA662EF2D1C48, 0x2685DB30348656A4,
	0x31E874B355421BB7, 0x430192770E205503,
	0x57C0676C029895A7, 0x5353BD4091AA96DB,
	0x3D4D67696E51F820, 0x09915A53D8667BEE,
	0x014A1A8A93F20738, 0x0026670030160D5F,
	0x0003DAF47E8DFB21, 0x0000557CD1C5F797,
	0x000006634617B3FF, 0x0000006965E15B13,
	0x00000005DBEFB646, 0x0000000047E9AB38,
	0x0000000002F93038, 0x00000000001B2445,
	0x000000000000D5A7, 0x00000000000005AA,
	0x0000000000000021, 0x0000000000000000,
};

static inline int8_t
mkgauss_sign_512(prng *rng, uint8_t parity)
{
	uint16_t r_hi, p_hi;
	uint64_t r_lo, p_lo;
	uint8_t c, v, k, neg;

	/*
	 * We use 80 random bits to determine the value, by looking in the
	 * cumulative probability table. However, we only use 15 bits for r_hi so
	 * we can check if r_hi < p_hi holds by computing (r_hi - p_hi) >> 15. This
	 * way is better than doing a comparison to achieve constant-time
	 * execution.
	 */
	prng_get_80_bits(rng, &r_hi, &r_lo);

	/*
	 * Get the sign bit out of the lowest part
	 */
	neg = (uint8_t)(r_lo >> 63);

	/*
	 * Unset the sign bits in the unsigned ints for convenience in comparisons
	 * later on, as we can now use the highest bit of `a - b` to check if `a <
	 * b` or not for numbers `a, b`.
	 */
	r_hi &= ~((uint16_t)1u << 15);
	r_lo &= ~((uint64_t)1u << 63);

	v = 0;
	for (k = 10; k < 26; k += 2) {
		/*
		 * Constant-time for:
		 *     p_lo = gauss_lo[k + parity];
		 */
		p_lo = (gauss_lo_512[k] & (parity - 1u)) | ((gauss_lo_512[k + 1] & -parity));

		/*
		 * Add 1 iff r_lo < p_lo.
		 */
		v += (uint8_t)((uint64_t)(r_lo - p_lo) >> 63);
	}

	/*
	 * If r_hi > 0, set v to zero, otherwise leave v as is. This is a
	 * micro-optimization as p_hi would be zero for all k >= 5.
	 */
	v = v & -((r_hi - 1) >> 15);

	for (k = 0; k < 10; k += 2) {
		/*
		 * Constant-time for:
		 *     p_lo = gauss_lo_512[k + parity];
		 *     p_hi = gauss_hi_512[k + parity];
		 */
		p_lo = (gauss_lo_512[k] & (parity - 1u)) | ((gauss_lo_512[k + 1] & -parity));
		p_hi = (gauss_hi_512[k] & (parity - 1u)) | ((gauss_hi_512[k + 1] & -parity));

		/*
		 * c = [[ r_lo < p_lo ]]
		 */
		c = (uint8_t)((uint64_t)(r_lo - p_lo) >> 63);

		/*
		 * Constant-time code to add 1 to v iff
		 *     r_hi < p_hi or (r_hi == p_hi and c is true)
		 * holds.
		 */
		c = (uint8_t)((uint16_t)(r_hi - p_hi - c) >> 15);
		v += c;
	}

	/*
	 * Multiply by two and apply the change in support:
	 * If parity = 0, then v = 0,2,4,...
	 * If parity = 1, then v = 1,3,5,...
	 */
	v = (v << 1) | parity;

	/*
	 * Apply the sign ('neg' flag). If neg = 0, this has no effect.
	 * However, if neg = 1, this changes v into -v = (~v) + 1.
	 */
	v = (v ^ -neg) + neg;
	return *(int8_t *)&v;
}

static inline int8_t
mkgauss_sign_1024(prng *rng, uint8_t parity)
{
	uint16_t r_hi, p_hi;
	uint64_t r_lo, p_lo;
	uint8_t c, v, k, neg;

	/*
	 * We use 80 random bits to determine the value, by looking in the
	 * cumulative probability table. However, we only use 15 bits for r_hi so
	 * we can check if r_hi < p_hi holds by computing (r_hi - p_hi) >> 15. This
	 * way is better than doing a comparison to achieve constant-time
	 * execution.
	 */
	prng_get_80_bits(rng, &r_hi, &r_lo);

	/*
	 * Get the sign bit out of the lowest part
	 */
	neg = (uint8_t)(r_lo >> 63);

	/*
	 * Unset the sign bits in the unsigned ints for convenience in comparisons
	 * later on, as we can now use the highest bit of `a - b` to check if `a <
	 * b` or not for numbers `a, b`.
	 */
	r_hi &= ~((uint16_t)1u << 15);
	r_lo &= ~((uint64_t)1u << 63);

	v = 0;
	for (k = 10; k < 26; k += 2) {
		/*
		 * Constant-time for:
		 *     p_lo = gauss_lo[k + parity];
		 */
		p_lo = (gauss_lo_1024[k] & (parity - 1u)) | ((gauss_lo_1024[k + 1] & -parity));

		/*
		 * Add 1 iff r_lo < p_lo.
		 */
		v += (uint8_t)((uint64_t)(r_lo - p_lo) >> 63);
	}

	/*
	 * If r_hi > 0, set v to zero, otherwise leave v as is. This is a
	 * micro-optimization as p_hi would be zero for all k >= 5.
	 */
	v = v & -((r_hi - 1) >> 15);

	for (k = 0; k < 10; k += 2) {
		/*
		 * Constant-time for:
		 *     p_lo = gauss_lo_1024[k + parity];
		 *     p_hi = gauss_hi_1024[k + parity];
		 */
		p_lo = (gauss_lo_1024[k] & (parity - 1u)) | ((gauss_lo_1024[k + 1] & -parity));
		p_hi = (gauss_hi_1024[k] & (parity - 1u)) | ((gauss_hi_1024[k + 1] & -parity));

		/*
		 * c = [[ r_lo < p_lo ]]
		 */
		c = (uint8_t)((uint64_t)(r_lo - p_lo) >> 63);

		/*
		 * Constant-time code to add 1 to v iff
		 *     r_hi < p_hi or (r_hi == p_hi and c is true)
		 * holds.
		 */
		c = (uint8_t)((uint16_t)(r_hi - p_hi - c) >> 15);
		v += c;
	}

	/*
	 * Multiply by two and apply the change in support:
	 * If parity = 0, then v = 0,2,4,...
	 * If parity = 1, then v = 1,3,5,...
	 */
	v = (v << 1) | parity;

	/*
	 * Apply the sign ('neg' flag). If neg = 0, this has no effect.
	 * However, if neg = 1, this changes v into -v = (~v) + 1.
	 */
	v = (v ^ -neg) + neg;
	return *(int8_t *)&v;
}

static inline int8_t
mkgauss_sign(prng *rng, uint8_t parity, unsigned logn)
{
	return (logn == 10)
		? mkgauss_sign_1024(rng, parity)
		: mkgauss_sign_512(rng, parity);
}

// =============================================================================

/*
 * The hash of a message is a point in {0,1}^{2n}, so it consists of 2n bits,
 * which is n/4 bytes.
 * Therefore, the first component ranges from byte 0 to byte n / 8 - 1,
 * and the second component ranges from byte n / 8 to n / 4 - 1.
 * However, if logn <= 3, the second component is just byte 1.
 */
#define SECOND_HASH(h, logn) \
	((h) + ((logn) <= 3 ? 1u : 1u << ((logn) - 3)))

/*
 * When flag = 1, this does nothing.
 * When flag = 0, this replaces s by h - s.
 */
static void
conditional_flip(uint16_t flag, int16_t *s, const uint8_t *h, unsigned logn)
{
	size_t n, u, v;
	int16_t value;
	uint8_t hash;

	n = MKN(logn);

	if (logn <= 3) {
		for (u = 0; u < n; u ++) {
			value = (int16_t)((h[0] >> u) & 1) - 2 * s[u];
			s[u] += value & (flag - 1);
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				value = (int16_t)(hash & 1) - 2 * s[u];
				s[u] += value & (flag - 1);
				hash >>= 1;
			}
		}
	}
}

static void
hash_to_fft(fpr *p, const uint8_t *h, unsigned logn)
{
	size_t n, u, v;
	uint8_t hash;

	n = MKN(logn);

	if (logn <= 3) {
		for (u = 0; u < n; u ++) {
			p[u] = fpr_of((h[0] >> u) & 1);
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				p[u] = fpr_of((hash >> v) & 1);
			}
		}
	}

	Zf(FFT)(p, logn);
}

/*
 * Returns the vector s = (h - noise) / 2 which gives a lattice point s that is
 * close to h / 2, assuming that noise has a small norm. This function assumes
 * that noise is in coefficient representation and integral (up to rounding
 * errors).
 */
static void
noise_to_lattice(int16_t *s, const uint8_t *h, const fpr *noise, unsigned logn)
{
	size_t n, u, v;
	uint8_t hash;
	int64_t x;

	n = MKN(logn);
	if (logn <= 3) {
		for (u = 0; u < n; u ++) {
			x = fpr_rint(noise[u]);
			s[u] = (int16_t)((int64_t)((h[0] >> u) & 1) - x) / 2;
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				x = fpr_rint(noise[u]);
				s[u] = (int16_t)((int64_t)(hash & 1) - x) / 2;
				hash >>= 1;
			}
		}
	}
}

static inline void
construct_basis(
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	fpr *restrict bf, fpr *restrict bg, fpr *restrict bF, fpr *restrict bG,
	unsigned logn)
{
	Zf(int8_to_fft)(bf, f, logn);
	Zf(int8_to_fft)(bg, g, logn);
	Zf(int8_to_fft)(bF, F, logn);

	if (G == NULL) {
		size_t u, hn;

		/*
		 * Compute (in FFT representation) G = (1 + gF) / f.
		 */
		hn = MKN(logn) >> 1;
		Zf(poly_prod_fft)(bG, bg, bF, logn);
		for (u = 0; u < hn; u++) {
			bG[u] = fpr_add(bG[u], fpr_one);
		}
		Zf(poly_div_fft)(bG, bf, logn);
	} else {
		Zf(int8_to_fft)(bG, G, logn);
	}
}

/*
 * Sample a vector (x0, x1) congruent to (t0, t1) = B * (h0, h1) (mod 2) from a
 * discrete Gaussian with lattice coset 2Z^{2n} + (t0, t1) and standard
 * deviation 1.278. Returns whether or not (x0, x1) has a squared l2-norm less
 * than the allowed bound for a signature.
 *
 * If this squared l2-norm is less than the allowed bound, the values in
 * (x0, x1) are put in FFT representation.
 */
static int
sample_short(prng *rng, fpr *restrict x0, fpr *restrict x1,
	const fpr *restrict bf, const fpr *restrict bg,
	const fpr *restrict bF, const fpr *restrict bG,
	const uint8_t *restrict h, unsigned logn)
{
	size_t n, u;
	int32_t norm, z;

	n = MKN(logn);
	norm = 0;

	hash_to_fft(x0, h, logn);
	hash_to_fft(x1, SECOND_HASH(h, logn), logn);

	/*
	 * Set the target vector to (t0, t1) = B * (h0, h1), i.e.:
	 *     t0 = f h0 + F h1,
	 *     t1 = g h0 + G h1.
	 */
	Zf(poly_matmul_fft)(bf, bF, bg, bG, x0, x1, logn);

	/*
	 * Sample and write the result in (x0, x1). Gaussian smoothing is used to
	 * not reveal information on the secret basis.
	 */
	Zf(iFFT)(x0, logn);
	Zf(iFFT)(x1, logn);

	for (u = 0; u < n; u ++) {
		z = fpr_rint(x0[u]) & 1;
		z = mkgauss_sign(rng, z, logn);
		x0[u] = fpr_of(z);
		norm += z*z;
	}
	for (u = 0; u < n; u ++) {
		z = fpr_rint(x1[u]) & 1;
		z = mkgauss_sign(rng, z, logn);
		x1[u] = fpr_of(z);
		norm += z*z;
	}

	Zf(FFT)(x0, logn);
	Zf(FFT)(x1, logn);

	/*
	 * Test whether the l2-norm of (x0, x1) is below the given bound. The
	 * code below uses only 32-bit operations to compute the squared norm,
	 * since the max. value is 2n * 128^2 <= 2^24 (when logn <= 9).
	 * For a large enough verification margin, it is unlikely that the
	 * norm of the gaussian (x0, x1) is too large.
	 */
	return (uint32_t)norm <= L2BOUND(logn);
}


/*
 * The following are different sign functions for uncompressed HAWK or HAWK.
 * If a lattice point is generated that is too far away from (h0, h1) / 2, s0
 * and s1 are untouched and 0 is returned; the caller should then try again.
 * Otherwise, 1 is returned and (s0, s1) contains a valid signature for (h0,
 * h1), except for a miniscule probability in HAWK that decompression gives a
 * different s0.
 *
 * All signing functions use a fast PRNG for gaussian sampling during signing,
 * that is seeded with the SHAKE256 context.
 */

/* see inner.h */
int
Zf(uncompressed_sign)(inner_shake256_context *rng,
	int16_t *restrict s0, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	size_t n;
	fpr *x0, *x1, *bf, *bg, *bF, *bG;
	uint16_t flag;
	int norm_okay;
	prng p;

	n = MKN(logn);
	bf = (fpr *)tmp;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	x0 = bG + n;
	x1 = x0 + n;

	Zf(prng_init)(&p, rng);

	construct_basis(f, g, F, G, bf, bg, bF, bG, logn);

	norm_okay = sample_short(&p, x0, x1, bf, bg, bF, bG, h, logn);

	/*
	 * Compute (s0, s1) = ((h0, h1) - B^{-1} (x0, x1)) / 2, so
	 *
	 *     s0 = (h0 - (x0 * G + x1 (-F))) / 2,
	 *     s1 = (h1 - (x0 * (-g) + x1 f)) / 2.
	 */
	Zf(poly_neg)(x0, logn);
	Zf(poly_matmul_fft)(bG, bF, bg, bf, x0, x1, logn);
	Zf(poly_neg)(x0, logn);
	Zf(iFFT)(x0, logn);
	Zf(iFFT)(x1, logn);

	noise_to_lattice(s0, h, x0, logn);
	noise_to_lattice(s1, SECOND_HASH(h, logn), x1, logn);

	flag = (uint16_t)Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn);
	conditional_flip(flag, s0, h, logn);
	conditional_flip(flag, s1, SECOND_HASH(h, logn), logn);

	return norm_okay;
}

/* see inner.h */
int
Zf(sign_dyn)(inner_shake256_context *rng, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
#if HAWK_RECOVER_CHECK
	size_t n, u;
#else
	size_t n;
#endif
	fpr *x0, *x1, *bf, *bg, *bF, *bG;
	uint16_t flag;
	int norm_okay;
	prng p;

	n = MKN(logn);
	bf = (fpr *)tmp;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	x0 = bG + n;
	x1 = x0 + n;

	Zf(prng_init)(&p, rng);

	construct_basis(f, g, F, G, bf, bg, bF, bG, logn);

	norm_okay = sample_short(&p, x0, x1, bf, bg, bF, bG, h, logn);

#if HAWK_RECOVER_CHECK
	/*
	 * Compute *twice* the rounding error, which is given by:
	 *
	 *     (f* x0 + g* x1) / (f* f + g* g).
	 *
	 * If the above quantity is in the (-1, 1)^n box, simple rounding works for
	 * recovering s0. Otherwise, reject this signature.
	 *
	 * Currently, this check is NOT performed, as the probability of a failure
	 * to happen here is (heuristically) less than 2^{-105}.
	 */
	Zf(poly_add_muladj_fft)(bF, x0, x1, bf, bg, logn);
	Zf(poly_invnorm2_fft)(bG, bf, bg, logn);
	Zf(poly_mul_autoadj_fft)(bF, bG, logn);
	Zf(iFFT)(bF, logn);

	for (u = 0; u < n; u++) {
		if (!fpr_lt(fpr_neg(fpr_almost_one), bF[u]) || !fpr_lt(bF[u], fpr_almost_one)) {
			return 0;
		}
	}
#endif

	/*
	 * Compute s1 in (s0, s1) = ((h0, h1) - B^{-1} (x0, x1)) / 2, so
	 *
	 *     s1 = (h1 - (x0 * (-g) + x1 f)) / 2.
	 */
	Zf(poly_neg)(x0, logn);
	Zf(poly_add_mul_fft)(bF, x0, x1, bg, bf, logn);
	Zf(iFFT)(bF, logn);
	noise_to_lattice(s1, SECOND_HASH(h, logn), bF, logn);

	flag = (uint16_t)Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn);
	conditional_flip(flag, s1, SECOND_HASH(h, logn), logn);

	if (logn <= 6 && !Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn))
		return 0;
	return norm_okay;
}

/* see inner.h */
int
Zf(sign)(inner_shake256_context *rng, int16_t *restrict s1,
	const fpr *restrict expanded_seckey, const uint8_t *restrict h,
	unsigned logn, uint8_t *restrict tmp)
{
#if HAWK_RECOVER_CHECK
	size_t n, u;
	const fpr *bf, *bg, *bF, *bG, *invq00;
#else
	size_t n;
	const fpr *bf, *bg, *bF, *bG;
#endif
	fpr *x0, *x1, *res;
	uint16_t flag;
	int norm_okay;
	prng p;

	n = MKN(logn);

	bf = expanded_seckey;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;

	x0 = (fpr *)tmp;
	x1 = x0 + n;
	res = x1 + n;

	Zf(prng_init)(&p, rng);

	norm_okay = sample_short(&p, x0, x1, bf, bg, bF, bG, h, logn);

#if HAWK_RECOVER_CHECK
	/*
	 * Compute *twice* the rounding error, which is given by:
	 *
	 *     (f* x0 + g* x1) / (f* f + g* g).
	 *
	 * If the above quantity is in the (-1, 1)^n box, simple rounding works for
	 * recovering s0. Otherwise, reject this signature.
	 *
	 * Currently, this check is NOT performed, as the probability of a failure
	 * to happen here is (heuristically) less than 2^{-105}.
	 */
	invq00 = bG + n;

	Zf(poly_add_muladj_fft)(res, x0, x1, bf, bg, logn);
	Zf(poly_mul_autoadj_fft)(res, invq00, logn);
	Zf(iFFT)(res, logn);

	for (u = 0; u < n; u++) {
		if (!fpr_lt(fpr_neg(fpr_almost_one), res[u]) || !fpr_lt(res[u], fpr_almost_one)) {
			return 0;
		}
	}
#endif

	/*
	 * Compute s1 in (s0, s1) = ((h0, h1) - B^{-1} (x0, x1)) / 2, so
	 *
	 *     s1 = (h1 - (x0 * (-g) + x1 f)) / 2.
	 */
	Zf(poly_neg)(x0, logn);
	Zf(poly_add_mul_fft)(res, x0, x1, bg, bf, logn);
	Zf(iFFT)(res, logn);
	noise_to_lattice(s1, SECOND_HASH(h, logn), res, logn);

	flag = (uint16_t)Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn);
	conditional_flip(flag, s1, SECOND_HASH(h, logn), logn);

	if (logn <= 6 && !Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn))
		return 0;
	return norm_okay;
}


/******************************************************************************
 * NTT part
 *****************************************************************************/
#define FERMAT_P (65537)

/*
 * p = 2^16 + 1.
 * g = 3623 which is the smallest element g of order 2048 such that g^64 == 2.
 * ig = g^{-1} = 12174
 */

static const uint32_t GMb[1024] = {
	1, 256, 16, 4096, 4, 1024, 64, 16384, 2, 512, 32, 8192, 8, 2048, 128, 32768, 4080, 61425, 65280, 65282, 16320, 49089, 64509, 64517, 8160, 57313, 65023, 65027, 32640, 32641, 63481, 63497, 4938, 18925, 13471, 40652, 19752, 10163, 53884, 31534, 9876, 37850, 26942, 15767, 39504, 20326, 42231, 63068, 27181, 11414, 41674, 51550, 43187, 45656, 35622, 9589, 54362, 22828, 17811, 37563, 20837, 25775, 5707, 19178, 59963, 14870, 41890, 41309, 43241, 59480, 36486, 34162, 54389, 29740, 18243, 17081, 20945, 53423, 7435, 2787, 64956, 47875, 56241, 45093, 63213, 60426, 28353, 49298, 64375, 30213, 46945, 24649, 60889, 55315, 56706, 33059, 1128, 26620, 18048, 32698, 4512, 40943, 6655, 65255, 2256, 53240, 36096, 65396, 9024, 16349, 13310, 64973, 14650, 14791, 37789, 40045, 58600, 59164, 20082, 29106, 29300, 29582, 10041, 14553, 51663, 52791, 40164, 58212, 21417, 43181, 14987, 35526, 20131, 41650, 59948, 11030, 42834, 20825, 29974, 5515, 40262, 17763, 54359, 22060, 20539, 15024, 939, 43773, 16619, 60096, 3756, 44018, 41078, 30048, 1878, 22009, 33238, 54655, 7512, 22499, 45965, 35917, 14533, 50376, 52786, 12594, 58132, 4893, 26393, 6297, 29066, 35215, 40035, 25188, 50727, 9786, 35843, 628, 49192, 10048, 12298, 2512, 157, 40192, 6149, 1256, 32847, 20096, 24596, 5024, 314, 14847, 30056, 26507, 22137, 30890, 54687, 40491, 23011, 58023, 60112, 53014, 44274, 61780, 43837, 15445, 46022, 50509, 8753, 12510, 8974, 3549, 35012, 50040, 35896, 14196, 17506, 25020, 17948, 7098, 4487, 34543, 6255, 28392, 40760, 14177, 62327, 30221, 31966, 56708, 52697, 55347, 15983, 28354, 59117, 60442, 63932, 47879, 39857, 45157, 33431, 38526, 10600, 26583, 2650, 23030, 42400, 40795, 1325, 11515, 21200, 53166, 5300, 46060, 19263, 16053, 18729, 10423, 37516, 35694, 9379, 41692, 18990, 11702, 37458, 20846, 9495, 5851, 18758, 17847, 37980, 23404, 63715, 57864, 36385, 8306, 58249, 34845, 14466, 33224, 61893, 50191, 7233, 16612, 50961, 4153, 28932, 911, 11095, 22229, 46446, 27979, 44380, 23379, 54710, 46379, 22190, 44458, 27355, 55958, 23223, 46758, 43883, 27221, 47070, 56649, 32213, 54403, 57206, 29985, 63315, 21001, 28603, 47761, 64426, 43269, 48875, 59970, 61093, 42002, 4995, 33517, 14383, 11976, 19980, 2994, 57532, 47904, 9990, 1497, 28766, 23952, 39960, 5988, 49527, 30271, 63130, 39178, 27025, 37015, 55909, 25638, 42563, 16986, 60723, 12819, 54050, 8493, 46281, 51276, 19589, 33972, 23398, 26021, 46683, 23114, 28055, 38547, 55658, 26919, 46796, 52042, 27829, 46228, 56110, 11557, 45779, 53838, 41968, 61277, 16118, 62914, 36798, 48497, 64472, 55045, 18399, 57017, 32236, 60291, 8059, 31457, 63407, 44553, 32553, 10369, 62089, 34830, 64675, 41476, 51745, 8246, 65106, 20738, 58641, 4123, 63813, 17415, 37953, 16492, 38278, 34155, 22615, 22184, 22038, 5546, 24923, 23199, 11019, 2773, 45230, 44368, 44076, 11092, 49846, 46398, 49990, 17725, 13396, 21452, 3349, 5363, 53584, 20271, 34443, 35450, 26792, 42904, 6698, 10726, 41631, 40542, 8056, 30689, 63359, 32265, 32224, 57219, 56825, 63523, 16112, 61378, 61181, 64530, 64448, 48901, 48113, 61509, 21531, 6828, 16811, 43711, 20587, 27312, 1707, 43770, 43062, 13656, 33622, 21885, 41174, 54624, 3414, 22003, 26900, 5015, 37178, 14703, 42063, 20060, 17638, 58812, 53800, 10030, 8819, 29406, 18589, 40120, 35276, 52087, 19064, 30646, 42876, 31577, 10719, 57047, 40430, 60771, 38128, 61292, 20215, 63154, 21438, 48557, 15323, 56005, 54238, 56621, 15827, 53955, 20341, 29873, 63308, 19209, 42939, 47705, 31654, 42373, 40682, 59746, 61079, 38418, 3623, 9970, 57968, 28446, 14492, 39880, 35261, 48247, 7246, 19940, 50399, 56892, 28984, 14223, 4985, 30957, 36015, 44660, 51944, 59190, 12986, 47566, 11165, 40149, 6493, 23783, 38351, 52843, 25972, 29595, 22330, 14761, 64310, 13573, 45905, 20557, 60629, 54292, 52546, 16691, 63083, 27146, 26273, 41114, 55721, 43047, 39555, 33382, 40189, 64612, 53191, 50737, 29682, 61837, 16153, 6337, 14841, 63687, 40845, 35937, 59364, 58137, 32306, 12674, 56331, 2596, 49315, 41536, 28713, 10384, 649, 35070, 47125, 5192, 33093, 17535, 57426, 20768, 1298, 4603, 57758, 40223, 6610, 53735, 34421, 29818, 26440, 18329, 49979, 14909, 13220, 41933, 3305, 59636, 52880, 36658, 23450, 39333, 47515, 39495, 28263, 26258, 58986, 26906, 46900, 13129, 29493, 13453, 56526, 52516, 52435, 53812, 57517, 44064, 2754, 49654, 33457, 45182, 11016, 2005, 49497, 22591, 5508, 33771, 1377, 24827, 22032, 4010, 63520, 7944, 33265, 61567, 57469, 31776, 1986, 49657, 61503, 15888, 993, 57597, 49401, 63552, 3972, 33777, 28302, 36242, 59610, 55576, 47671, 13894, 41829, 25693, 56604, 6947, 53683, 45615, 29805, 27788, 18121, 51386, 1678, 36346, 26848, 57240, 6712, 14310, 41855, 32349, 3356, 7155, 53696, 48943, 13424, 28620, 18173, 64698, 30392, 46986, 27513, 30869, 56031, 56870, 44515, 57939, 60784, 28435, 55026, 61738, 46525, 48203, 23493, 50341, 35931, 23156, 50600, 42811, 12650, 27087, 5789, 40170, 6325, 46312, 35663, 20085, 25300, 54174, 11578, 14803, 57748, 37663, 6450, 12775, 34381, 19578, 25800, 51100, 49959, 9789, 12900, 25550, 3225, 39156, 51600, 36663, 18619, 47800, 35756, 43893, 8939, 60126, 11950, 44498, 37238, 30063, 5975, 22249, 17878, 54715, 23900, 23459, 8137, 51425, 64655, 36356, 32548, 9089, 62009, 14350, 16274, 37313, 63773, 7175, 65096, 18178, 58481, 28700, 24372, 13217, 62267, 14861, 31951, 52868, 52457, 59444, 48744, 26434, 58997, 29722, 63902, 40199, 39377, 53351, 18131, 53946, 27948, 11155, 6987, 19173, 46255, 44620, 36262, 42355, 55896, 22310, 13974, 38346, 26973, 23703, 23004, 56231, 40379, 47715, 26479, 28313, 30442, 59786, 46008, 46925, 15221, 29893, 52958, 56626, 60884, 54035, 7336, 42980, 51839, 32310, 29344, 40846, 10745, 63703, 14672, 20423, 38141, 64620, 58688, 16155, 21490, 61869, 8673, 57567, 7694, 3554, 34692, 33657, 30776, 14216, 17346, 49597, 15388, 7108, 3847, 1777, 61552, 28432, 61397, 54289, 64834, 16643, 48977, 20545, 62725, 1035, 57257, 43041, 64131, 33286, 32417, 41090, 59913, 2070, 31613, 31877, 47049, 51273, 60915, 61971, 57122, 8481, 63226, 63754, 28561, 37009, 56293, 58405, 48707, 16962, 4224, 32752, 2047, 65273, 16896, 65471, 8188, 64481, 8448, 65504, 4094, 65009, 33792, 65405, 16376, 63425, 38456, 14186, 25463, 30365, 22750, 56744, 36315, 55923, 11375, 28372, 50926, 60730, 45500, 47951, 7093, 46309, 4902, 9709, 12895, 24270, 19608, 38836, 51580, 31543, 9804, 19418, 25790, 48540, 39216, 12135, 37623, 63086, 35039, 56952, 36328, 59251, 9082, 31197, 14238, 40393, 4541, 48367, 7119, 52965, 18164, 62394, 28476, 15249, 22923, 35495, 39083, 43624, 26155, 10906, 25258, 43422, 45846, 5453, 12629, 21711, 52310, 21812, 50516, 21307, 17783, 30395, 22380, 27561, 5595, 56043, 23983, 44707, 35566, 60790, 44760, 55122, 11190, 46549, 47966, 23877, 5181, 15596, 17359, 52925, 20724, 62384, 3899, 15089, 10362, 31192, 34718, 40313, 41448, 59231, 7798, 30178, 58411, 10780, 17058, 41406, 37033, 43120, 2695, 34550, 51285, 21560, 34116, 17275, 8529, 20703, 5390, 3563, 24348, 7073, 61883, 47631, 31855, 28292, 50921, 59450, 48696, 14146, 58229, 29725, 63710, 56584, 36305, 53363
};

static const uint32_t iGMb[1024] = {
	1, 65281, 61441, 65521, 49153, 65473, 64513, 65533, 32769, 65409, 63489, 65529, 57345, 65505, 65025, 65535, 2040, 2056, 32896, 32897, 510, 514, 8224, 57377, 1020, 1028, 16448, 49217, 255, 257, 4112, 61457, 46359, 59830, 39762, 44700, 27974, 47726, 42709, 11175, 55948, 29915, 19881, 22350, 13987, 23863, 54123, 38356, 2469, 23306, 45211, 26033, 49770, 38595, 27687, 55661, 34003, 11653, 55374, 45785, 24885, 52066, 46612, 60599, 7325, 25373, 12746, 13874, 50984, 55496, 35955, 36237, 36431, 45455, 6373, 6937, 25492, 27748, 50746, 50887, 564, 52227, 49188, 56513, 141, 29441, 12297, 63281, 282, 58882, 24594, 61025, 32839, 47489, 38917, 64409, 32478, 8831, 10222, 4648, 40888, 18592, 35324, 1162, 16239, 37184, 5111, 2324, 20444, 9296, 17662, 581, 62750, 58102, 12114, 44592, 48456, 47294, 35797, 11148, 31375, 29051, 6057, 22296, 24228, 23647, 50667, 5574, 49484, 46274, 19477, 60237, 12371, 44337, 54022, 64212, 24742, 23137, 42507, 62887, 38954, 54937, 27011, 32106, 20380, 25680, 17658, 1605, 5095, 6420, 37183, 49554, 10190, 12840, 8829, 33571, 35316, 3210, 51360, 24777, 37145, 59282, 30994, 61050, 58439, 47589, 40517, 48031, 51341, 29641, 15497, 30525, 61988, 56563, 53027, 56784, 15028, 19515, 50092, 21700, 3757, 21263, 12523, 5425, 7514, 42526, 25046, 10850, 34647, 43400, 39030, 35481, 50690, 65223, 60513, 40941, 45441, 32690, 64281, 59388, 25345, 65380, 63025, 53239, 55489, 16345, 64909, 29694, 55751, 14810, 40349, 25502, 30322, 36471, 59240, 39144, 60644, 7405, 52943, 12751, 15161, 51004, 29620, 19572, 43038, 58025, 10882, 32299, 43528, 63659, 35489, 24459, 21519, 61781, 5441, 48918, 21764, 64598, 50513, 44998, 43477, 11178, 47774, 25275, 60022, 35563, 44712, 22703, 54507, 5589, 23887, 45406, 30011, 50550, 22356, 44120, 27119, 4458, 5791, 24855, 23164, 33883, 17832, 22598, 46328, 2229, 35664, 45196, 11582, 49710, 8916, 11299, 9532, 50214, 16980, 44099, 2383, 45322, 4245, 27409, 4766, 25107, 8490, 54818, 33960, 22661, 34891, 46473, 13450, 30261, 25417, 46948, 36131, 56718, 55507, 11737, 6725, 47899, 45477, 23474, 50834, 28359, 60522, 38637, 43534, 62123, 10913, 24363, 43652, 31915, 51881, 22475, 21767, 63830, 38225, 44950, 21826, 48726, 58709, 44006, 4028, 17424, 16636, 1089, 1007, 4356, 4159, 49425, 2014, 8712, 8318, 33313, 33272, 2178, 34848, 57481, 24995, 23906, 54811, 58839, 22633, 38745, 30087, 31094, 45266, 11953, 60174, 62188, 44085, 52141, 47812, 15547, 19139, 15691, 54445, 21461, 21169, 20307, 62764, 54518, 42338, 40614, 59991, 43499, 43353, 42922, 31382, 27259, 49045, 27584, 48122, 1724, 61414, 6896, 44799, 431, 57291, 13792, 24061, 862, 30707, 3448, 55168, 32984, 20984, 2130, 34080, 57478, 5246, 33301, 8520, 47138, 10492, 1065, 17040, 28739, 2623, 49419, 4260, 23569, 11699, 19758, 53980, 9427, 19309, 37708, 13495, 18741, 38618, 9879, 26990, 37482, 42423, 18854, 39516, 42139, 31565, 45948, 14261, 19256, 57044, 11487, 52718, 4814, 48551, 22974, 39899, 9628, 28522, 38512, 26359, 2407, 35266, 16010, 59549, 25577, 41585, 36771, 64040, 55547, 17633, 8005, 62543, 45557, 53561, 51154, 32020, 60542, 23535, 4444, 5567, 16662, 22268, 1111, 17776, 36934, 44536, 2222, 35552, 8331, 11134, 33324, 8888, 18467, 38316, 21654, 18779, 42314, 9579, 38182, 21079, 43347, 19158, 10827, 42158, 21157, 37558, 19091, 43308, 54442, 64626, 36605, 61384, 14576, 48925, 58304, 15346, 3644, 32313, 51071, 30692, 7288, 57231, 29152, 7673, 1822, 42133, 27557, 47690, 46779, 59686, 56042, 44691, 28079, 53835, 46547, 23845, 56158, 29843, 28021, 55114, 46808, 12174, 29232, 8953, 1827, 35812, 7308, 51391, 16841, 6087, 14616, 37245, 33682, 17906, 3654, 58464, 41189, 61974, 60147, 44834, 57008, 48262, 31421, 43977, 14252, 30987, 62842, 22417, 28504, 24131, 48479, 54757, 7126, 35359, 57739, 6306, 24089, 25224, 30819, 34345, 55175, 50448, 61638, 3153, 44813, 12612, 48178, 49941, 60356, 41660, 17571, 18988, 54347, 10415, 20777, 4747, 29971, 20830, 41554, 9494, 59942, 37976, 43157, 35142, 47754, 44230, 15021, 43725, 13227, 43826, 52908, 60084, 19691, 22115, 40279, 54631, 39382, 21913, 26454, 30042, 42614, 50288, 37061, 3143, 47373, 12572, 58418, 17170, 60996, 25144, 51299, 34340, 56455, 6286, 29209, 8585, 30498, 2451, 27914, 53402, 26321, 16997, 39747, 46119, 55733, 33994, 13957, 26701, 45929, 41267, 52642, 55828, 60635, 19228, 58444, 17586, 20037, 4807, 14611, 37165, 54162, 9614, 29222, 8793, 42787, 35172, 40074, 51351, 27081, 2112, 49161, 132, 31745, 528, 61443, 33, 57089, 1056, 57349, 66, 48641, 264, 63490, 32785, 61313, 48575, 16830, 7132, 9244, 28528, 36976, 1783, 2311, 57056, 8415, 3566, 4622, 14264, 18488, 33660, 33924, 63467, 5624, 24447, 33120, 32251, 1406, 22496, 8280, 64502, 2812, 44992, 16560, 48894, 703, 11248, 4140, 37105, 3985, 63760, 61690, 58429, 50149, 15940, 48191, 51321, 34761, 31880, 30845, 61983, 57843, 7970, 56864, 3668, 44047, 49382, 6849, 917, 27396, 45114, 50865, 1834, 54792, 24691, 36193, 33227, 13698, 22557, 58201, 11502, 4653, 8911, 12579, 35644, 50316, 18612, 19529, 5751, 35095, 37224, 39058, 17822, 25158, 9306, 42533, 41834, 38564, 27191, 51563, 43227, 9641, 23182, 29275, 20917, 19282, 46364, 58550, 54382, 37589, 11591, 47406, 12186, 26160, 25338, 1635, 35815, 6540, 39103, 16793, 6093, 13080, 12669, 33586, 50676, 3270, 52320, 41165, 36837, 7056, 47359, 441, 58362, 1764, 28224, 49263, 51187, 3528, 56448, 32989, 29181, 882, 14112, 57400, 42078, 41637, 10822, 47659, 43288, 59562, 35474, 28299, 21039, 53587, 5411, 56598, 21644, 29781, 17737, 46918, 28874, 13937, 26381, 62312, 39987, 52637, 55748, 15578, 14437, 39737, 45959, 31156, 52762, 59087, 27874, 7789, 50734, 53959, 11363, 40237, 45452, 29874, 19225, 59212, 25367, 59748, 38450, 52887, 22726, 14937, 42381, 29606, 15196, 42044, 17334, 19012, 3799, 10511, 37102, 4753, 7598, 21022, 8667, 9506, 34668, 38024, 18551, 35145, 839, 47364, 36917, 52113, 16594, 11841, 58382, 62181, 33188, 23682, 51227, 58825, 8297, 38689, 29191, 63859, 14151, 47416, 37749, 35732, 19922, 11854, 58590, 8933, 39844, 23708, 51643, 17866, 9961, 5927, 29295, 37235, 31760, 61565, 1985, 16136, 7940, 64544, 49649, 4034, 15880, 63551, 33761, 8068, 3970, 32272, 57593, 2017, 61527, 43505, 40710, 64160, 31766, 60029, 42946, 16040, 63532, 54521, 20355, 32080, 15883, 62783, 21473, 8020, 11725, 13102, 13021, 9011, 52084, 36044, 52408, 18637, 38631, 6551, 39279, 37274, 26042, 18022, 26204, 42087, 28879, 12657, 5901, 62232, 23604, 52317, 50628, 15558, 47208, 39097, 35719, 31116, 11802, 58927, 25314, 7779, 60934, 64239, 44769, 8111, 48002, 32444, 60345, 18412, 30467, 64888, 55153, 36824, 24001, 16222, 62941, 9206, 52863, 33231, 7400, 6173, 29600, 24692, 1850, 50696, 59200, 49384, 3700, 35855, 14800, 12346, 925, 25348, 32155, 25982, 22490, 9816, 24423, 39264, 38391, 2454, 48846, 12991, 11245, 4908, 44980, 19632, 51964, 1227, 50776, 43207, 35942, 39565, 12694, 27186, 41754, 59044, 25388, 54372, 17971, 52551, 6347, 13593, 20877, 29522, 34580, 60552, 51314, 36553, 8645, 15138, 45597, 58291, 17290, 30276, 25657, 51045, 37091, 7569, 55567, 61914	
};

static uint32_t
mf_conv_small(int x)
{
	/*
	 * If x < 0, the cast to uint32_t will set the high bit to 1.
	 */
	uint32_t y;

	y = (uint32_t)x;
	y += FERMAT_P & -(y >> 31);
	return y;
}

/* see inner.h */
static int32_t
mf_conv_signed(uint32_t x)
{
	/*
	 * If x > P/2, subtract P.
	 */
	int32_t y;

	y = (int32_t)x;
	y -= (int32_t)(FERMAT_P & -(((FERMAT_P >> 1) - x) >> 31));
	return y;
}

/* see inner.h */
static uint32_t
mf_add(uint32_t x, uint32_t y)
{
	/*
	 * We compute x + y - q. If the result is negative, then the
	 * high bit will be set, and 'd >> 31' will be equal to 1;
	 * thus '-(d >> 31)' will be an all-one pattern. Otherwise,
	 * it will be an all-zero pattern. In other words, this
	 * implements a conditional addition of q.
	 */
	uint32_t d;

	d = x + y - FERMAT_P;
	d += FERMAT_P & -(d >> 31);
	return d;
}

/* see inner.h */
static uint32_t
mf_sub(uint32_t x, uint32_t y)
{
	/*
	 * As in mq_add(), we use a conditional addition to ensure the
	 * result is in the 0..q-1 range.
	 */
	uint32_t d;

	d = x - y;
	d += FERMAT_P & -(d >> 31);
	return d;
}

/*
 * Division by 2 modulo q. Operand must be in the 0..q-1 range.
 */
static inline uint32_t
mf_rshift1(uint32_t x)
{
	x += FERMAT_P & -(x & 1);
	return (x >> 1);
}

static uint32_t
mf_mul(uint32_t x, uint32_t y)
{
	uint32_t z;

	z = x * y;
	z = (z & 0xFFFF) - (z >> 16);
	z += FERMAT_P & -(z >> 31);

	/*
	 * Warning: x = y = 2^16 will overflow here, giving (incorrectly) z = 0.
	 * Instead, we should have z = 1 in this case.
	 */
	z += (x & y) >> 16;

	return z;
}

#define mf_sqr(x) mf_mul((x), (x))

static uint32_t
mf_div(uint32_t x, uint32_t y)
{
	/*
	 * We invert y by computing y^(p-2) mod p.
	 * We use the following addition chain, which provides exponents (binary):
	 *
	 * e0 = 1
	 * e1 = 10
	 * e2 = 11
	 * e3 = 1100
	 * e4 = 1111
	 * e5 = 11110000
	 * e6 = 11111111
	 * e7 = 1111111100000000
	 * e8 = 1111111111111111 = p-2
	 */
	uint32_t y1, y2, y3, y4, y5, y6, y7, y8;

	y1 = mf_sqr(y);
	y2 = mf_mul(y1, y);
	y3 = mf_sqr(mf_sqr(y2));
	y4 = mf_mul(y3, y2);
	y5 = mf_sqr(mf_sqr(mf_sqr(mf_sqr(y4))));
	y6 = mf_mul(y5, y4);
	y7 = mf_sqr(mf_sqr(mf_sqr(mf_sqr(mf_sqr(mf_sqr(mf_sqr(mf_sqr(y6))))))));
	y8 = mf_mul(y7, y6);

	return mf_mul(x, y8);
}

/* see inner.h */
static void
mf_NTT(uint32_t *a, unsigned logn)
{
	size_t n, t, m;

	n = MKN(logn);
	t = n;
	for (m = 1; m < n; m <<= 1) {
		size_t ht, i, j1;

		ht = t >> 1;
		for (i = 0, j1 = 0; i < m; i ++, j1 += t) {
			size_t j, j2;
			uint32_t s;

			s = GMb[m + i];
			j2 = j1 + ht;
			for (j = j1; j < j2; j ++) {
				uint32_t u, v;

				u = a[j];
				v = mf_mul(a[j + ht], s);
				a[j] = (uint32_t)mf_add(u, v);
				a[j + ht] = (uint32_t)mf_sub(u, v);
			}
		}
		t = ht;
	}
}

/* see inner.h */
static void
mf_iNTT(uint32_t *a, unsigned logn)
{
	size_t n, t, m;
	uint32_t ni;

	n = MKN(logn);
	t = 1;
	m = n;
	while (m > 1) {
		size_t hm, dt, i, j1;

		hm = m >> 1;
		dt = t << 1;
		for (i = 0, j1 = 0; i < hm; i ++, j1 += dt) {
			size_t j, j2;
			uint32_t s;

			j2 = j1 + t;
			s = iGMb[hm + i];
			for (j = j1; j < j2; j ++) {
				uint32_t u, v, w;

				u = a[j];
				v = a[j + t];
				a[j] = (uint32_t)mf_add(u, v);
				w = mf_sub(u, v);
				a[j + t] = (uint32_t) mf_mul(w, s);
			}
		}
		t = dt;
		m = hm;
	}

	/*
	 * To complete the inverse NTT, we must now divide all values by
	 * n (the vector size). We thus need the inverse of n, i.e. we
	 * need to divide 1 by 2 logn times. But we also want it in
	 * Montgomery representation, i.e. we also want to multiply it
	 * by R = 2^16. In the common case, this should be a simple right
	 * shift. The loop below is generic and works also in corner cases;
	 * its computation time is negligible.
	 */
	ni = 1;
	for (m = n; m > 1; m >>= 1) {
		ni = mf_rshift1(ni);
	}
	for (m = 0; m < n; m ++) {
		a[m] = (uint32_t)mf_mul(a[m], ni);
	}
}

/* see inner.h */
static void
mf_int8_to_NTT(uint32_t *restrict p, const int8_t *restrict f, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u++) {
		p[u] = mf_conv_small(f[u]);
	}
	mf_NTT(p, logn);
}

/* see inner.h */
static void
mf_poly_div(uint32_t *f, uint32_t *g, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u ++) {
		f[u] = (uint32_t)mf_div(f[u], g[u]);
	}
}

/* see inner.h */
static int
mf_is_invertible(int8_t *f, unsigned logn, uint8_t *restrict tmp)
{
	uint32_t *p;
	size_t n, u;
	int res;

	n = MKN(logn);
	p = (uint32_t*)tmp;
	res = 1;
	mf_int8_to_NTT(p, f, logn);
	for (u = 0; u < n; u++) {
		// res &= (p[u] > 0)
		res &= 1U - ((p[u] - 1U) >> 15);
	}
	return res;
}

/* see inner.h */
static void
mf_NTRU(
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	uint32_t *restrict bf, uint32_t *restrict bg,
	uint32_t *restrict bF, uint32_t *restrict bG, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);

	mf_int8_to_NTT(bf, f, logn);
	mf_int8_to_NTT(bg, g, logn);
	mf_int8_to_NTT(bF, F, logn);

	if (G == NULL) {
		/*
		 * Compute (in NTT representation) G = (1 + gF) / f.
		 */
		for (u = 0; u < n; u++) {
			bG[u] = mf_mul(bg[u], bF[u]);
			bG[u] = mf_add(1, bG[u]);
		}
		mf_poly_div(bG, bf, logn);
	} else {
		mf_int8_to_NTT(bG, G, logn);
	}
}

/******************************************************************************
 * Now functions that use the above
 *****************************************************************************/

static void
hash_to_ntt(uint32_t *p, const uint8_t *h, unsigned logn)
{
	size_t n, u, v;
	uint8_t hash;

	n = MKN(logn);

	if (logn <= 3) {
		for (u = 0; u < n; u ++) {
			p[u] = (h[0] >> u) & 1;
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				p[u] = (hash >> v) & 1;
			}
		}
	}

	mf_NTT(p, logn);
}

static int
sample_short_NTT(prng *rng, uint32_t *restrict x0, uint32_t *restrict x1,
	const uint32_t *restrict bf, const uint32_t *restrict bg,
	const uint32_t *restrict bF, const uint32_t *restrict bG,
	const uint8_t *restrict h, unsigned logn)
{
	size_t n, u;
	int32_t norm, center, offset, latcoord;

	n = MKN(logn);
	norm = 0;

	hash_to_ntt(x0, h, logn);
	hash_to_ntt(x1, SECOND_HASH(h, logn), logn);

	/*
	 * Set the target vector to (t0, t1) = B * (h0, h1), i.e.:
	 *     t0 = f h0 + F h1,
	 *     t1 = g h0 + G h1.
	 */

	for (u = 0; u < n; u ++) {
		uint32_t res0, res1;

		res0 = mf_add(mf_mul(bf[u], x0[u]), mf_mul(bF[u], x1[u]));
		res1 = mf_add(mf_mul(bg[u], x0[u]), mf_mul(bG[u], x1[u]));
		x0[u] = res0;
		x1[u] = res1;
	}

	/*
	 * Sample and write the result in (x0, x1). Gaussian smoothing is used to
	 * not reveal information on the secret basis.
	 */
	mf_iNTT(x0, logn);
	mf_iNTT(x1, logn);

	/*
	 * For the NTT, we store the lattice point B*s close to target t = B*h/2 in
	 * memory, so a signature is given by calculating B^{-1} (B*s). In the FFT
	 * case, we store the doubled distance from the lattice point to target.
	 * Here it is convenient to store the lattice point straight on, as the
	 * division by two works best as soon as possible.
	 */
	for (u = 0; u < n; u ++) {
		center = mf_conv_signed(x0[u]);
		offset = mkgauss_sign(rng, center & 1, logn);
		latcoord = (center - offset) / 2;
		x0[u] = mf_conv_small(latcoord);
		norm += offset * offset;
	}
	for (u = 0; u < n; u ++) {
		center = mf_conv_signed(x1[u]);
		offset = mkgauss_sign(rng, center & 1, logn);
		latcoord = (center - offset) / 2;
		x1[u] = mf_conv_small(latcoord);
		norm += offset * offset;
	}

	mf_NTT(x0, logn);
	mf_NTT(x1, logn);

	/*
	 * Test whether the l2-norm of (x0, x1) is below the given bound. The
	 * code below uses only 32-bit operations to compute the squared norm,
	 * since the max. value is 2n * 128^2 <= 2^24 (when logn <= 9).
	 * For a large enough verification margin, it is unlikely that the
	 * norm of the gaussian (x0, x1) is too large.
	 */
	return (uint32_t)norm <= L2BOUND(logn);
}

/* see inner.h */
int
Zf(uncompressed_sign_NTT)(inner_shake256_context *rng,
	int16_t *restrict s0, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u;
	uint32_t flag, *bf, *bg, *bF, *bG, *x0, *x1;
	int norm_okay;
	prng p;

	n = MKN(logn);

	bf = (uint32_t *)tmp;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	x0 = bG + n;
	x1 = x0 + n;

	Zf(prng_init)(&p, rng);

	mf_NTRU(f, g, F, G, bf, bg, bF, bG, logn);

	norm_okay = sample_short_NTT(&p, x0, x1, bf, bg, bF, bG, h, logn);

	/*
	 * Compute (s0, s1) = ((h0, h1) - B^{-1} (x0, x1)) / 2, so
	 *
	 *     s0 = (h0 - (x0 * G + x1 (-F))) / 2,
	 *     s1 = (h1 - (x0 * (-g) + x1 f)) / 2.
	 */

	for (u = 0; u < n; u++) {
		uint32_t z0, z1;
		z0 = mf_sub(mf_mul(bG[u], x0[u]), mf_mul(bF[u], x1[u]));
		z1 = mf_sub(mf_mul(bf[u], x1[u]), mf_mul(bg[u], x0[u]));

		x0[u] = z0;
		x1[u] = z1;
	}

	mf_iNTT(x0, logn);
	mf_iNTT(x1, logn);

	/*
	 * The polynomial x0 is stored at s0, so conversion to int16_t is done
	 * automatically. Normalize s0 elements into the [-q/2..q/2] range.
	 * Do similarly for s1/x1.
	 */
	for (u = n; u -- > 0; ) {
		s0[u] = mf_conv_signed(x0[u]);
		s1[u] = mf_conv_signed(x1[u]);
	}

	flag = (uint32_t)Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn);
	conditional_flip(flag, s0, h, logn);
	conditional_flip(flag, s1, SECOND_HASH(h, logn), logn);
	return norm_okay;
}

/* see inner.h */
int
Zf(sign_NTT)(inner_shake256_context *rng, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u;
	uint32_t flag, *bf, *bg, *bF, *bG, *x0, *x1;
	int norm_okay;
	prng p;

	n = MKN(logn);

	bf = (uint32_t *)tmp;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	x0 = bG + n;
	x1 = x0 + n;

	Zf(prng_init)(&p, rng);

	mf_NTRU(f, g, F, G, bf, bg, bF, bG, logn);

	norm_okay = sample_short_NTT(&p, x0, x1, bf, bg, bF, bG, h, logn);

	/*
	 * Compute s1 in (s0, s1) = ((h0, h1) - B^{-1} (x0, x1)) / 2, so
	 *
	 *     s1 = (h1 - (x0 * (-g) + x1 f)) / 2.
	 */
	for (u = 0; u < n; u++) {
		x1[u] = mf_sub(mf_mul(bf[u], x1[u]), mf_mul(bg[u], x0[u]));
	}
	mf_iNTT(x1, logn);

	/*
	 * The polynomial x1 is stored at s1, so conversion to int16_t is done
	 * automatically. Normalize s1 elements into the [-q/2..q/2] range.
	 */
	for (u = 0; u < n; u ++) {
		s1[u] = mf_conv_signed(x1[u]);
	}

	flag = (uint32_t)Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn);
	conditional_flip(flag, s1, SECOND_HASH(h, logn), logn);

	return norm_okay;
}

/* =================================================================== */

/* see inner.h */
void
Zf(expand_seckey)(fpr *restrict expanded_seckey,
	const int8_t *f, const int8_t *g, const int8_t *F, unsigned logn)
{
	size_t n;
	fpr *bf, *bg, *bF, *bG;

	n = MKN(logn);
	bf = expanded_seckey;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;

	/*
	 * We load the secret key elements directly into the 2x2 matrix B.
	 */
	construct_basis(f, g, F, NULL, bf, bg, bF, bG, logn);

#if HAWK_RECOVER_CHECK
	Zf(poly_invnorm2_fft)(bG + n, bf, bg, logn);
#endif
}
