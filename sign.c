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
static const uint16_t gauss_hi[10] = {
	0x580B, 0x35F9,
	0x1D34, 0x0DD7,
	0x05B7, 0x020C,
	0x00A2, 0x002B,
	0x000A, 0x0001,
};

static const uint64_t gauss_lo[26] = {
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

static inline int8_t
mkgauss_sign(prng *rng, uint8_t parity)
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
		p_lo = (gauss_lo[k] & (parity - 1u)) | ((gauss_lo[k + 1] & -parity));

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
		 *     p_lo = gauss_lo[k + parity];
		 *     p_hi = gauss_hi[k + parity];
		 */
		p_lo = (gauss_lo[k] & (parity - 1u)) | ((gauss_lo[k + 1] & -parity));
		p_hi = (gauss_hi[k] & (parity - 1u)) | ((gauss_hi[k + 1] & -parity));

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
				p[u] = fpr_of(hash & 1);
				hash >>= 1;
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
		 * Compute G = (1 + gF) / f, where all polynomials are in FFT
		 * representation.
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
		z = mkgauss_sign(rng, z);
		x0[u] = fpr_of(z);
		norm += z*z;
	}
	for (u = 0; u < n; u ++) {
		z = fpr_rint(x1[u]) & 1;
		z = mkgauss_sign(rng, z);
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
	return (uint32_t)norm <= Zf(l2bound)[logn];
}

/*
 * The following are helper functions for the sign-functions. If a lattice
 * point is generated that is too far away from (h0, h1) / 2, s0 and s1 are
 * untouched and 0 is returned; the caller should then try again. Otherwise, 1
 * is returned and (s0, s1) contain a valid signature for (h0, h1).
 */

/* helper for Zf(sign_dyn) */
static int
do_sign_dyn(prng *rng, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u;
	fpr *x0, *x1, *bf, *bg, *bF, *bG;
	uint16_t flag;

	n = MKN(logn);
	bf = (fpr *)tmp;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	x0 = bG + n;
	x1 = x0 + n;

	construct_basis(f, g, F, G, bf, bg, bF, bG, logn);

	if (!sample_short(rng, x0, x1, bf, bg, bF, bG, h, logn)) {
		return 0;
	}

	/*
	 * Compute *twice* the rounding error, which is given by:
	 *
	 *     (f* x0 + g* x1) / (f* f + g* g).
	 *
	 * If the above quantity is in the (-1, 1)^n box, simple rounding works for
	 * recovering s0. Otherwise, reject this signature.
	 */
	Zf(poly_add_muladj_fft)(bF, x0, x1, bf, bg, logn);
	Zf(poly_invnorm2_fft)(bG, bf, bg, logn);
	Zf(poly_mul_autoadj_fft)(bF, bG, logn);
	Zf(iFFT)(bF, logn);

	for (u = 0; u < n; u++) {
		if (!fpr_lt(fpr_neg(fpr_one), bF[u]) || !fpr_lt(bF[u], fpr_one)) {
			return 0;
		}
	}

	/*
	 * Compute s1 in (s0, s1) = ((h0, h1) - B^{-1} (x0, x1)) / 2, so
	 *
	 *     s1 = (h1 - (x0 * (-g) + x1 f)) / 2.
	 */
	Zf(poly_neg)(bg, logn);
	Zf(poly_add_mul_fft)(bF, x0, x1, bg, bf, logn);
	Zf(iFFT)(bF, logn);
	noise_to_lattice(s1, SECOND_HASH(h, logn), bF, logn);

	flag = (uint16_t)Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn);
	conditional_flip(flag, s1, SECOND_HASH(h, logn), logn);

	return 1;
}

/* helper for Zf(sign) */
static int
do_sign(prng *rng, int16_t *restrict s1,
	const fpr *restrict expanded_seckey, const uint8_t *restrict h,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u;
	const fpr *bf, *bg, *bF, *bG, *invq00;
	fpr *x0, *x1, *res;
	uint16_t flag;

	n = MKN(logn);

	bf = expanded_seckey;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	invq00 = bG + n;

	x0 = (fpr *)tmp;
	x1 = x0 + n;
	res = x1 + n;

	if (!sample_short(rng, x0, x1, bf, bg, bF, bG, h, logn)) {
		return 0;
	}

	/*
	 * Compute *twice* the rounding error, which is given by:
	 *
	 *     (f* x0 + g* x1) / (f* f + g* g).
	 *
	 * If the above quantity is in the (-1, 1)^n box, simple rounding works for
	 * recovering s0. Otherwise, reject this signature.
	 */
	Zf(poly_add_muladj_fft)(res, x0, x1, bf, bg, logn);
	Zf(poly_mul_autoadj_fft)(res, invq00, logn);
	Zf(iFFT)(res, logn);

	for (u = 0; u < n; u++) {
		if (!fpr_lt(fpr_neg(fpr_one), res[u]) || !fpr_lt(res[u], fpr_one)) {
			return 0;
		}
	}

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

	return 1;
}


static inline void
int8_to_ntt(uint16_t *restrict p, const int8_t *restrict f, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u++) {
		p[u] = Zf(mq_conv_small)(f[u]);
	}
	Zf(mq_NTT)(p, logn);
}

/* helper for Zf(sign_simple) */
static int
do_sign_simple(prng *rng,
	int16_t *restrict s0, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u, v, w;
	uint8_t h0, h1;
	uint16_t flag, *bf, *bg, *bF, *bG, *x0, *x1;
	int32_t norm, z;

	n = MKN(logn);
	norm = 0;

	bf = (uint16_t *)tmp;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	x0 = (uint16_t *)s0;
	x1 = (uint16_t *)s1;

	int8_to_ntt(bf, f, logn);
	int8_to_ntt(bg, g, logn);
	int8_to_ntt(bF, F, logn);

	if (G == NULL) {
		/*
		 * Compute G = (1 + gF) / f, where all polynomials are in NTT
		 * representation.
		 */
		for (u = 0; u < n; u++) {
			bG[u] = Zf(mq_mul)(bg[u], bF[u]);
			bG[u] = Zf(mq_add)(1, bG[u]);
		}
		Zf(mq_poly_div)(bG, bf, logn);
	} else {
		int8_to_ntt(bG, G, logn);
	}

	/*
	 * Put the hash (h0, h1) inside x0, x1.
	 */
	if (logn <= 3) {
		for (u = 0; u < n; u ++) {
			x0[u] = (h[0] >> u) & 1;
			x1[u] = (h[1] >> u) & 1;
		}
	} else {
		for (u = w = 0; u < n; w ++) {
			h0 = h[w];
			h1 = h[n/8 + w];
			for (v = 0; v < 8; v ++, u ++) {
				x0[u] = h0 & 1;
				x1[u] = h1 & 1;
				h0 >>= 1;
				h1 >>= 1;
			}
		}
	}

	Zf(mq_NTT)(x0, logn);
	Zf(mq_NTT)(x1, logn);

	Zf(mq_poly_tomonty)(x0, logn);
	Zf(mq_poly_tomonty)(x1, logn);

	/*
	 * Set the target vector to (t0, t1) = B * (h0, h1), i.e.:
	 *     t0 = f h0 + F h1,
	 *     t1 = g h0 + G h1.
	 */
	for (u = 0; u < n; u ++) {
		uint32_t res0, res1;

		res0 = Zf(mq_add)(Zf(mq_montymul)(bf[u], x0[u]),
			Zf(mq_montymul)(bF[u], x1[u]));
		res1 = Zf(mq_add)(Zf(mq_montymul)(bg[u], x0[u]),
			Zf(mq_montymul)(bG[u], x1[u]));
		x0[u] = res0;
		x1[u] = res1;
	}

	/*
	 * Sample and write the result in (x0, x1). Gaussian smoothing is used to
	 * not reveal information on the secret basis.
	 */
	Zf(mq_iNTT)(x0, logn);
	Zf(mq_iNTT)(x1, logn);

	for (u = 0; u < n; u ++) {
		z = Zf(mq_conv_signed)(x0[u]);
		z = mkgauss_sign(rng, z & 1);
		x0[u] = Zf(mq_conv_small)(z);
		norm += z*z;
	}
	for (u = 0; u < n; u ++) {
		z = Zf(mq_conv_signed)(x1[u]);
		z = mkgauss_sign(rng, z & 1);
		x1[u] = Zf(mq_conv_small)(z);
		norm += z*z;
	}

	/*
	 * Test whether the l2-norm of (x0, x1) is below the given bound. The
	 * code below uses only 32-bit operations to compute the squared norm,
	 * since the max. value is 2n * 128^2 <= 2^24 (when logn <= 9).
	 * For a large enough verification margin, it is unlikely that the
	 * norm of the gaussian (x0, x1) is too large.
	 */
	if ((uint32_t)norm > Zf(l2bound)[logn]) {
		return 0;
	}

	/*
	 * Norm of (x0, x1) is acceptable.
	 */
	Zf(mq_NTT)(x0, logn);
	Zf(mq_NTT)(x1, logn);

	Zf(mq_poly_tomonty)(x0, logn);
	Zf(mq_poly_tomonty)(x1, logn);

	/*
	 * Compute (s0, s1) = ((h0, h1) - B^{-1} (x0, x1)) / 2, so
	 *
	 *     s0 = (h0 - (x0 * G + x1 (-F))) / 2,
	 *     s1 = (h1 - (x0 * (-g) + x1 f)) / 2.
	 */

	for (u = 0; u < n; u++) {
		uint16_t z0, z1;
		z0 = Zf(mq_sub)(Zf(mq_montymul(bG[u], x0[u])),
			Zf(mq_montymul)(bF[u], x1[u]));
		z1 = Zf(mq_sub)(Zf(mq_montymul)(bf[u], x1[u]),
			Zf(mq_montymul)(bg[u], x0[u]));

		x0[u] = z0;
		x1[u] = z1;
	}
	Zf(mq_iNTT)(x0, logn);
	Zf(mq_iNTT)(x1, logn);

	/*
	 * The polynomial x0 is stored at s0, so conversion to int16_t is done
	 * automatically. Normalize s0 elements into the [-q/2..q/2] range.
	 * Do similarly for s1/x1.
	 */
	for (u = 0; u < n; u ++) {
		s0[u] = Zf(mq_conv_signed)(x0[u]);
		s1[u] = Zf(mq_conv_signed)(x1[u]);
	}

	n = MKN(logn);
	if (logn <= 3) {
		for (u = 0; u < n; u ++) {
			s0[u] = (((h[0] >> u) & 1) - s0[u]) / 2;
			s1[u] = (((h[1] >> u) & 1) - s1[u]) / 2;
		}
	} else {
		for (u = w = 0; u < n; w ++) {
			h0 = h[w];
			h1 = h[n / 8 + w];
			for (v = 0; v < 8; v ++, u ++) {
				s0[u] = ((h0 & 1) - s0[u]) / 2;
				s1[u] = ((h1 & 1) - s1[u]) / 2;
				h0 >>= 1;
				h1 >>= 1;
			}
		}
	}

	flag = (uint16_t)Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn);
	conditional_flip(flag, s0, h, logn);
	conditional_flip(flag, s1, SECOND_HASH(h, logn), logn);

	return 1;
}

/* helper for Zf(sign_NTT) */
static int
do_sign_NTT(prng *rng, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u, v, w;
	uint8_t h0, h1;
	uint16_t flag, *bf, *bg, *bF, *bG, *x0, *x1;
	int32_t norm, z;

	n = MKN(logn);
	norm = 0;

	bf = (uint16_t *)tmp;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	x0 = bG + n;
	x1 = (uint16_t *)s1;

	int8_to_ntt(bf, f, logn);
	int8_to_ntt(bg, g, logn);
	int8_to_ntt(bF, F, logn);

	if (G == NULL) {
		/*
		 * Compute G = (1 + gF) / f, where all polynomials are in NTT
		 * representation.
		 */
		for (u = 0; u < n; u++) {
			bG[u] = Zf(mq_mul)(bg[u], bF[u]);
			bG[u] = Zf(mq_add)(1, bG[u]);
		}
		Zf(mq_poly_div)(bG, bf, logn);
	} else {
		int8_to_ntt(bG, G, logn);
	}

	/*
	 * Put the hash (h0, h1) inside x0, x1.
	 */
	if (logn <= 3) {
		for (u = 0; u < n; u ++) {
			x0[u] = (h[0] >> u) & 1;
			x1[u] = (h[1] >> u) & 1;
		}
	} else {
		for (u = w = 0; u < n; w ++) {
			h0 = h[w];
			h1 = h[n/8 + w];
			for (v = 0; v < 8; v ++, u ++) {
				x0[u] = h0 & 1;
				x1[u] = h1 & 1;
				h0 >>= 1;
				h1 >>= 1;
			}
		}
	}

	Zf(mq_NTT)(x0, logn);
	Zf(mq_NTT)(x1, logn);

	Zf(mq_poly_tomonty)(x0, logn);
	Zf(mq_poly_tomonty)(x1, logn);

	/*
	 * Set the target vector to (t0, t1) = B * (h0, h1), i.e.:
	 *     t0 = f h0 + F h1,
	 *     t1 = g h0 + G h1.
	 */
	for (u = 0; u < n; u ++) {
		uint32_t res0, res1;

		res0 = Zf(mq_add)(Zf(mq_montymul)(bf[u], x0[u]),
			Zf(mq_montymul)(bF[u], x1[u]));
		res1 = Zf(mq_add)(Zf(mq_montymul)(bg[u], x0[u]),
			Zf(mq_montymul)(bG[u], x1[u]));
		x0[u] = res0;
		x1[u] = res1;
	}

	/*
	 * Sample and write the result in (x0, x1). Gaussian smoothing is used to
	 * not reveal information on the secret basis.
	 */
	Zf(mq_iNTT)(x0, logn);
	Zf(mq_iNTT)(x1, logn);

	for (u = 0; u < n; u ++) {
		z = Zf(mq_conv_signed)(x0[u]);
		z = mkgauss_sign(rng, z & 1);
		x0[u] = Zf(mq_conv_small)(z);
		norm += z*z;
	}
	for (u = 0; u < n; u ++) {
		z = Zf(mq_conv_signed)(x1[u]);
		z = mkgauss_sign(rng, z & 1);
		x1[u] = Zf(mq_conv_small)(z);
		norm += z*z;
	}

	/*
	 * Test whether the l2-norm of (x0, x1) is below the given bound. The
	 * code below uses only 32-bit operations to compute the squared norm,
	 * since the max. value is 2n * 128^2 <= 2^24 (when logn <= 9).
	 * For a large enough verification margin, it is unlikely that the
	 * norm of the gaussian (x0, x1) is too large.
	 */
	if ((uint32_t)norm > Zf(l2bound)[logn]) {
		return 0;
	}

	/*
	 * Norm of (x0, x1) is acceptable.
	 */
	Zf(mq_NTT)(x0, logn);
	Zf(mq_NTT)(x1, logn);

	Zf(mq_poly_tomonty)(x0, logn);
	Zf(mq_poly_tomonty)(x1, logn);

	/*
	 * Compute s1 in (s0, s1) = ((h0, h1) - B^{-1} (x0, x1)) / 2, so
	 *
	 *     s1 = (h1 - (x0 * (-g) + x1 f)) / 2.
	 */

	for (u = 0; u < n; u++) {
		x1[u] = Zf(mq_sub)(Zf(mq_montymul)(bf[u], x1[u]),
			Zf(mq_montymul)(bg[u], x0[u]));
	}
	Zf(mq_iNTT)(x1, logn);

	/*
	 * The polynomial x1 is stored at s1, so conversion to int16_t is done
	 * automatically. Normalize s1 elements into the [-q/2..q/2] range.
	 */
	for (u = 0; u < n; u ++) {
		s1[u] = Zf(mq_conv_signed)(x1[u]);
	}

	n = MKN(logn);
	if (logn <= 3) {
		for (u = 0; u < n; u ++) {
			s1[u] = (((h[1] >> u) & 1) - s1[u]) / 2;
		}
	} else {
		for (u = w = 0; u < n; w ++) {
			h1 = h[n / 8 + w];
			for (v = 0; v < 8; v ++, u ++) {
				s1[u] = ((h1 & 1) - s1[u]) / 2;
				h1 >>= 1;
			}
		}
	}

	flag = (uint16_t)Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn);
	conditional_flip(flag, s1, SECOND_HASH(h, logn), logn);

	return 1;
}

/* =================================================================== */

/* see inner.h */
void
Zf(expand_seckey)(fpr *restrict expanded_seckey,
	const int8_t *f, const int8_t *g, const int8_t *F, unsigned logn)
{
	size_t n;
	fpr *bf, *bg, *bF, *bG, *invq00;

	n = MKN(logn);
	bf = expanded_seckey;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	invq00 = bG + n;

	/*
	 * We load the private key elements directly into the 2x2 matrix B.
	 */
	construct_basis(f, g, F, NULL, bf, bg, bF, bG, logn);
	Zf(poly_invnorm2_fft)(invq00, bf, bg, logn);
}

/*
 * Use a fast PRNG for gaussian sampling during signing.
 */

/* see inner.h */
void
Zf(sign_dyn)(inner_shake256_context *rng, int16_t *restrict sig,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	prng p;
	do {
		Zf(prng_init)(&p, rng);
	} while (!do_sign_dyn(&p, sig, f, g, F, G, h, logn, tmp));
}

/* see inner.h */
void
Zf(sign_NTT)(inner_shake256_context *rng, int16_t *restrict sig,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	prng p;
	do {
		Zf(prng_init)(&p, rng);
	} while (!do_sign_NTT(&p, sig, f, g, F, G, h, logn, tmp));
}

/* see inner.h */
void
Zf(sign_simple)(inner_shake256_context *rng,
	int16_t *restrict s0, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	prng p;
	do {
		Zf(prng_init)(&p, rng);
	} while (!do_sign_simple(&p, s0, s1, f, g, F, G, h, logn, tmp));
}

/* see inner.h */
void
Zf(sign)(inner_shake256_context *rng, int16_t *restrict sig,
	const fpr *restrict expanded_seckey, const uint8_t *restrict h,
	unsigned logn, uint8_t *restrict tmp)
{
	prng p;
	do {
		Zf(prng_init)(&p, rng);
	} while (!do_sign(&p, sig, expanded_seckey, h, logn, tmp));
}

