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

/* To generate the values in the table below, run the following code:

#include<bits/stdc++.h>
int main() {
	long double sigma = 1.292, mu = 0, p63 = powl(2, 63), table[100], csum[100];
	unsigned long long results[2][20] = {};
	for (int i = 0; i < 2; i++, mu += 0.5) {
		for (int x = 100; x --> 0; ) {
			table[x] = expl(-0.5 * (x-mu)*(x-mu) / sigma / sigma);
			csum[x] = table[x];
			if (x < 99) csum[x] += csum[x+1];
		}
		results[i][0] = llroundl((1.0+i) * p63 * table[i] / (csum[i] + csum[1]));
		for (int x = 19; --x >= 1; )
			results[i][x] = llroundl(p63 * csum[1+i+x] / csum[1+i]);
	}
	for (int x = 0; x < 20; x++)
		printf("\t%19lluu, %19lluu\n", results[0][x], results[1][x]);
}

 */

/*
 * Table below incarnates two discrete Gaussian distribution:
 *    D(x) = exp(-((x - mu)^2)/(2*sigma^2))
 * where sigma = 1.292 and mu is 0 or 1 / 2.
 * Element 0 of the first table is P(x = 0) and 2*P(x = 1) in the second table.
 * For k > 0, element k is P(x >= k+1 | x > 0) in the first table, and
 * P(x >= k+2 | x > 1) in the second table.
 * For constant-time principle, mu = 0 is in the even indices and
 * mu = 1 / 2 is in the odd indices.
 * Probabilities are scaled up by 2^63.
 */
static const uint64_t gauss_1292[24] = {
	2847982254933138603u, 5285010687306232178u,
	3115855658194614154u, 2424313226695581870u,
	 629245045388085487u,  372648834165936922u,
	  73110737927091842u,   32559817584178793u,
	   4785625785139312u,    1592210133688742u,
	    174470148146634u,      43209976786070u,
	      3520594834759u,        647780323462u,
	        39186846585u,          5350987999u,
	          240149359u,            24322099u,
	             809457u,               60785u,
	               1500u,                  83u,
	                  2u,                   0u
};

/*
 * Sample an integer with parity equal to double_mu from a discrete Gaussian
 * distribution with support 2\ZZ + double_mu, mean 0 and sigma 2 * 1.292.
 * That is, an integer x (== double_mu mod 2) is chosen with probability
 * proportional to:
 *
 *     exp(- x^2 / (8 1.292^2)).
 */
static inline int
mkgauss_1292(prng *rng, uint8_t double_mu)
{
	uint64_t r;
	uint32_t v, k, neg;
	int32_t w;

	/*
	 * We use two random 64-bit values. First value is used to generate the
	 * non-zero value. Second value decides on whether the generated value is 0
	 * and the sign of the value. For constant-time code we have to read the
	 * complete table.
	 */

	r = prng_get_u64(rng) & ~((uint64_t)1 << 63);
	v = 1;
	for (k = 1; k < 12; k ++) {
		/*
		 * Add 1 if r < gauss_1292[2 * k + double_mu].
		 */
		v += (r - gauss_1292[2 * k + double_mu]) >> 63;
	}

	/*
	 * First value:
	 *  - flag 'neg' is randomly selected to be 0 or 1.
	 *  - if r/2^63 <= P(X == 0), then set v to zero.
	 */
	r = prng_get_u64(rng);
	neg = (uint32_t)(r >> 63);
	r &= ~((uint64_t)1 << 63);
	v &= -((gauss_1292[double_mu] - r) >> 63);

	/*
	 * We apply the sign ('neg' flag).
	 * If mu = 0 change v to v if neg = 0 and -v if neg = 1.
	 * If mu = 1/2, change v to v + 1 if neg = 0 and -v if neg = 1.
	 */
	v = (v ^ -neg) + neg + (~neg & double_mu);

	/*
	 * Now, transform the support of the sampler from Z to 2Z - double_mu, i.e.
	 * for mu = 1/2, we have -1, 1 with equal likelihood and -2, 2 with equal
	 * likelihood, etc.
	 */
	w = *(int32_t *)&v;
	return 2 * w - (int) double_mu;
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

static void
hash_to_fft(fpr *p, const uint8_t *h, unsigned logn)
{
	size_t n, u, v;
	uint8_t hash;

	n = MKN(logn);

	if (logn <= 3) {
		for (v = 0; v < n; v ++) {
			p[v] = fpr_of((h[0] >> v) & 1);
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
		for (v = 0; v < n; v ++) {
			x = fpr_rint(noise[v]);
			s[v] = (int16_t)((int64_t)((h[0] >> v) & 1) - x) / 2;
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
 * deviation 1.292. Returns whether or not (x0, x1) has a squared l2-norm less
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
		z = mkgauss_1292(rng, z);
		x0[u] = fpr_of(z);
		norm += z*z;
	}
	for (u = 0; u < n; u ++) {
		z = fpr_rint(x1[u]) & 1;
		z = mkgauss_1292(rng, z);
		x1[u] = fpr_of(z);
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
	Zf(FFT)(x0, logn);
	Zf(FFT)(x1, logn);

	return 1;
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
	uint16_t *bf, *bg, *bF, *bG, *x0, *x1;
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
		for (v = 0; v < n; v ++) {
			x0[v] = (h[0] >> v) & 1;
			x1[v] = (h[1] >> v) & 1;
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
		z = mkgauss_1292(rng, z & 1);
		x0[u] = Zf(mq_conv_small)(z);
		norm += z*z;
	}
	for (u = 0; u < n; u ++) {
		z = Zf(mq_conv_signed)(x1[u]);
		z = mkgauss_1292(rng, z & 1);
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
		for (v = 0; v < n; v ++) {
			s0[v] = (((h[0] >> v) & 1) - s0[v]) / 2;
			s1[v] = (((h[1] >> v) & 1) - s1[v]) / 2;
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
	uint16_t *bf, *bg, *bF, *bG, *x0, *x1;
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
		for (v = 0; v < n; v ++) {
			x0[v] = (h[0] >> v) & 1;
			x1[v] = (h[1] >> v) & 1;
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
		z = mkgauss_1292(rng, z & 1);
		x0[u] = Zf(mq_conv_small)(z);
		norm += z*z;
	}
	for (u = 0; u < n; u ++) {
		z = Zf(mq_conv_signed)(x1[u]);
		z = mkgauss_1292(rng, z & 1);
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
		for (v = 0; v < n; v ++) {
			s1[v] = (((h[1] >> v) & 1) - s1[v]) / 2;
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

