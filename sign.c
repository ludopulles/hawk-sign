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

#include <assert.h>
#include "inner.h"

// =============================================================================

/*
 * This number indicates the maximum l2-norm that is allowed for the small
 * error (from a valid signature) that is chosen around the target point during
 * signature generation.
 * The coefficients of this error are distributed according to a discrete
 * gaussian over the integers that have the same parity as the target
 * coefficient and with standard deviation 2 sigma_sig. To compute these
 * values, use:
 *
 *     l2bound(logn) = floor( (verif_margin * 2 sigma_sig)^2 * 2n ).
 *
 * Here, we have taken verif_margin = 1.1 and sigma_sig = 1.292.
 */
static const uint32_t l2bound[10] = {
	0 /* unused */, 32, 64, 129, 258, 517, 1034, 2068, 4136, 8273
};

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
static const uint64_t gauss_1292[26] = {
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
	                  2u,                   0u,
	                  0u,                   0u,
};

/*
 * Sample a integer with parity equal to double_mu from a Gaussian distribution
 * with standard deviation 2 * 1.292, i.e. a value x == parity (mod 2) is
 * chosen with probability proportional to
 *
 *     exp(- x^2 / (8 1.292^2)).
 *
 * Note: The RNG must be ready for extraction (already flipped).
 */
static inline int
mkgauss_1292(prng *rng, uint8_t double_mu)
{
	/* We use two random 64-bit values. First value decides on whether the
	 * generated value is 0, and, if not, the sign of the value. Second
	 * random 64-bit word is used to generate the non-zero value.
	 *
	 * For constant-time code we have to read the complete table. Currently
	 * sampling takes up most time of the signing process, so this is
	 * something that should be optimized.
	 */
	uint64_t r;
	uint32_t f, v, k, neg;

	/*
	 * First value:
	 *  - flag 'neg' is randomly selected to be 0 or 1.
	 *  - flag 'f' is set to 1 if the generated value is zero,
	 *    or set to 0 otherwise.
	 */
	r = prng_get_u64(rng);

	neg = (uint32_t)(r >> 63);
	r &= ~((uint64_t)1 << 63);
	f = (uint32_t)((r - gauss_1292[double_mu]) >> 63);

	/*
	 * We produce a new random 63-bit integer r, and go over
	 * the array, starting at index 1. The first time an array element with
	 * value less than or equal to r, v is set to the 'right' value, and f
	 * is set to 1 so v only gets a value once.
	 */
	v = double_mu;

	r = prng_get_u64(rng);

	r &= ~((uint64_t)1 << 63);
	for (k = 1; k < 13; k ++) {
		uint32_t t;
		t = (uint32_t)((r - gauss_1292[2 * k + double_mu]) >> 63) ^ 1;
		v |= (k << 1) & -(t & (f ^ 1));
		f |= t;
	}

	/*
	 * We apply the sign ('neg' flag).
	 */
	v = (v ^ -neg) + neg;

	/*
	 * Return v interpreted as a signed integer
	 */
	return *(int32_t *)&v;
}

// =============================================================================

/*
 * Convert an integer polynomial (with small values) into the
 * representation with complex numbers.
 * Also works when r and t overlap, since the loop goes in decreasing order.
 */
static void
smallints_to_fpr(fpr *r, const int8_t *t, unsigned logn)
{
	size_t u;

	u = MKN(logn);
	while (u --> 0) {
		r[u] = fpr_of(t[u]);
	}
}

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
hash_to_fpr(fpr *p, const uint8_t *h, unsigned logn)
{
	size_t n, u, v;
	uint8_t hash;

	n = MKN(logn);

	if (logn <= 3) {
		hash = *h;
		for (v = 0; v < n; v ++) {
			p[v] = fpr_of(hash & 1);
			hash >>= 1;
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

	n = MKN(logn);
	if (logn <= 3) {
		hash = *h;
		for (v = 0; v < n; v ++) {
			s[v] = (int16_t) ((int64_t)(hash & 1) - fpr_rint(noise[v])) / 2;
			hash >>= 1;
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				s[u] = (int16_t) ((int64_t)(hash & 1) - fpr_rint(noise[u])) / 2;
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
	smallints_to_fpr(bf, f, logn);
	smallints_to_fpr(bg, g, logn);
	smallints_to_fpr(bF, F, logn);
	Zf(FFT)(bf, logn);
	Zf(FFT)(bg, logn);
	Zf(FFT)(bF, logn);

	if (G == NULL) {
		size_t u, hn;

		/*
		 * Compute G = (1 + gF) / f, where all polynomials are in FFT representation.
		 */
		hn = MKN(logn) >> 1;
		Zf(poly_prod_fft)(bG, bg, bF, logn);
		for (u = 0; u < hn; u++) {
			bG[u] = fpr_add(bG[u], fpr_one);
		}
		Zf(poly_div_fft)(bG, bf, logn);
	} else {
		smallints_to_fpr(bG, G, logn);
		Zf(FFT)(bG, logn);
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

	hash_to_fpr(x0, h, logn);
	hash_to_fpr(x1, SECOND_HASH(h, logn), logn);
	Zf(FFT)(x0, logn);
	Zf(FFT)(x1, logn);

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
	if ((uint32_t)norm > l2bound[logn]) {
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
 * Compute signature of (h0, h1): a vector (s0, s1) that is close to
 * (h0, h1) / 2 wrt Q. Here, s0 is not returned.
 * 1 is returned iff (s0, s1) is a valid signature and s0 can be recovered from
 * s1 with simple rounding.
 *
 * Note: tmp[] must have space for at least 48 * 2^logn bytes.
 */
static int
do_sign_dyn(prng *rng, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u;
	fpr *x0, *x1, *bf, *bg, *bF, *bG;

	n = MKN(logn);
	bf = (fpr *) tmp;
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

/*
 * Compute signature of h: a vector (s0, s1) that is close to (h0, h1) / 2 wrt Q.
 * Here, s0 is not returned.
 * 1 is returned iff (s0, s1) is a valid signature and s0 can be recovered from
 * s1 with simple rounding.
 *
 * Note: tmp[] must have space for at least 24 * 2^logn bytes.
 */
static int
do_sign(prng *rng, int16_t *restrict s1,
	const fpr *restrict expanded_seckey, const uint8_t *restrict h,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u;
	const fpr *bf, *bg, *bF, *bG, *biq00;
	fpr *x0, *x1, *res;

	n = MKN(logn);

	bf = expanded_seckey;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	biq00 = bG + n;

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
	Zf(poly_mul_autoadj_fft)(res, biq00, logn);
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

/*
 * Compute signature of (h0, h1): a vector (s0, s1) that is close to
 * (h0, h1) / 2 wrt Q. If 1 is returned, (s0, s1) is a valid signature.
 *
 * Note: tmp[] must have space for at least 50 * 2^logn bytes
 */
static int
do_complete_sign(prng *rng,
	int16_t *restrict s0, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	size_t n;
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
	 * Calculate B^{-1} (x0, x1) = ((-G) x0 + F x1, g x0 + (-f) x1).
	 */
	Zf(poly_neg)(bF, logn);
	Zf(poly_neg)(bg, logn);
	Zf(poly_matmul_fft)(bG, bF, bg, bf, x0, x1, logn);

	/*
	 * Extract the signature from x0, t1
	 */
	Zf(iFFT)(x0, logn);
	Zf(iFFT)(x1, logn);

	noise_to_lattice(s0, h, x0, logn);
	noise_to_lattice(s1, SECOND_HASH(h, logn), x1, logn);
	return 1;
}

/* =================================================================== */

/* see inner.h */
void
Zf(expand_seckey)(fpr *restrict expanded_seckey,
	const int8_t *f, const int8_t *g, const int8_t *F, unsigned logn)
{
	size_t n;
	fpr *bf, *bg, *bF, *bG, *biq00;

	n = MKN(logn);
	bf = expanded_seckey;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	biq00 = bG + n;

	/*
	 * We load the private key elements directly into the 2x2 matrix B.
	 */
	construct_basis(f, g, F, NULL, bf, bg, bF, bG, logn);
	Zf(poly_invnorm2_fft)(biq00, bf, bg, logn);
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
Zf(complete_sign)(inner_shake256_context *rng,
	int16_t *restrict s0, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	prng p;
	do {
		Zf(prng_init)(&p, rng);
	} while (!do_complete_sign(&p, s0, s1, f, g, F, G, h, logn, tmp));
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

