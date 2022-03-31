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
 * where sigma = 1.292 and mu is 0 or 1/2.
 * Element 0 of the first table is P(x = 0) and 2*P(x = 1) in the second table.
 * For k > 0, element k is P(x >= k+1 | x > 0) in the first table, and
 * P(x >= k+2 | x > 1) in the second table.
 * For constant-time principle, mu = 0 is in the even indices and
 * mu = 1/2 is in the odd indices.
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
	for (size_t u = MKN(logn); u --> 0; ) {
		r[u] = fpr_of(t[u]);
	}
}

/*
 * Sample a vector (x0, x1) that is congruent to (t0, t1) = B * (h0, h1)^t modulo 2 from
 * a Discrete Gaussian with lattice coset 2Z^{2n} + (t0, t1) and standard deviation 1.292.
 * Returns whether or not (x0, x1) has a squared l2-norm less than bound.
 *
 * Note: tmp[] must have space for at least 48 * 2^logn bytes.
 */
static int
sample_short(prng *rng, int8_t *restrict x0, int8_t *restrict x1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const int8_t *restrict h0, const int8_t *restrict h1,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u;
	int32_t norm, z;
	fpr *t0, *t1, *t2, *t3, *t4, *t5;

	n = MKN(logn);
	norm = 0;
	t0 = (fpr *)tmp;
	t1 = t0 + n;
	t2 = t1 + n;
	t3 = t2 + n;
	t4 = t3 + n;
	t5 = t4 + n;

	/*
	 * Set the target vector to B * (h0, h1)^T.
	 */
	for (u = 0; u < n; u++) {
		t0[u] = fpr_of(h0[u] & 1);
		t1[u] = fpr_of(h1[u] & 1);
		t2[u] = fpr_of(f[u] & 1);
		t3[u] = fpr_of(g[u] & 1);
		t4[u] = fpr_of(F[u] & 1);
		t5[u] = fpr_of(G[u] & 1);
	}

	Zf(FFT)(t0, logn);
	Zf(FFT)(t1, logn);
	Zf(FFT)(t2, logn);
	Zf(FFT)(t3, logn);
	Zf(FFT)(t4, logn);
	Zf(FFT)(t5, logn);

	Zf(poly_mul_fft)(t2, t0, logn); // f h0
	Zf(poly_mul_fft)(t3, t0, logn); // g h0
	Zf(poly_mul_fft)(t4, t1, logn); // F h1
	Zf(poly_mul_fft)(t5, t1, logn); // G h1
	Zf(poly_add)(t2, t4, logn); // f h0 + F h1
	Zf(poly_add)(t3, t5, logn); // g h0 + G h1

	Zf(iFFT)(t2, logn);
	Zf(iFFT)(t3, logn);

	/*
	 * Sample and write the result in (x0, x1). Gaussian smoothing is used to
	 * not reveal information on the secret basis.
	 */
	for (u = 0; u < n; u ++) {
		z = fpr_rint(t2[u]) & 1;
		z = mkgauss_1292(rng, z);
		x0[u] = (int8_t) z;
		norm += z*z;
	}
	for (u = 0; u < n; u ++) {
		z = fpr_rint(t3[u]) & 1;
		z = mkgauss_1292(rng, z);
		x1[u] = (int8_t) z;
		norm += z*z;
	}

	/*
	 * Test whether the l2-norm of (x0, x1) is below the given bound. The
	 * code below uses only 32-bit operations to compute the squared norm,
	 * since the max. value is 2n * 128^2 <= 2^24 (when logn <= 9).
	 * For a large enough verification margin, it is unlikely that the
	 * norm of the gaussian (x0, x1) is too large.
	 */
	return (uint32_t)norm <= l2bound[logn];
}

/*
 * Compute signature of (h0, h1): a vector (s0, s1) that is close to (h0, h1)/2 wrt Q.
 * Here, s0 is not returned.
 * When 0 is returned, signing failed.
 * When 1 is returned, there exists some s0 for which (s0, s1) is a signature
 * but reconstructing s0 might still fail.
 *
 * Note: tmp[] must have space for at least 50 * 2^logn bytes.
 */
static int
do_sign(prng *rng, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const int8_t *restrict h0, const int8_t *restrict h1,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u;
	int8_t *x0, *x1;
	fpr *t0, *t1, *t2;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;
	t2 = t1 + n;
	x0 = (int8_t *)tmp;
	x1 = x0 + n;

	if (!sample_short(rng, x0, x1, f, g, F, G, h0, h1, logn,
			(uint8_t *)(x1 + n))) {
		return 0;
	}

	/*
	 * Note that (x0, x1) == B (h0, h1) (mod 2) so we obtain a lattice point
	 * (s0, s1) that is close to (h0, h1)/2 wrt Q by calculating:
	 *     (s0, s1) = ((h0, h1) - B^{-1} (x0, x1)) / 2.
	 */
	smallints_to_fpr(t1, x1, logn);
	smallints_to_fpr(t2, x0, logn);
	// Now override x0, x1
	smallints_to_fpr(t0, f, logn);
	Zf(FFT)(t0, logn);
	Zf(FFT)(t1, logn);
	Zf(FFT)(t2, logn);
	Zf(poly_mul_fft)(t0, t1, logn);
	smallints_to_fpr(t1, g, logn);
	Zf(FFT)(t1, logn);
	Zf(poly_mul_fft)(t1, t2, logn);
	Zf(poly_sub)(t0, t1, logn); // s1 = x1 f - x0 g.

	Zf(iFFT)(t0, logn);
	for (u = 0; u < n; u ++) {
		s1[u] = (int16_t) (h1[u] - fpr_rint(t0[u])) / 2;
	}

	return 1;
}

/*
 * Compute signature of (h0, h1): a vector (s0, s1) that is close to (h0, h1)/2 wrt Q.
 * Here, s0 is not returned.
 * 1 is returned iff (s0, s1) is a valid signature and s0 can be recovered from
 * s1 with simple rounding.
 *
 * Note: tmp[] must have space for at least 50 * 2^logn bytes.
 */
static int
do_guaranteed_sign(prng *rng, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const fpr *restrict q00, const int8_t *restrict h0, const int8_t *restrict h1,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u;
	int8_t *x0, *x1;
	fpr *tf, *tg, *tx0, *tx1, *t4;

	n = MKN(logn);
	tf = (fpr *)tmp;
	tg = tf + n;
	tx0 = tg + n;
	tx1 = tx0 + n;
	t4 = tx1 + n;

	x0 = (int8_t *)tmp;
	x1 = x0 + n;

	if (!sample_short(rng, x0, x1, f, g, F, G, h0, h1, logn,
			(uint8_t*)(x1 + n))) {
		return 0;
	}

	/*
	 * Note that (x0, x1) == B (h, 0) (mod 2) so we obtain a lattice point
	 * (s0, s1) that is close to (h0, h1)/2 wrt Q by calculating:
	 *     (s0, s1) = ((h, 0) - B^{-1} (x0, x1)) / 2.
	 */
	smallints_to_fpr(tx0, x0, logn);
	smallints_to_fpr(tx1, x1, logn);
	// Now override x0, x1:
	smallints_to_fpr(tf, f, logn);
	smallints_to_fpr(tg, g, logn);

	Zf(FFT)(tx0, logn);
	Zf(FFT)(tx1, logn);
	Zf(FFT)(tf, logn);
	Zf(FFT)(tg, logn);

	Zf(poly_add_muladj_fft)(t4, tx0, tx1, tf, tg, logn);
	Zf(poly_div_autoadj_fft)(t4, q00, logn);
	Zf(iFFT)(t4, logn); // err / 2 = (f^* x0 + g^* x1) / q00

	// If err / 2 is not in (-.5,.5)^n, simple rounding will fail
	for (u = 0; u < n; u++) {
		if (!fpr_lt(fpr_neg(fpr_one), t4[u]) || !fpr_lt(t4[u], fpr_one)) {
			return 0;
		}
	}

	Zf(poly_mul_fft)(tf, tx1, logn);
	Zf(poly_mul_fft)(tg, tx0, logn);
	Zf(poly_sub)(tf, tg, logn); // s1 = x1 f - x0 g.

	Zf(iFFT)(tf, logn);
	for (u = 0; u < n; u ++) {
		s1[u] = (int16_t) (h1[u] - fpr_rint(tf[u])) / 2;
	}

	return 1;
}

/*
 * Compute signature of h: a vector (s0, s1) that is close to (h0, h1)/2 wrt Q.
 * Here, s0 is not returned.
 * 1 is returned iff (s0, s1) is a valid signature and s0 can be recovered from
 * s1 with simple rounding.
 *
 * Note: tmp[] must have space for at least 26 * 2^logn bytes.
 */
static int
do_fft_sign(prng *rng, int16_t *restrict s1,
	const fpr *restrict expanded_seckey, const uint8_t *restrict h,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u, v, w;
	int8_t *x0, *x1;
	const fpr *tf, *tg, *tF, *tG, *tiq00;
	fpr *tx0, *tx1, *tres;
	int32_t norm, z;


	n = MKN(logn);
	norm = 0;

	tf = expanded_seckey;
	tg = tf + n;
	tF = tg + n;
	tG = tF + n;
	tiq00 = tG + n;

	tx0 = (fpr *)tmp;
	tx1 = tx0 + n;
	tres = tx1 + n;
	x0 = (int8_t *)(tres + n);
	x1 = x0 + n;

	/*
	 * Set the target vector to (h0, h1) B (row-notation).
	 */
	for (u = 0, w = 0; w < n; u ++ ) {
		uint8_t h0, h1;

		h0 = h[u];
		h1 = h[n / 8 + u];
		for (v = 0; v < 8; v ++, w ++) {
			tx0[w] = fpr_of(h0 & 1);
			tx1[w] = fpr_of(h1 & 1);
			h0 >>= 1;
			h1 >>= 1;
		}
	}

	/*
	 * Sample (x0, x1) according to a discrete gaussian distribution on 2
	 * Z^{2n} + (h0, h1) B with standard deviation 2 sigma_{sig}. Gaussian
	 * smoothing is used to not reveal information on the secret basis.
	 */
	Zf(FFT)(tx0, logn);
	Zf(FFT)(tx1, logn);

	/*
	 * First, calculate the 0th component of (h0, h1) B, and sample x0.
	 */
	Zf(poly_add_mul_fft)(tres, tx0, tx1, tf, tF, logn);
	Zf(iFFT)(tres, logn);

	for (u = 0; u < n; u ++) {
		z = fpr_rint(tres[u]) & 1;
		z = mkgauss_1292(rng, z);
		x0[u] = (int8_t) z;
		norm += z*z;
	}

	/*
	 * Second, calculate the 1th component of (h0, h1) B, and sample x1.
	 */
	Zf(poly_add_mul_fft)(tres, tx0, tx1, tg, tG, logn);
	Zf(iFFT)(tres, logn);

	for (u = 0; u < n; u ++) {
		z = fpr_rint(tres[u]) & 1;
		z = mkgauss_1292(rng, z);
		x1[u] = (int8_t) z;
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
		// Norm is too large, so signature would not be valid.
		return 0;
	}

	smallints_to_fpr(tx0, x0, logn);
	smallints_to_fpr(tx1, x1, logn);
	Zf(FFT)(tx0, logn);
	Zf(FFT)(tx1, logn);

	/*
	 * Calculate the rounding errors that occur when we want to recover s0 from
	 * s1 during verification of this signature. These errors should be of
	 * absolute value at most 1/2.
	 * The errors are calculated by:
	 *
	 *     (f^* x0 + g^* x1) / (f^* f + g^* g).
	 *
	 * If err / 2 is not in (-.5,.5)^n, verification will reconstruct a
	 * different s0 so the signature may not be valid anymore.
	 */
	Zf(poly_add_muladj_fft)(tres, tx0, tx1, tf, tg, logn);
	Zf(poly_mul_autoadj_fft)(tres, tiq00, logn);
	Zf(iFFT)(tres, logn);

	for (u = 0; u < n; u++) {
		if (!fpr_lt(fpr_neg(fpr_one), tres[u]) || !fpr_lt(tres[u], fpr_one)) {
			return 0;
		}
	}

	/*
	 * In row-notation, (x0, x1) == (h0, h1) B (mod 2) holds while it also has
	 * a small norm since the discrete gaussian distribution has small standard
	 * deviation. Thus, (s0, s1) := ((h0, h1) - (x0, x1) B^{-1}) / 2 is a
	 * lattice point that is close to (h0, h1)/2 with respect to the quadratic
	 * form Q.
	 *
	 * Note: B = [[f, g], [F, G]] with det. 1 so B^{-1} = [[G, -g], [-F, f]]
	 * yielding s1 = (h1 - (-x0 g + x1 f)) / 2.
	 */
	Zf(poly_neg)(tx0, logn);
	Zf(poly_add_mul_fft)(tres, tx0, tx1, tg, tf, logn);
	Zf(iFFT)(tres, logn);

	for (u = 0, w = 0; w < n; u ++) {
		uint8_t h1;

		h1 = h[n / 8 + u];
		for (v = 0; v < 8; v ++, w ++) {
			s1[w] = ((int16_t)(h1 & 1) - fpr_rint(tres[w])) / 2;
			h1 >>= 1;
		}
	}

	return 1;
}

/*
 * Compute signature of (h0, h1): a vector (s0, s1) that is close to
 * (h0, h1)/2 wrt Q. If 1 is returned, (s0, s1) is a valid signature.
 *
 * Note: tmp[] must have space for at least 50 * 2^logn bytes
 */
static int
do_complete_sign(prng *rng,
	int16_t *restrict s0, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const int8_t *restrict h0, const int8_t *restrict h1, unsigned logn,
	uint8_t *restrict tmp)
{
	size_t n, u;
	int8_t *x0, *x1;
	fpr *t0, *t1, *t2, *t3, *t4;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;
	t2 = t1 + n;
	t3 = t2 + n;
	t4 = t3 + n;
	x0 = (int8_t *)tmp;
	x1 = x0 + n;

	if (!sample_short(rng, x0, x1, f, g, F, G, h0, h1, logn,
			(uint8_t *)(x1 + n))) {
		return 0;
	}

	/*
	 * Get the signature corresponding to that tiny vector, i.e.
	 * s = x * B^{-1}. Thus s0 = x0 G - x1 F and s1 = -x0 g + x1 f.
	 */
	smallints_to_fpr(t2, x0, logn);
	smallints_to_fpr(t3, x1, logn);
	// Now override x0, x1
	smallints_to_fpr(t0, G, logn);
	smallints_to_fpr(t1, f, logn);
	Zf(FFT)(t0, logn);
	Zf(FFT)(t1, logn);
	Zf(FFT)(t2, logn);
	Zf(FFT)(t3, logn);

	Zf(poly_mul_fft)(t0, t2, logn); // t0 = x0 G
	Zf(poly_mul_fft)(t1, t3, logn); // t1 = x1 f

	smallints_to_fpr(t4, F, logn);
	Zf(FFT)(t4, logn);
	Zf(poly_mul_fft)(t3, t4, logn);
	Zf(poly_sub)(t0, t3, logn); // t0 = x0 G - x1 F

	smallints_to_fpr(t4, g, logn);
	Zf(FFT)(t4, logn);
	Zf(poly_mul_fft)(t2, t4, logn);
	Zf(poly_sub)(t1, t2, logn); // t1 = x1 f - x0 g

	/*
	 * Extract the signature from t0, t1
	 */
	Zf(iFFT)(t0, logn);
	Zf(iFFT)(t1, logn);

	for (u = 0; u < n; u ++) {
		s0[u] = (int16_t) (h0[u] - fpr_rint(t0[u])) / 2;
		s1[u] = (int16_t) (h1[u] - fpr_rint(t1[u])) / 2;
	}
	// Now (s0, s1) is a lattice that is close to (h0, h1)/2 wrt Q.
	return 1;
}

/* =================================================================== */
/*
 * Use a fast PRNG for gaussian sampling during signing.
 */

/* see inner.h */
void
Zf(sign)(inner_shake256_context *rng, int16_t *restrict sig,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const int8_t *restrict h0, const int8_t *restrict h1,
	unsigned logn, uint8_t *restrict tmp)
{
	prng p;
	do {
		Zf(prng_init)(&p, rng);
	} while (!do_sign(&p, sig, f, g, F, G, h0, h1, logn, tmp));
}

/* see inner.h */
void
Zf(guaranteed_sign)(inner_shake256_context *rng, int16_t *restrict sig,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const fpr *restrict q00, const int8_t *restrict h0,
	const int8_t *restrict h1, unsigned logn,
	uint8_t *restrict tmp)
{
	prng p;
	do {
		Zf(prng_init)(&p, rng);
	} while (!do_guaranteed_sign(&p, sig, f, g, F, G, q00, h0, h1, logn, tmp));
}

/* see inner.h */
void
Zf(complete_sign)(inner_shake256_context *rng,
	int16_t *restrict s0, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const int8_t *restrict h0, const int8_t *restrict h1,
	unsigned logn, uint8_t *restrict tmp)
{
	prng p;
	do {
		Zf(prng_init)(&p, rng);
	} while (!do_complete_sign(&p, s0, s1, f, g, F, G, h0, h1, logn, tmp));
}


/* see inner.h */
void
Zf(expand_seckey)(fpr *restrict expanded_seckey,
	const int8_t *f, const int8_t *g, const int8_t *F, unsigned logn)
{
	size_t n, hn, u;
	fpr *bf, *bg, *bF, *bG, *biq00;

	n = MKN(logn);
	hn = n >> 1;

	bf = expanded_seckey;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	biq00 = bG + n;

	/*
	 * We load the private key elements directly into the 2x2 matrix B.
	 */
	smallints_to_fpr(bf, f, logn);
	smallints_to_fpr(bg, g, logn);
	smallints_to_fpr(bF, F, logn);

	/*
	 * Compute the FFT for the key elements
	 */
	Zf(FFT)(bf, logn);
	Zf(FFT)(bg, logn);
	Zf(FFT)(bF, logn);

	/*
	 * Compute G = (1 + gF) / f
	 */
	Zf(poly_prod_fft)(bG, bg, bF, logn);
	for (u = 0; u < hn; u++)
		bG[u] = fpr_add(bG[u], fpr_one);
	Zf(poly_div_fft)(bG, bf, logn);

	Zf(poly_invnorm2_fft)(biq00, bf, bg, logn);
}

/* see inner.h */
void
Zf(fft_sign)(inner_shake256_context *rng, int16_t *restrict sig,
	const fpr *restrict expanded_seckey, const uint8_t *restrict h,
	unsigned logn, uint8_t *restrict tmp)
{
	prng p;
	do {
		Zf(prng_init)(&p, rng);
	} while (!do_fft_sign(&p, sig, expanded_seckey, h, logn, tmp));
}

