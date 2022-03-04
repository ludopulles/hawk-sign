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

/** To generate the values in the table below, run the following code:

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
 * Generate a random value with a Gaussian distribution centered on double_mu/2.
 * The RNG must be ready for extraction (already flipped).
 *
 * Distribution has standard deviation 1.292 sqrt(512/N).
 */
static int
mkgauss(void *samp_ctx, unsigned logn, uint8_t double_mu)
{
	unsigned u, g;
	int val;

	sampler_context *sc = (sampler_context *)samp_ctx;

	g = 1U << (9 - logn);
	val = 0;
	for (u = 0; u < g; u ++) {
		/*
		 * Each iteration generates one value with the
		 * Gaussian distribution for N = 512.
		 *
		 * We use two random 64-bit values. First value
		 * decides on whether the generated value is 0, and,
		 * if not, the sign of the value. Second random 64-bit
		 * word is used to generate the non-zero value.
		 *
		 * For constant-time code we have to read the complete
		 * table. This has negligible cost, compared with the
		 * remainder of the keygen process (solving the NTRU
		 * equation).
		 */
		uint64_t r;
		uint32_t f, v, k, neg;

		/*
		 * First value:
		 *  - flag 'neg' is randomly selected to be 0 or 1.
		 *  - flag 'f' is set to 1 if the generated value is zero,
		 *    or set to 0 otherwise.
		 */
		r = prng_get_u64(&sc->p);
		neg = (uint32_t)(r >> 63);
		r &= ~((uint64_t)1 << 63);
		f = (uint32_t)((r - gauss_1292[double_mu]) >> 63);

		/*
		 * We produce a new random 63-bit integer r, and go over
		 * the array, starting at index 1. We store in v the
		 * index of the first array element which is not greater
		 * than r, unless the flag f was already 1.
		 */
		v = 0;
		r = prng_get_u64(&sc->p);
		r &= ~((uint64_t)1 << 63);
		for (k = 1; k < 13; k ++) {
			uint32_t t;
			t = (uint32_t)((r - gauss_1292[2 * k + double_mu]) >> 63) ^ 1;
			v |= k & -(t & (f ^ 1));
			f |= t;
		}

		/*
		 * We apply the sign ('neg' flag). If the value is zero and mu = 0,
		 * the sign has no effect. Moreover, if mu = 1/2 and neg=0, add one.
		 */
		v = (v ^ -neg) + neg + (~neg & double_mu);

		/*
		 * Generated value is added to val.
		 */
		val += *(int32_t *)&v;
	}
	return val;
}

// =============================================================================

/*
 * Convert an integer polynomial (with small values) into the
 * representation with complex numbers.
 */
static void
smallints_to_fpr(fpr *r, const int8_t *t, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u ++) {
		r[u] = fpr_of(t[u]);
	}
}

/*
 * Sample a vector (x0, x1) that is congruent to (h*f, h*g) modulo 2 from
 * a Discrete Gaussian with lattice coset 2Z^{2n} + (h*f, h*g) and
 * standard deviation 1/isigma_sig.
 * Returns whether or not (x0, x1) has a squared l2-norm less than bound.
 * Note: tmp must have a size of at least 4n bytes.
 */
static int
sample_short(void *samp_ctx, int8_t *restrict x0, int8_t *restrict x1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict hm, fpr isigma_sig, uint32_t bound,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u;
	int32_t norm, z;
	fpr *t0, *t1, *t2;

	n = MKN(logn);
	norm = 0;
	t0 = (fpr *)tmp;
	t1 = t0 + n;
	t2 = t1 + n;

	/*
	 * Set the target vector to [hm, 0] * B (hm is the hashed message).
	 */
	for (u = 0; u < n; u++) {
		t0[u] = fpr_of(f[u] & 1);
		t1[u] = fpr_of(g[u] & 1);
		t2[u] = fpr_of(hm[u] & 1);
	}

	Zf(FFT)(t0, logn);
	Zf(FFT)(t1, logn);
	Zf(FFT)(t2, logn);

	Zf(poly_mul_fft)(t0, t2, logn);
	Zf(poly_mul_fft)(t1, t2, logn);

	Zf(iFFT)(t0, logn);
	Zf(iFFT)(t1, logn);

	/*
	 * Sample and write the result in (x0,x1). Gaussian smoothing is used
	 * to not reveal information on the secret basis.
	 */
	for (u = 0; u < n; u ++) {
		z = fpr_rint(t0[u]) & 1;
		// z = 2*Zf(sampler)(samp_ctx, fpr_half(fpr_of(z)), isigma_sig) - z;
		z = 2*mkgauss(samp_ctx, logn, z) - z;
		x0[u] = (int8_t) z;
		norm += z*z;
		// printf("%d ", z);
	}
	for (u = 0; u < n; u ++) {
		z = fpr_rint(t1[u]) & 1;
		// z = 2*Zf(sampler)(samp_ctx, fpr_half(fpr_of(z)), isigma_sig) - z;
		z = 2*mkgauss(samp_ctx, logn, z) - z;
		x1[u] = (int8_t) z;
		norm += z*z;
		// printf("%d ", z);
	}

	/*
	 * Test whether the l2-norm of (x0, x1) is below the given bound. The
	 * code below uses only 32-bit operations to compute the squared norm,
	 * since the max. value is 2n * 128^2 <= 2^24 (when logn <= 9).
	 * For a large enough verification margin, it is unlikely that the
	 * norm of the gaussian (x0, x1) is too large.
	 */
	return (uint32_t)norm <= bound;
}

/*
 * Compute signature of hm: a vector (s0, s1) that is close to (hm/2, 0) wrt Q.
 * Here, s0 is not returned.
 * When 0 is returned, signing failed.
 * When 1 is returned, there exists some s0 for which (s0, s1) is a signature
 * but reconstructing s0 might still fail.
 *
 * tmp must have room for at least 24 * 2^logn bytes
 */
static int
do_sign(void *samp_ctx, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict hm, fpr isigma_sig, uint32_t bound,
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

	if (!sample_short(samp_ctx, x0, x1, f, g, hm, isigma_sig, bound, logn,
			(uint8_t *)(x1 + n))) {
		printf("ERR\n");
		return 0;
	}

	/*
	 * Note that (x0, x1) == B (h, 0) (mod 2) so we obtain a lattice point
	 * (s0, s1) that is close to (hm/2, 0) wrt Q by calculating:
	 *     (s0, s1) = ((h, 0) - B^{-1} (x0, x1)) / 2.
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
		s1[u] = (int16_t) -fpr_rint(t0[u]);
		// shouldn't happen, except when FFT had rounding issues
		if (s1[u] & 1) { exit(1); return 0; }
		s1[u] /= 2;
	}

	return 1;
}

/*
 * Compute signature of hm: a vector (s0, s1) that is close to (hm/2, 0) wrt Q.
 * Here, s0 is not returned.
 * 1 is returned iff (s0, s1) is a valid signature and s0 can be recovered from
 * s1 with simple rounding.
 *
 * tmp must have room for at least 40 * 2^logn bytes
 */
static int
do_guaranteed_sign(void *samp_ctx, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const fpr *restrict q00, const int8_t *restrict hm, fpr isigma_sig,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp)
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

	if (!sample_short(samp_ctx, x0, x1, f, g, hm, isigma_sig, bound, logn,
			(uint8_t*)(x1 + n))) {
		return 0;
	}

	/*
	 * Note that (x0, x1) == B (h, 0) (mod 2) so we obtain a lattice point
	 * (s0, s1) that is close to (hm/2, 0) wrt Q by calculating:
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
	Zf(iFFT)(t4, logn); // err = (f^* x0 + g^* x1) / q00

	// If err is not in (-.5,.5)^n, simple rounding will fail
	for (u = 0; u < n; u++) {
		if (!fpr_lt(fpr_neg(fpr_onehalf), t4[u])
			|| !fpr_lt(t4[u], fpr_onehalf)) return 0;
	}

	Zf(poly_mul_fft)(tf, tx1, logn);
	Zf(poly_mul_fft)(tg, tx0, logn);
	Zf(poly_sub)(tf, tg, logn); // s1 = x1 f - x0 g.

	Zf(iFFT)(tf, logn);
	for (u = 0; u < n; u ++) {
		s1[u] = (int16_t) -fpr_rint(tf[u]);
		// shouldn't happen, except when FFT had rounding issues
		if (s1[u] & 1) return 0;
		s1[u] /= 2;
	}

	return 1;
}

/*
 * Compute signature of hm: a vector (s0, s1) that is close to (hm/2, 0) wrt Q.
 * If 1 is returned, (s0, s1) is a valid signature.
 *
 * tmp must have room for at least 40 * 2^logn bytes
 */
static int
do_complete_sign(void *samp_ctx,
	int16_t *restrict s0, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const int8_t *restrict hm, fpr isigma_sig, uint32_t bound,
	unsigned logn, uint8_t *restrict tmp)
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

	if (!sample_short(samp_ctx, x0, x1, f, g, hm, isigma_sig, bound, logn,
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
		s0[u] = (int16_t)fpr_rint(t0[u]);
		// shouldn't happen, except when FFT had rounding issues
		if ((s0[u] ^ hm[u]) & 1) return 0;
		s0[u] = ((hm[u] & 1) - s0[u]) / 2;
	}
	for (u = 0; u < n; u ++) {
		s1[u] = (int16_t) -fpr_rint(t1[u]);
		// shouldn't happen, except when FFT had rounding issues
		if (s1[u] & 1) return 0;
		s1[u] /= 2;
	}
	// Now (s0, s1) is a lattice that is close to (h/2, 0) wrt Q.
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
	const int8_t *restrict hm, fpr isigma_sig, uint32_t bound,
	unsigned logn, uint8_t *restrict tmp)
{
	sampler_context spc;
	spc.sigma_min = fpr_sigma_min[logn];
	do {
		Zf(prng_init)(&spc.p, rng);
	} while (!do_sign(&spc, sig, f, g, hm, isigma_sig,
			bound, logn, tmp));
}

/* see inner.h */
void
Zf(guaranteed_sign)(inner_shake256_context *rng, int16_t *restrict sig,
	const int8_t *restrict f, const int8_t *restrict g,
	const fpr *restrict q00, const int8_t *restrict hm, fpr isigma_sig,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp)
{
	sampler_context spc;
	spc.sigma_min = fpr_sigma_min[logn];
	do {
		Zf(prng_init)(&spc.p, rng);
	} while (!do_guaranteed_sign((void *)&spc, sig, f, g, q00, hm,
			isigma_sig, bound, logn, tmp));
}

/* see inner.h */
void
Zf(complete_sign)(inner_shake256_context *rng,
	int16_t *restrict s0, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const int8_t *restrict hm, fpr isigma_sig, uint32_t bound,
	unsigned logn, uint8_t *restrict tmp)
{
	sampler_context spc;
	spc.sigma_min = fpr_sigma_min[logn];
	do {
		Zf(prng_init)(&spc.p, rng);
	} while (!do_complete_sign((void *)&spc, s0, s1, f, g, F, G,
			hm, isigma_sig, bound, logn, tmp));
}
