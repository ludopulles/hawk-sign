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

	for (u = 0; u < n; u++) {
		x0[u] = fpr_rint(t0[u]) & 1;
		x1[u] = fpr_rint(t1[u]) & 1;
	}

	/*
	 * Sample and write the result in (x0,x1). Gaussian smoothing is used
	 * to not reveal information on the secret basis.
	 */
	for (u = 0; u < n; u ++) {
		z = 2*Zf(sampler)(samp_ctx, fpr_half(fpr_of(x0[u])), isigma_sig) - x0[u];
		x0[u] = (int8_t) z;
		norm += z*z;
	}
	for (u = 0; u < n; u ++) {
		z = 2*Zf(sampler)(samp_ctx, fpr_half(fpr_of(x1[u])), isigma_sig) - x1[u];
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
		if (s1[u] & 1) return 0;
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
