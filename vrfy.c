/*
 * Hawk signature verification.
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
 * Add to a polynomial its own adjoint. This function works only in FFT
 * representation.
 */
void
Zf(poly_addselfadj_fft)(fpr *a, unsigned logn)
{
	size_t hn, u;

	hn = MKN(logn) >> 1;
	for (u = 0; u < hn; u ++)
		a[u] = fpr_double(a[u]);
	for (u = 0; u < hn; u ++)
		a[u + hn] = fpr_zero;
}

/*
 * Add polynomial b to polynomial a, where b is autoadjoint. Both a and b are
 * in FFT representation. Since b is autoadjoint, all its FFT coefficients are
 * real, and the array b contains only N/2 elements.
 */
void
Zf(poly_add_autoadj_fft)(fpr *a, fpr *b, unsigned logn)
{
	size_t hn, u;

	hn = MKN(logn) >> 1;
	for (u = 0; u < hn; u ++)
		a[u] = fpr_add(a[u], b[u]);
}

/* =================================================================== */

/* see inner.h */
int
Zf(complete_verify)(const int8_t *restrict h0, const int8_t *restrict h1,
	const int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp)
{
	size_t u, hn, n;
	fpr *t0, *t1, trace;

	n = MKN(logn);
	hn = n >> 1;
	t0 = (fpr *)tmp;
	t1 = t0 + n;

	/**
	 * Put (t0, t1) = (2 s0 - h0, 2 s1 - h1) in FFT representation so we can
	 * calculate the l2-norm wrt Q.
	 */
	for (u = 0; u < n; u ++) {
		t0[u] = fpr_of(2 * s0[u] - (h0[u] & 1));
	}
	for (u = 0; u < n; u ++) {
		t1[u] = fpr_of(2 * s1[u] - (h1[u] & 1));
	}

	// Takes 50% of time (2us of 4us):
	Zf(FFT)(t0, logn);
	Zf(FFT)(t1, logn);

	trace = fpr_zero;

	// Calculate normalized trace( (t0, t1)^* Q (t0, t1) )
	// Note that the trace is by definition the sum under its n embeddings
	for (u = 0; u < hn; u ++) {
		fpr a_re, a_im, b_re, b_im, ab_re, ab_im;

		a_re = t0[u];
		a_im = t0[u + hn];
		b_re = t1[u];
		b_im = t1[u + hn];

		fpr norm0 = fpr_add(fpr_sqr(a_re), fpr_sqr(a_im));
		fpr norm1 = fpr_add(fpr_sqr(b_re), fpr_sqr(b_im));

		// contribution of t0 Q00 t0*
		trace = fpr_add(trace, fpr_mul(norm0, q00[u]));
		// contribution of t1 Q11 t1*
		trace = fpr_add(trace, fpr_mul(norm1, q11[u]));

		// ab = t1 t0*
		ab_re = fpr_add(fpr_mul(a_re, b_re), fpr_mul(a_im, b_im));
		ab_im = fpr_sub(fpr_mul(a_re, b_im), fpr_mul(a_im, b_re));

		// contribution of t1 q10 t0* + t0 q01 t1*
		trace = fpr_add(trace, fpr_double(
			fpr_sub(fpr_mul(ab_re, q10[u]), fpr_mul(ab_im, q10[u + hn]))
		));
	}

	/*
	 * Note: only n/2 embeddings are stored, because they come in pairs.
	 */
	trace = fpr_double(trace);
	/*
	 * Renormalize, so we get the geometric norm of (t0, t1) w.r.t Q.
	 */
	trace = fpr_div(trace, fpr_of(n));

	/*
	 * First check whether the norm is in range [0, 2^31).
	 * Signature is valid iff
	 *     `Tr(s* Q s) / n (=Tr(x^* x)/n = sum_i x_i^2) <= bound`.
	 */
	return fpr_lt(fpr_zero, trace) && fpr_lt(trace, fpr_ptwo31m1)
		&& (uint32_t)fpr_rint(trace) <= bound;
}

/* see inner.h */
int
Zf(verify_simple_rounding)(const int8_t *restrict h0, const int8_t *restrict h1,
	int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp)
{
	size_t u, n;
	fpr *t0;

	n = MKN(logn);
	t0 = (fpr *)tmp;

	for (u = 0; u < n; u ++) {
		t0[u] = fpr_half(fpr_of((h1[u] & 1) - 2 * s1[u]));
	}

	Zf(FFT)(t0, logn);
	Zf(poly_mul_fft)(t0, q10, logn); // (h1/2 - s1) q10
	Zf(poly_div_autoadj_fft)(t0, q00, logn); // (h1/2 - s1) q10/q00
	Zf(iFFT)(t0, logn);

	// Recover s0 with s0 = round(h0/2 + (h1/2 - s1) q10 / q00)
	for (u = 0; u < n; u ++) {
		s0[u] = fpr_rint(fpr_add(fpr_half(fpr_of(h0[u] & 1)), t0[u]));
	}

	return Zf(complete_verify)(h0, h1, s0, s1, q00, q10, q11, bound, logn, tmp);
}

/* see inner.h */
int
Zf(verify_nearest_plane)(const int8_t *restrict h0, const int8_t *restrict h1,
	int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp)
{
	/*
	 * This works better than simple rounding.
	 * Reconstruct s0, by running Babai's NP algorithm with target
	 *     h0 / 2 + (h1 / 2 - s1) * q10 / q00.
	 */

	size_t n, u;
	fpr *t0, *t1;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;
	for (u = 0; u < n; u ++) {
		t0[u] = fpr_half(fpr_of(h0[u] & 1));
	}
	for (u = 0; u < n; u ++) {
		t1[u] = fpr_half(fpr_of((h1[u] & 1) - 2 * s1[u]));
	}

	Zf(FFT)(t0, logn);
	Zf(FFT)(t1, logn);
	Zf(poly_mul_fft)(t1, q10, logn);
	Zf(poly_div_fft)(t1, q00, logn);
	Zf(poly_add)(t0, t1, logn); // t0 = h0 / 2 + (h1 / 2 - s1) q10/q00

	memcpy(t1, q00, n * sizeof(fpr));
	// Run Babai with target t0 and Gram-matrix q00.
	Zf(ffNearestPlane_dyn)(t0, t1, logn, t1 + n);
	Zf(iFFT)(t0, logn);
	for (u = 0; u < n; u ++) {
		s0[u] = fpr_rint(t0[u]);
	}

	/**
	 * Now run the casual verification.
	 */
	return Zf(complete_verify)(h0, h1, s0, s1, q00, q10, q11, bound, logn, tmp);
}
