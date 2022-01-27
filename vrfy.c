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

int
Zf(verify)(const int8_t *restrict hm,
	int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp)
{
	size_t u, n;
	fpr *t0, *t1, *t2, *t3, trace;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;
	t2 = t1 + n;
	t3 = t2 + n;

	// multiply s1 by 2
	for (u = 0; u < n; u ++) {
		t0[u] = fpr_of(2 * s1[u]);
	}
	for (u = 0; u < n; u ++) {
		t2[u] = fpr_of(hm[u] & 1);
	}

	// Compute s0 = h%2 + 2 round(-q10 s1 / (2 q00) - h%2)
	Zf(FFT)(t0, logn);
	Zf(FFT)(t2, logn);
	// copy s1 for later.
	memcpy(t1, t0, n * sizeof *t0);

	// Zf(poly_mulconst)(t0, fpr_onehalf, logn);
	Zf(poly_neg)(t0, logn);
	Zf(poly_mul_fft)(t0, q10, logn); // -q10 s1
	// Note: q00 is self adjoint
	Zf(poly_div_autoadj_fft)(t0, q00, logn); // -s1 q10/q00
	Zf(poly_sub)(t0, t2, logn); // -s1 q10/q00 - h%2
	Zf(iFFT)(t0, logn);

	for (u = 0; u < n; u ++) {
		s0[u] = fpr_rint(fpr_half(t0[u]));
	}

	// Currently in memory: s0, s1 (in FFT representation)
	for (u = 0; u < n; u ++) {
		t0[u] = fpr_of(2 * s0[u] + (hm[u] & 1));
	}
	Zf(FFT)(t0, logn);

	// Currently in memory: s0, s1, s1, s0 (in FFT representation)
	memcpy(t2, t1, n * sizeof *t0);
	memcpy(t3, t0, n * sizeof *t0);

	// Compute s0 q00 s0* + s0 q01 s1* + s1 q10 s0* + s1 q11 s1*
	Zf(poly_mulselfadj_fft)(t2, logn);
	Zf(poly_mulselfadj_fft)(t3, logn);
	Zf(poly_mul_autoadj_fft)(t2, q11, logn); // t2 = s1 q11 s1*
	Zf(poly_mul_autoadj_fft)(t3, q00, logn); // t3 = s0 q00 s0*
	Zf(poly_muladj_fft)(t1, t0, logn); // t1 = s1 s0*
	Zf(poly_mul_fft)(t1, q10, logn); // t1 = s1 q10 s0*

	Zf(poly_addselfadj_fft)(t1, logn); // t1 = s1 q10 s0* + s0 q01 s1*
	Zf(poly_add_autoadj_fft)(t1, t2, logn);
	Zf(poly_add_autoadj_fft)(t1, t3, logn);

	trace = fpr_zero;
	for (u = 0; u < n/2; u ++) {
		trace = fpr_add(trace, t1[u]);
	}

	// note: only n/2 embeddings are stored,
	// the others are simply the conjugate embeddings.
	// TODO: this can be optimized in the verif_bound, cancelling with 2 in (2d).
	trace = fpr_double(trace);

	/*
	 * Signature is valid if and only if `v` is short enough.
	 */
	return fpr_lt(trace, fpr_of(bound * n));
}


// Verifies a signature given by (s0, s1) instead of only s1, so it does not
// have to reconstruct s0.
int
Zf(complete_verify)(const int8_t *restrict hm,
	const int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp)
{
	size_t u, n;
	fpr *t0, *t1, *t2, *t3, trace;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;
	t2 = t1 + n;
	t3 = t2 + n;

	for (u = 0; u < n; u ++) {
		t0[u] = fpr_of(2 * s0[u] + (hm[u] & 1));
	}
	for (u = 0; u < n; u ++) {
		t1[u] = fpr_of(2 * s1[u]);
	}

	Zf(FFT)(t0, logn);
	Zf(FFT)(t1, logn);

	// Currently in memory: s0, s1, s1, s0 (in FFT representation)
	memcpy(t2, t1, n * sizeof *t0);
	memcpy(t3, t0, n * sizeof *t0);

	// Compute s0 q00 s0* + s0 q01 s1* + s1 q10 s0* + s1 q11 s1*
	Zf(poly_mulselfadj_fft)(t2, logn);
	Zf(poly_mulselfadj_fft)(t3, logn);
	Zf(poly_mul_autoadj_fft)(t2, q11, logn); // t2 = s1 q11 s1*
	Zf(poly_mul_autoadj_fft)(t3, q00, logn); // t3 = s0 q00 s0*
	Zf(poly_muladj_fft)(t1, t0, logn); // t1 = s1 s0*
	Zf(poly_mul_fft)(t1, q10, logn); // t1 = s1 q10 s0*

	Zf(poly_addselfadj_fft)(t1, logn); // t1 = s1 q10 s0* + s0 q01 s1*
	Zf(poly_add_autoadj_fft)(t1, t2, logn);
	Zf(poly_add_autoadj_fft)(t1, t3, logn);

	trace = fpr_zero;
	for (u = 0; u < n/2; u ++) {
		trace = fpr_add(trace, t1[u]);
	}

	// note: only n/2 embeddings are stored,
	// the others are simply the conjugate embeddings.
	trace = fpr_double(trace);

	/*
	 * Signature is valid if and only if `Tr(s* Q s) / n = Tr(x* x)/ <= bound`.
	 */
	return fpr_lt(trace, fpr_of(bound * n));
}


