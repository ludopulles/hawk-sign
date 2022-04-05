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
 * The hash of a message has two parts: (h0, h1), where each is a polynomial of
 * length n with coefficients in {0,1}. However, these bits are collected into
 * bytes for saving RAM usage. This function returns the index in the hash
 * array at which the bits start for h1.
 */
#define SECOND_HASH(h, logn) \
	((h) + ((logn) <= 3 ? 1u : 1u << ((logn) - 3)))

/**
 * If s != NULL, set the polynomial p equal to h - 2 * s.
 * Otherwise, set p equal to h.
 * Returns the polynomial in FFT format.
 */
static void
hash_to_fft(fpr *p, const uint8_t *h, const int16_t *s, unsigned logn)
{
	size_t n, u, v;
	uint8_t hash;

	n = MKN(logn);
	if (logn <= 3) {
		hash = h[0];
		for (v = 0; v < n; v ++) {
			p[v] = fpr_of((hash & 1) - (s == NULL ? 0 : 2 * s[v]));
			hash >>= 1;
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				p[u] = fpr_of((hash & 1) - (s == NULL ? 0 : 2 * s[u]));
				hash >>= 1;
			}
		}
	}
	Zf(FFT)(p, logn);
}

/*
 * Returns whether squared geometric norm of (t0, t1) with respect to Q is at
 * most the specified allowed l2-bound. The t0, t1 are assumed to be in FFT
 * representation.
 */
static int
has_short_trace(fpr *t0, fpr *t1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	unsigned logn)
{
	size_t u, hn;
	fpr trace;

	hn = MKN(logn) >> 1;
	trace = fpr_zero;

	/*
	 * Calculate trace( (t0, t1)^* Q (t0, t1) ) / n, and determine if this is
	 * short. The trace of a polynomial in FFT representation is the sum of its
	 * `n` complex embeddings, and since half of the embeddings are stored, we
	 * only need to sum the real value for each pair of complex embeddings.
	 * 
	 * Add the contribution of t0 q00 t0*.
	 */
	for (u = 0; u < hn; u ++) {
		fpr norm_t0;

		norm_t0 = fpr_add(fpr_sqr(t0[u]), fpr_sqr(t0[u + hn]));
		trace = fpr_add(trace, fpr_mul(norm_t0, q00[u]));
	}

	/*
	 * Add the contribution of t1 q11 t1*.
	 */
	for (u = 0; u < hn; u ++) {
		fpr norm_t1;

		norm_t1 = fpr_add(fpr_sqr(t1[u]), fpr_sqr(t1[u + hn]));
		trace = fpr_add(trace, fpr_mul(norm_t1, q11[u]));
	}

	/*
	 * Add the contribution of t1 q10 t0* + t0 q01 t1*.
	 */
	for (u = 0; u < hn; u ++) {
		fpr re, im;

		// t0* t1
		re = fpr_add(fpr_mul(t0[u], t1[u]), fpr_mul(t0[u + hn], t1[u + hn]));
		im = fpr_sub(fpr_mul(t0[u], t1[u + hn]), fpr_mul(t0[u + hn], t1[u]));

		trace = fpr_add(trace, fpr_double(
			fpr_sub(fpr_mul(re, q10[u]), fpr_mul(im, q10[u + hn]))
		));
	}

	/*
	 * Renormalize to obtain the squared geometric norm of (t0, t1) w.r.t Q.
	 */
	trace = fpr_div(trace, fpr_of(hn));

	/*
	 * Check whether the norm is in range [0, 2^31). Signature is valid iff
	 * squared norm of (t0, t1) w.r.t. Q is at most bound.
	 */
	return !fpr_lt(trace, fpr_zero)
		&& fpr_lt(trace, fpr_ptwo31m1)
		&& (uint32_t)fpr_rint(trace) <= Zf(l2bound)[logn];
}

/* =================================================================== */

/* see inner.h */
void
Zf(complete_pubkey)(const int16_t *restrict iq00, const int16_t *restrict iq10,
	fpr *restrict q00, fpr *restrict q10, fpr *restrict q11, unsigned logn)
{
	size_t u, n, hn;

	n = MKN(logn);
	hn = n >> 1;

	// doing this in reverse, allows iq00, iq10 to overlap with the begin
	// of q00.
	for (u = n; u -- > 0; ) {
		q10[u] = fpr_of(iq10[u]);
	}
	for (u = n; u -- > 0; ) {
		q00[u] = fpr_of(iq00[u]);
	}

	Zf(FFT)(q00, logn);
	Zf(FFT)(q10, logn);

	/*
	 * Reconstruct q11 using q11 = (1 + q10 adj(q10)) / q00.
	 */
	Zf(poly_prod_selfadj_fft)(q11, q10, logn);
	for (u = 0; u < hn; u ++) {
		q11[u] = fpr_add(q11[u], fpr_one);
	}
	Zf(poly_div_autoadj_fft)(q11, q00, logn);
}

/* see inner.h */
int
Zf(complete_verify)(const uint8_t *restrict h,
	const int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n;
	fpr *t0, *t1;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;

	/*
	 * Put (t0, t1) = (h0 - 2 s0, h1 - 2 s1) in FFT representation.
	 */
	hash_to_fft(t0, h, s0, logn);
	hash_to_fft(t1, SECOND_HASH(h, logn), s1, logn);
	return has_short_trace(t0, t1, q00, q10, q11, logn);
}

/* see inner.h */
int
Zf(verify_simple_rounding)(const uint8_t *restrict h,
	int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u, v, w;
	uint8_t h0;
	fpr *t0;

	n = MKN(logn);
	t0 = (fpr *)tmp;

	hash_to_fft(t0, SECOND_HASH(h, logn), s1, logn);
	Zf(poly_mul_fft)(t0, q10, logn);
	Zf(poly_div_autoadj_fft)(t0, q00, logn);
	Zf(iFFT)(t0, logn);

	/*
	 * Recover s0 with s0 = round(h0 / 2 + (h1 / 2 - s1) q10 / q00).
	 * Put (t0, t1) = (h0 - 2 * s0, h1 - 2 * s1) in FFT representation.
	 */
	if (logn <= 3) {
		h0 = h[0];
		for (v = 0; v < n; v ++) {
			s0[v] = fpr_rint(fpr_half(fpr_add(fpr_of(h0 & 1), t0[v])));
			h0 >>= 1;
		}
	} else {
		for (u = 0, w = 0; w < n; u ++) {
			h0 = h[u];
			for (v = 0; v < 8; v ++, w ++) {
				s0[w] = fpr_rint(fpr_half(fpr_add(fpr_of(h0 & 1), t0[w])));
				h0 >>= 1;
			}
		}
	}

	return Zf(complete_verify)(h, s0, s1, q00, q10, q11, logn, tmp);
}

/* see inner.h */
int
Zf(verify_nearest_plane)(const uint8_t *restrict h,
	int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	unsigned logn, uint8_t *restrict tmp)
{
	/*
	 * This works slightly better than simple rounding, but is also slower.
	 * Reconstruct s0, by running Babai's NP algorithm with target
	 *     h0 / 2 + (h1 / 2 - s1) * q10 / q00.
	 */

	size_t n;
	fpr *t0, *t1;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;

	hash_to_fft(t0, h, NULL, logn);
	hash_to_fft(t1, SECOND_HASH(h, logn), s1, logn);

	Zf(poly_mul_fft)(t1, q10, logn);
	Zf(poly_div_autoadj_fft)(t1, q00, logn);
	Zf(poly_add)(t0, t1, logn);
	Zf(poly_mulconst)(t0, fpr_onehalf, logn);

	memcpy(t1, q00, n * sizeof(fpr));
	/*
	 * Run Babai with target t0 and Gram-matrix q00.
	 */
	Zf(ffNearestPlane_dyn)(t0, t1, logn, t1 + n);
	Zf(fft_to_int16)(s0, t0, logn);

	/*
	 * Now run the casual verification.
	 */
	return Zf(complete_verify)(h, s0, s1, q00, q10, q11, logn, tmp);
}

/* see inner.h */
int
Zf(verify_simple_rounding_fft)(const uint8_t *restrict h,
	const int16_t *restrict s1, const fpr *restrict q00,
	const fpr *restrict q10, const fpr *restrict q11,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t u, v, w, n;
	uint8_t h0;
	int16_t s0w;
	fpr *t0, *t1;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;

	hash_to_fft(t1, SECOND_HASH(h, logn), s1, logn);
	Zf(poly_prod_fft)(t0, t1, q10, logn);
	Zf(poly_div_autoadj_fft)(t0, q00, logn);
	Zf(iFFT)(t0, logn);

	/*
	 * Recover s0 with s0 = round(h0 / 2 + (h1 / 2 - s1) q10 / q00).
	 * Put (t0, t1) = (h0 - 2 * s0, h1 - 2 * s1) in FFT representation.
	 */
	if (logn <= 3) {
		h0 = h[0];
		for (v = 0; v < n; v ++) {
			s0w = fpr_rint(fpr_half(fpr_add(fpr_of(h0 & 1), t0[v])));
			t0[v] = fpr_of((h0 & 1) - 2 * s0w);
			h0 >>= 1;
		}
	} else {
		for (u = 0, w = 0; w < n; u ++) {
			h0 = h[u];
			for (v = 0; v < 8; v ++, w ++) {
				s0w = fpr_rint(fpr_half(fpr_add(fpr_of(h0 & 1), t0[w])));
				t0[w] = fpr_of((h0 & 1) - 2 * s0w);
				h0 >>= 1;
			}
		}
	}
	Zf(FFT)(t0, logn);

	return has_short_trace(t0, t1, q00, q10, q11, logn);
}

