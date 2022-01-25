/*
 * ffLDL-tree taken from Falcon's sign.c code.
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

/* =================================================================== */
/*
 * Binary case:
 *   N = 2^logn
 *   phi = X^N+1
 */

/*
 * Get the size of the LDL tree for an input with polynomials of size
 * 2^logn. The size is expressed in the number of elements.
 */
static inline unsigned
ffLDL_treesize(unsigned logn)
{
	/*
	 * For logn = 0 (polynomials are constant), the tree is empty as the leaves
	 * do not require any information to be stored. Otherwise, the tree node
	 * has size 2^logn, and has two child trees for size logn-1 each. Thus,
	 * treesize s() must fulfill these two relations:
	 *
	 *   s(0) = s(1) = 0
	 *   s(logn) = (2^logn) + 2*s(logn-1)
	 */
	return (logn - 1) << logn;
}

/*
 * Inner function for ffLDL_fft(). It expects the matrix to be both
 * auto-adjoint and quasicyclic; also, it uses the source operands
 * as modifiable temporaries.
 *
 * tmp[] must have room for at least one polynomial.
 */
static void
ffLDL_fft_inner(fpr *restrict tree,
	fpr *restrict g0, fpr *restrict g1, unsigned logn, fpr *restrict tmp)
{
	size_t n, hn;
	fpr *tree0, *tree1;

	if (logn == 0)
		return;

	n = MKN(logn);
	hn = n >> 1;
	tree0 = tree + n;
	tree1 = tree + n + ffLDL_treesize(logn - 1);

	/*
	 * The LDL decomposition yields L (which is written in the tree)
	 * and the diagonal of D. Since d00 = g0, we just write d11
	 * into tmp.
	 */
	Zf(poly_LDLmv_fft)(tmp, tree, g0, g1, g0, logn);

	/*
	 * Split d00 (currently in g0) and d11 (currently in tmp). We
	 * reuse g0 and g1 as temporary storage spaces:
	 *   d00 splits into g1, g1+hn
	 *   d11 splits into g0, g0+hn
	 */
	Zf(poly_split_fft)(g1, g1 + hn, g0, logn);
	Zf(poly_split_fft)(g0, g0 + hn, tmp, logn);

	/*
	 * Each split result is the first row of a new auto-adjoint
	 * quasicyclic matrix for the next recursive step.
	 */
	ffLDL_fft_inner(tree0, g1, g1 + hn, logn - 1, tmp);
	ffLDL_fft_inner(tree1, g0, g0 + hn, logn - 1, tmp);
}

/* see inner.h */
void
Zf(ffLDL_fft)(fpr *restrict tree, const fpr *restrict q00, unsigned logn,
	fpr *restrict tmp)
{
	size_t n, hn;
	fpr *g0, *g1;

	if (logn <= 1)
		return;

	n = MKN(logn);
	hn = n >> 1;
	g0 = tmp;
	g1 = tmp + hn;
	tmp += n;

	Zf(poly_split_fft)(g0, g1, q00, logn);
	ffLDL_fft_inner(tree, g0, g1, logn - 1, tmp);
}

/*
 * Perform Fast Fourier Sampling for target vector t and LDL tree T.
 * tmp[] must have size for at least two polynomials of size 2^logn.
 */
static void
ffNearestPlane_inner(fpr *restrict z0, fpr *restrict z1,
	const fpr *restrict tree, const fpr *restrict t0, const fpr *restrict t1,
	unsigned logn, fpr *restrict tmp)
{
	size_t n, hn;
	const fpr *tree0, *tree1;

	/*
	 * Normal end of recursion is for logn == 0.
	 * Inline the last two recursion levels to get better performance here.
	 */
	if (logn == 0) {
		/*
		 * Simple rounding when n=1, since Z[\zeta] = Z[i] has two orthogonal
		 * basis vectors so Babai == simple rounding in this case.
		 */
		z0[0] = fpr_of(fpr_rint(t0[0]));
		z1[0] = fpr_of(fpr_rint(t1[0]));
		return;
	}

	/*
	 * General recursive case (logn >= 1).
	 */
	n = MKN(logn);
	hn = n >> 1;
	tree0 = tree + n;
	tree1 = tree + n + ffLDL_treesize(logn - 1);

	/*
	 * We split t1 into z1 (reused as temporary storage), then do
	 * the recursive invocation, with output in tmp. We finally
	 * merge back into z1.
	 */
	Zf(poly_split_fft)(z1, z1 + hn, t1, logn);
	ffNearestPlane_inner(tmp, tmp + hn,
		tree1, z1, z1 + hn, logn - 1, tmp + n);
	Zf(poly_merge_fft)(z1, tmp, tmp + hn, logn);

	/*
	 * Compute tb0 = t0 + (t1 - z1) * L. Value tb0 ends up in tmp[].
	 */
	memcpy(tmp, t1, n * sizeof *t1);
	Zf(poly_sub)(tmp, z1, logn);
	Zf(poly_mul_fft)(tmp, tree, logn);
	Zf(poly_add)(tmp, t0, logn);

	/*
	 * Second recursive invocation.
	 */
	Zf(poly_split_fft)(z0, z0 + hn, tmp, logn);
	ffNearestPlane_inner(tmp, tmp + hn,
		tree0, z0, z0 + hn, logn - 1, tmp + n);
	Zf(poly_merge_fft)(z0, tmp, tmp + hn, logn);
}

/*
 * Perform Fast Fourier Sampling for target vector t and LDL tree T.
 * tmp[] must have size for at least two polynomials of size 2^logn.
 */
void
Zf(ffNearestPlane)(fpr *restrict z, const fpr *restrict tree,
	const fpr *restrict t, unsigned logn, fpr *restrict tmp)
{
	size_t n, hn;
	fpr *t0, *t1;

	n = MKN(logn);
	hn = n >> 1;

	t0 = tmp;
	t1 = t0 + hn;
	tmp += n;

	Zf(poly_split_fft)(t0, t1, t, logn);
	ffNearestPlane_inner(z, z + hn, tree, t0, t1, logn - 1, tmp);
	memcpy(t0, z, n * sizeof(fpr));
	Zf(poly_merge_fft)(z, t0, t1, logn);
}

/*
 * Convert an integer polynomial (with small values) into the
 * representation with complex numbers.
 */
static void
smallints_to_fpr(fpr *r, const int8_t *t, unsigned logn)
{
	size_t n = MKN(logn), u;
	for (u = 0; u < n; u ++) {
		r[u] = fpr_of(t[u]);
	}
}


void
Zf(ffBabai)(const int8_t *restrict f, const int8_t *restrict g,
	int8_t *restrict F, int8_t *restrict G, unsigned logn,
	fpr *restrict tree, fpr *restrict tmp)
{
	size_t n = MKN(logn);

	fpr *t = tmp, *z = t + n;
	fpr *_f = z + n, *_g = _f + n;
	fpr *_F = _g + n, *_G = _F + n;
	tmp += 6*n;

	smallints_to_fpr(_f, f, logn);
	smallints_to_fpr(_g, g, logn);
	smallints_to_fpr(_F, F, logn);
	smallints_to_fpr(_G, G, logn);
	Zf(FFT)(_f, logn);
	Zf(FFT)(_g, logn);
	Zf(FFT)(_F, logn);
	Zf(FFT)(_G, logn);

	// calculate (F adj(f) + G adj(g)) / (f adj(f) + g adj(g))
	Zf(poly_add_muladj_fft)(t, _F, _G, _f, _g, logn);
	Zf(poly_invnorm2_fft)(z, _f, _g, logn);
	Zf(poly_mul_autoadj_fft)(t, z, logn);

	// inverting an autoadj polynomial:
	for (size_t u = 0; u < n/2; u++)
		z[u] = fpr_inv(z[u]),
		z[u + n/2] = fpr_zero;

	Zf(ffLDL_fft)(tree, z, logn, tmp);
	Zf(ffNearestPlane)(z, tree, t, logn, tmp);

	/* memcpy(_F, z, n * sizeof(fpr));
	Zf(iFFT)(_F, logn);
	printf("z = ");
	for (size_t u = 0; u < n; u++) printf("%d ", fpr_rint(_F[u]));
	printf("\n"); */

	Zf(poly_mul_fft)(_f, z, logn);
	Zf(poly_mul_fft)(_g, z, logn);
	Zf(iFFT)(_f, logn);
	Zf(iFFT)(_g, logn);
	// subtract (f,g) from (F,G)

	for (size_t u = 0; u < n; u++) {
		F[u] -= fpr_rint(_f[u]);
		G[u] -= fpr_rint(_g[u]);
	}
}

/* =======================================================================
 * Do straight Babai, if you build the tree once for every Babai call,
 * so building tree is not more efficient.
 * tmp[] requires a size of at least 6 2^logn fpr values.
 */
// TODO: make this method use less memory, the result (z0, z1) can
// probably be stored in (t0, t1).
static void
ffBabai_inner(
	fpr *restrict z0, fpr *restrict z1,
	fpr *restrict g0, fpr *restrict g1,
	fpr *restrict t0, fpr *restrict t1,
	unsigned logn, fpr *restrict tmp)
{
	size_t n, hn;
	fpr *d11, *l10;

	/*
	 * Normal end of recursion is for logn == 0.
	 * Inline the last two recursion levels to get better performance here.
	 */
	if (logn == 0) {
		/*
		 * Simple rounding when n=1, since Z[\zeta] = Z[i] has two orthogonal
		 * basis vectors so Babai == simple rounding in this case.
		 */
		z0[0] = fpr_of(fpr_rint(t0[0]));
		z1[0] = fpr_of(fpr_rint(t1[0]));
		return;
	}

	n = MKN(logn);
	hn = n >> 1;

	d11 = tmp;
	l10 = d11 + n;
	tmp += 2*n;

	/*
	 * The LDL decomposition yields L and the diagonal of D. Note that d00 = g0.
	 */
	Zf(poly_LDLmv_fft)(d11, l10, g0, g1, g0, logn);

	/*
	 * Split d00 (currently in g0) and d11 (currently in tmp). We
	 * reuse g0 and g1 as temporary storage spaces:
	 *   d00 splits into g1, g1+hn
	 *   d11 splits into g0, g0+hn
	 */
	Zf(poly_split_fft)(g1, g1 + hn, g0, logn);
	Zf(poly_split_fft)(g0, g0 + hn, d11, logn);

	/*
	 * Each split result is the first row of a new auto-adjoint
	 * quasicyclic matrix for the next recursive step.
	 * We split t1 into z1 (reused as temporary storage), then do
	 * the recursive invocation, with output in tmp. We finally
	 * merge back into z1.
	 */
	Zf(poly_split_fft)(z1, z1 + hn, t1, logn);
	ffBabai_inner(tmp, tmp + hn, g0, g0 + hn, z1, z1 + hn, logn - 1, tmp + n);
	Zf(poly_merge_fft)(z1, tmp, tmp + hn, logn);

	/*
	 * Compute tb0 = t0 + (t1 - z1) * L. Value tb0 ends up in tmp[].
	 */
	memcpy(tmp, t1, n * sizeof *t1);
	Zf(poly_sub)(tmp, z1, logn);
	Zf(poly_mul_fft)(tmp, l10, logn);
	Zf(poly_add)(tmp, t0, logn);

	/*
	 * Second recursive invocation.
	 */
	Zf(poly_split_fft)(z0, z0 + hn, tmp, logn);
	ffBabai_inner(tmp, tmp + hn, g1, g1 + hn, z0, z0 + hn, logn - 1, tmp + n);
	Zf(poly_merge_fft)(z0, tmp, tmp + hn, logn);
}

// tmp[] requires a size of at least 6 2^logn fpr values.
void
Zf(ffBabai_reduce)(const fpr *restrict f, const fpr *restrict g,
	fpr *restrict F, fpr *restrict G, int8_t *restrict Fn,
	int8_t *restrict Gn, unsigned logn, fpr *tmp)
{
	size_t n, hn;
	fpr *t, *z, *g0, *g1;

	n = MKN(logn);
	hn = n >> 1;
	t = tmp;
	z = t + n;
	g0 = z + n;
	g1 = g0 + hn;

	Zf(poly_add_muladj_fft)(t, F, G, f, g, logn);
	Zf(poly_invnorm2_fft)(z, f, g, logn);
	Zf(poly_mul_autoadj_fft)(t, z, logn);

	// inverting an autoadj polynomial:
	for (size_t u = 0; u < hn; u++)
		z[u] = fpr_inv(z[u]);
	for (size_t u = hn; u < n; u++)
		z[u] = fpr_zero;

	// t = (F adj(f) + G adj(g)) / (f adj(f) + g adj(g))
	// z = f adj(f) + g adj(g)

	// Execute Babai with Gram-matrix z, target t, and put the result in z.
	Zf(poly_split_fft)(g0, g1, z, logn);

	// currently, t is the target
	Zf(poly_split_fft)(z, z + hn, t, logn);
	ffBabai_inner(t, t + hn, g0, g1, z, z + hn, logn - 1, tmp + 3*n);
	Zf(poly_merge_fft)(z, t, t + hn, logn);
	// currently, z is the result

	// modify F and Fn:
	memcpy(t, z, n * sizeof *z);
	Zf(poly_mul_fft)(t, f, logn);
	Zf(poly_sub)(F, t, logn); // F -= z*f
	Zf(iFFT)(t, logn);
	for (size_t u = 0; u < n; u++) 
		Fn[u] -= fpr_rint(t[u]);

	// modify G and Gn:
	memcpy(t, z, n * sizeof *z);
	Zf(poly_mul_fft)(t, g, logn);
	Zf(poly_sub)(G, t, logn); // G -= z*f
	Zf(iFFT)(t, logn);
	for (size_t u = 0; u < n; u++) 
		Gn[u] -= fpr_rint(t[u]);

/*
	// If we want to print what we subtracted:
	Zf(iFFT)(z, logn);
	printf("z = ");
	for (size_t u = 0; u < n; u++) printf("%d ", fpr_rint(z[u]));
	printf("\n"); */
}

// tmp[] requires a size of at least 10 2^logn fpr values.
void
Zf(ffStraightBabai)(const int8_t *restrict f, const int8_t *restrict g,
	int8_t *restrict F, int8_t *restrict G,
	unsigned logn, fpr *restrict tmp)
{
	size_t n;
	fpr *_f, *_g, *_F, *_G;

	n = MKN(logn);

	_f = tmp;
	_g = _f + n;
	_F = _g + n;
	_G = _F + n;

	smallints_to_fpr(_f, f, logn);
	smallints_to_fpr(_g, g, logn);
	smallints_to_fpr(_F, F, logn);
	smallints_to_fpr(_G, G, logn);

	Zf(FFT)(_f, logn);
	Zf(FFT)(_g, logn);
	Zf(FFT)(_F, logn);
	Zf(FFT)(_G, logn);

	Zf(ffBabai_reduce)(_f, _g, _F, _G, F, G, logn, tmp + 4*n);
}

// TODO: test this code
