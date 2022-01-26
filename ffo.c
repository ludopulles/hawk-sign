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
 * Run Babai's nearest plane algorithm directly with the Gram matrix q00
 * of a single basis element (f,g).
 * Assumes t, g are both of length 2^logn
 *
 * tmp[] requires a size of at least 2 2^logn fpr values.
 */
static void
ffBabai_inner(fpr *restrict t, fpr *restrict g, unsigned logn, fpr *restrict tmp)
{
	size_t n, hn;

	/*
	 * Normal end of recursion is for logn == 1.
	 * Inline the last two recursion levels to get better performance here.
	 */
	if (logn == 1) {
		/*
		 * Simple rounding when n=2, since Z[\zeta] = Z[i] has two orthogonal
		 * basis vectors so Babai == simple rounding in this case.
		 */
		t[0] = fpr_of(fpr_rint(t[0]));
		t[1] = fpr_of(fpr_rint(t[1]));
		return;
	}

	n = MKN(logn);
	hn = n >> 1;

	Zf(poly_split_fft)(tmp, tmp + hn, g, logn);
	Zf(poly_split_fft)(g, g + hn, t, logn);
	// Memory layout: t: t; g: t_0, t_1; tmp: g_0, g_1.

	/*
	 * The LDL decomposition yields L and the diagonal of D. Note that D_00 = g_0.
	 */
	Zf(poly_LDLmv_fft)(t + hn, t, tmp, tmp + hn, tmp, logn - 1);
	// Memory layout: t: L_10, D_11; g: t_0, t_1; tmp: g_0 = D_00, g_1.

	memcpy(tmp + hn, g + hn, hn * sizeof *t);
	// Memory layout: t: L_10, D_11; g: t_0, t_1; tmp: g_0 = D_00, t_1.

	/*
	 * First recursive invocation: target is t_1 and gram matrix is D_11.
	 */
	ffBabai_inner(g + hn, t + hn, logn - 1, tmp + n);
	// Memory layout: t: L_10, ???; g: t_0, z_1; tmp: g_0 = D_00, t_1.

	/*
	 * Compute t'_0 = t_0 + (t_1 - z_1) * L_{10}, and put this value in t_0
	 */
	Zf(poly_sub)(tmp + hn, g + hn, logn - 1);
	Zf(poly_mul_fft)(tmp + hn, t, logn - 1);
	Zf(poly_add)(g, tmp + hn, logn - 1);
	// Memory layout: t: L_10, ???; g: t'_0, z_1; tmp: g_0 = D_00, (t_1-z_1)*L_10

	/*
	 * Second recursive invocation: target is t'_0 and gram matrix is D_00.
	 */
	ffBabai_inner(g, tmp, logn - 1, tmp + n);
	// Memory layout: t: L_10, ???; g: z_0, z_1; tmp: ???

	Zf(poly_merge_fft)(t, g, g + hn, logn);
	// Memory layout: t: z; g: z_0, z_1; tmp: ???

	// ffBabai_inner(z, g1, g1 + hn, tmp, logn - 1, tmp + n);
	// Zf(poly_split_fft)(z0, z0 + hn, tmp, logn);
	// ffBabai_inner(tmp, tmp + hn, g1, g1 + hn, z0, z0 + hn, logn - 1, tmp + n);
	// Zf(poly_merge_fft)(z0, tmp, tmp + hn, logn);
}

// tmp[] requires a size of at least 4 2^logn fpr values, so 32 2^logn bytes.
void
Zf(ffBabai_reduce)(const fpr *restrict f, const fpr *restrict g,
	fpr *restrict F, fpr *restrict G, int8_t *restrict Fn,
	int8_t *restrict Gn, unsigned logn, fpr *tmp)
{
	size_t n, hn, u;
	fpr *t, *q;

	n = MKN(logn);
	hn = n >> 1;
	t = tmp;
	q = t + n;

	Zf(poly_add_muladj_fft)(t, F, G, f, g, logn);
	Zf(poly_invnorm2_fft)(q, f, g, logn);
	Zf(poly_mul_autoadj_fft)(t, q, logn);

	// Since q is self-adjoint, invert real part
	for (u = 0; u < hn; u++)
		q[u] = fpr_inv(q[u]);
	// Alternative:
	// Zf(poly_add_muladj_fft)(q, f, g, f, g, logn);

	/**
	 * Now we have the target t and Gram matrix q:
	 * t = (F adj(f) + G adj(g)) / (f adj(f) + g adj(g))
	 * q = f adj(f) + g adj(g)
	 */

	// Execute Babai with Gram-matrix q, target t, and put the result in t.
	ffBabai_inner(t, q, logn, tmp + 2*n);

	// modify F and Fn:
	memcpy(q, t, n * sizeof *t);
	Zf(poly_mul_fft)(q, f, logn);
	Zf(poly_sub)(F, q, logn); // F -= k*f
	Zf(iFFT)(q, logn);
	for (u = 0; u < n; u++)
		Fn[u] -= fpr_rint(q[u]);

	// modify G and Gn:
	memcpy(q, t, n * sizeof *t);
	Zf(poly_mul_fft)(q, g, logn);
	Zf(poly_sub)(G, q, logn); // G -= k*f
	Zf(iFFT)(q, logn);
	for (u = 0; u < n; u++)
		Gn[u] -= fpr_rint(q[u]);

/*	// If we want to print what we subtracted:
	Zf(iFFT)(t, logn);
	printf("z = ");
	for (u = 0; u < n; u++) printf("%d ", fpr_rint(t[u]));
	printf("\n"); */
}

// tmp[] requires a size of at least 8 2^logn fpr values, so 64 2^logn bytes.
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
