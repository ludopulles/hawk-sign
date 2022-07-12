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
 * @author   Ludo Pulles <ludo.pulles@cwi.nl>
 */

#include "inner.h"

/* =================================================================== */
/*
 * Binary case:
 *   N = 2^logn
 *   phi = X^N+1
 */

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
	tree1 = tree + n + LDL_TREESIZE(logn - 1);

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
 * Perform Fast Fourier Nearest Plane for target vector t and LDL tree T.
 * tmp[] must have size for at least two polynomials of size 2^logn.
 */
static void
ffNearestPlane_inner(fpr *restrict z0, fpr *restrict z1,
	const fpr *restrict tree, const fpr *restrict t0,
	const fpr *restrict t1, unsigned logn, fpr *restrict tmp)
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
	tree1 = tree + n + LDL_TREESIZE(logn - 1);

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

/* see inner.h */
void
Zf(ffNearestPlane_tree)(fpr *restrict z, const fpr *restrict tree,
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

/* see inner.h */
void
Zf(ffNearestPlane_dyn)(fpr *restrict t, fpr *restrict g, unsigned logn, fpr *restrict tmp)
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
	Zf(ffNearestPlane_dyn)(g + hn, t + hn, logn - 1, tmp + n);
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
	Zf(ffNearestPlane_dyn)(g, tmp, logn - 1, tmp + hn);
	// Memory layout: t: L_10, ???; g: z_0, z_1; tmp: ???

	Zf(poly_merge_fft)(t, g, g + hn, logn);
	// Memory layout: t: z; g: z_0, z_1; tmp: ???
}

/* see inner.h */
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
	for (u = 0; u < hn; u++) {
		q[u] = fpr_inv(q[u]);
		q[u + hn] = fpr_zero;
	}
	// Alternative:
	// Zf(poly_add_muladj_fft)(q, f, g, f, g, logn);

	/*
	 * Now execute Babai with target t and Gram matrix q, where
	 *     t = (F adj(f) + G adj(g)) / (f adj(f) + g adj(g)),
	 *     q = f adj(f) + g adj(g).
	 */
	Zf(ffNearestPlane_dyn)(t, q, logn, tmp + 2*n);

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
}
