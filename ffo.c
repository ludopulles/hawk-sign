/*
 * Hawk ffNearestPlane implementation to be used in key generation.
 * ffLDL-tree is taken from Falcon's sign.c code.
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2022 Hawk Project
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
 * @author   Ludo Pulles <ludo.pulles@cwi.nl>
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
