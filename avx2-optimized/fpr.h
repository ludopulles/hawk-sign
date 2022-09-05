/*
 * Floating-point operations.
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


/* ====================================================================== */

#include <math.h>

/*
 * We wrap the native 'double' type into a structure so that the C compiler
 * complains if we inadvertently use raw arithmetic operators on the 'fpr'
 * type instead of using the inline functions below. This should have no
 * extra runtime cost, since all the functions below are 'inline'.
 */
typedef struct { double v; } fpr;

static inline fpr
FPR(double v)
{
	fpr x;

	x.v = v;
	return x;
}

static inline fpr
fpr_of(int64_t i)
{
	return FPR((double)i);
}

static const fpr fpr_zero = { 0.0 };
static const fpr fpr_one = { 1.0 };
static const fpr fpr_two = { 2.0 };
static const fpr fpr_onehalf = { 0.5 };
static const fpr fpr_ptwo31 = { 2147483648.0 };
static const fpr fpr_ptwo31m1 = { 2147483647.0 };
static const fpr fpr_mtwo31m1 = { -2147483647.0 };
static const fpr fpr_ptwo63m1 = { 9223372036854775807.0 };
static const fpr fpr_mtwo63m1 = { -9223372036854775807.0 };
static const fpr fpr_ptwo63 = { 9223372036854775808.0 };
static const fpr fpr_almost_one = { 0.98 };

static inline int64_t
fpr_rint(fpr x)
{
	/*
	 * We do not want to use llrint() since it might be not
	 * constant-time.
	 *
	 * Suppose that x >= 0. If x >= 2^52, then it is already an
	 * integer. Otherwise, if x < 2^52, then computing x+2^52 will
	 * yield a value that will be rounded to the nearest integer
	 * with exactly the right rules (round-to-nearest-even).
	 *
	 * In order to have constant-time processing, we must do the
	 * computation for both x >= 0 and x < 0 cases, and use a
	 * cast to an integer to access the sign and select the proper
	 * value. Such casts also allow us to find out if |x| < 2^52.
	 */
	int64_t sx, tx, rp, rn, m;
	uint32_t ub;

	sx = (int64_t)(x.v - 1.0);
	tx = (int64_t)x.v;
	rp = (int64_t)(x.v + 4503599627370496.0) - 4503599627370496;
	rn = (int64_t)(x.v - 4503599627370496.0) + 4503599627370496;

	/*
	 * If tx >= 2^52 or tx < -2^52, then result is tx.
	 * Otherwise, if sx >= 0, then result is rp.
	 * Otherwise, result is rn. We use the fact that when x is
	 * close to 0 (|x| <= 0.25) then both rp and rn are correct;
	 * and if x is not close to 0, then trunc(x-1.0) yields the
	 * appropriate sign.
	 */

	/*
	 * Clamp rp to zero if tx < 0.
	 * Clamp rn to zero if tx >= 0.
	 */
	m = sx >> 63;
	rn &= m;
	rp &= ~m;

	/*
	 * Get the 12 upper bits of tx; if they are not all zeros or
	 * all ones, then tx >= 2^52 or tx < -2^52, and we clamp both
	 * rp and rn to zero. Otherwise, we clamp tx to zero.
	 */
	ub = (uint32_t)((uint64_t)tx >> 52);
	m = -(int64_t)((((ub + 1) & 0xFFF) - 2) >> 31);
	rp &= m;
	rn &= m;
	tx &= ~m;

	/*
	 * Only one of tx, rn or rp (at most) can be non-zero at this
	 * point.
	 */
	return tx | rn | rp;
}

static inline int64_t
fpr_floor(fpr x)
{
	int64_t r;

	/*
	 * The cast performs a trunc() (rounding toward 0) and thus is
	 * wrong by 1 for most negative values. The correction below is
	 * constant-time as long as the compiler turns the
	 * floating-point conversion result into a 0/1 integer without a
	 * conditional branch or another non-constant-time construction.
	 * This should hold on all modern architectures with an FPU (and
	 * if it is false on a given arch, then chances are that the FPU
	 * itself is not constant-time, making the point moot).
	 */
	r = (int64_t)x.v;
	return r - (x.v < (double)r);
}

static inline int64_t
fpr_trunc(fpr x)
{
	return (int64_t)x.v;
}

static inline fpr
fpr_add(fpr x, fpr y)
{
	return FPR(x.v + y.v);
}

static inline fpr
fpr_sub(fpr x, fpr y)
{
	return FPR(x.v - y.v);
}

static inline fpr
fpr_neg(fpr x)
{
	return FPR(-x.v);
}

static inline fpr
fpr_half(fpr x)
{
	return FPR(x.v * 0.5);
}

static inline fpr
fpr_double(fpr x)
{
	return FPR(x.v + x.v);
}

static inline fpr
fpr_mul(fpr x, fpr y)
{
	return FPR(x.v * y.v);
}

static inline fpr
fpr_sqr(fpr x)
{
	return FPR(x.v * x.v);
}

static inline fpr
fpr_inv(fpr x)
{
	return FPR(1.0 / x.v);
}

static inline fpr
fpr_div(fpr x, fpr y)
{
	return FPR(x.v / y.v);
}

static inline int
fpr_lt(fpr x, fpr y)
{
	return x.v < y.v;
}

#define fpr_gm_tab   Zf(fpr_gm_tab)
extern const fpr fpr_gm_tab[];

#define fpr_p2_tab   Zf(fpr_p2_tab)
extern const fpr fpr_p2_tab[];

/* ====================================================================== */

