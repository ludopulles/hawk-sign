#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#include "../keygen.c"

/*
Table of average number of bits required to represent all the coefficients of
the polynomials f,g taken over 1000 samples:

Depth    avg  ( std) --> (x ints)
--------------------------------------------------------------------------------
Depth 0: 3.00 (0.00) --> (1 ints)
Depth 1: 8.05 (0.22) --> (1 ints)
Depth 2: 18.62 (0.49) --> (1 ints)
Depth 3: 38.67 (0.65) --> (2 ints)
Depth 4: 77.61 (1.31) --> (3 ints)
Depth 5: 153.50 (2.31) --> (6 ints)
Depth 6: 302.62 (3.95) --> (11 ints)
Depth 7: 597.95 (6.94) --> (21 ints)
Depth 8: 1184.35 (12.00) --> (41 ints)
Depth 9: 2368.50 (23.99) --> (82 ints)
*/

// Set large enough upper bound:
static const size_t __MAX_BL_LARGE__[9] = { 2, 2, 6, 8, 16, 32, 64, 128, 200 };

/*
 * Helper function to measure number of bits required to represent a number.
 */
static size_t
number_of_bits(uint32_t *fp, size_t sz)
{
	size_t sgn = (fp[sz-1] >> 30) & 1; // 0 = positive, 1 = negative
	while (sz > 0 && (fp[sz-1] == 0 || fp[sz-1] == 2147483647U)) sz--;

	if (sz == 0) return 1 + sgn; // 0 or -1 as value

	size_t res = (sz-1) * 31;
	for (size_t b = 31; b -- > 0; )
		if (((fp[sz-1] >> b) & 1) != sgn)
			return res + (b+1);
	fprintf(stderr, "Unexpected value: %d\n", (int) fp[sz-1]);
	assert(0);
}

static int
_solve_NTRU_intermediate(unsigned logn_top,
	const int8_t *f, const int8_t *g, size_t *unreduced_bits, unsigned depth,
	uint32_t *tmp)
{
	/*
	 * In this function, 'logn' is the log2 of the degree for
	 * this step. If N = 2^logn, then:
	 *  - the F and G values already in fk->tmp (from the deeper
	 *    levels) have degree N/2;
	 *  - this function should return F and G of degree N.
	 */
	unsigned logn;
	size_t n, hn, slen, dlen, llen, rlen, FGlen, u;
	uint32_t *Fd, *Gd, *Ft, *Gt, *ft, *gt, *t1;
	fpr *rt1, *rt2, *rt3, *rt4, *rt5;
	int scale_fg, minbl_fg, maxbl_fg, maxbl_FG, scale_k;
	uint32_t *x, *y;
	int32_t *k;
	const small_prime *primes;

	logn = logn_top - depth;
	n = (size_t)1 << logn;
	hn = n >> 1;

	/*
	 * slen = size for our input f and g; also size of the reduced
	 *        F and G we return (degree N)
	 *
	 * dlen = size of the F and G obtained from the deeper level
	 *        (degree N/2 or N/3)
	 *
	 * llen = size for intermediary F and G before reduction (degree N)
	 *
	 * We build our non-reduced F and G as two independent halves each,
	 * of degree N/2 (F = F0 + X*F1, G = G0 + X*G1).
	 */
	slen = MAX_BL_SMALL[depth];
	dlen = MAX_BL_SMALL[depth + 1];
	llen = __MAX_BL_LARGE__[depth];
	primes = PRIMES;

	/*
	 * Fd and Gd are the F and G from the deeper level.
	 */
	Fd = tmp;
	Gd = Fd + dlen * hn;

	/*
	 * Compute the input f and g for this level. Note that we get f
	 * and g in RNS + NTT representation.
	 */
	ft = Gd + dlen * hn;
	make_fg(ft, f, g, logn_top, depth, 1);

	/*
	 * Move the newly computed f and g to make room for our candidate
	 * F and G (unreduced).
	 */
	Ft = tmp;
	Gt = Ft + n * llen;
	t1 = Gt + n * llen;
	memmove(t1, ft, 2 * n * slen * sizeof *ft);
	ft = t1;
	gt = ft + slen * n;
	t1 = gt + slen * n;

	/*
	 * Move Fd and Gd _after_ f and g.
	 */
	memmove(t1, Fd, 2 * hn * dlen * sizeof *Fd);
	Fd = t1;
	Gd = Fd + hn * dlen;

	/*
	 * We reduce Fd and Gd modulo all the small primes we will need,
	 * and store the values in Ft and Gt (only n/2 values in each).
	 */
	for (u = 0; u < llen; u ++) {
		uint32_t p, p0i, R2, Rx;
		size_t v;
		uint32_t *xs, *ys, *xd, *yd;

		p = primes[u].p;
		p0i = modp_ninv31(p);
		R2 = modp_R2(p, p0i);
		Rx = modp_Rx((unsigned)dlen, p, p0i, R2);
		for (v = 0, xs = Fd, ys = Gd, xd = Ft + u, yd = Gt + u;
			v < hn;
			v ++, xs += dlen, ys += dlen, xd += llen, yd += llen)
		{
			*xd = zint_mod_small_signed(xs, dlen, p, p0i, R2, Rx);
			*yd = zint_mod_small_signed(ys, dlen, p, p0i, R2, Rx);
		}
	}

	/*
	 * We do not need Fd and Gd after that point.
	 */

	/*
	 * Compute our F and G modulo sufficiently many small primes.
	 */
	for (u = 0; u < llen; u ++) {
		uint32_t p, p0i, R2;
		uint32_t *gm, *igm, *fx, *gx, *Fp, *Gp;
		size_t v;

		/*
		 * All computations are done modulo p.
		 */
		p = primes[u].p;
		p0i = modp_ninv31(p);
		R2 = modp_R2(p, p0i);

		/*
		 * If we processed slen words, then f and g have been
		 * de-NTTized, and are in RNS; we can rebuild them.
		 */
		if (u == slen) {
			zint_rebuild_CRT(ft, slen, slen, n, primes, 1, t1);
			zint_rebuild_CRT(gt, slen, slen, n, primes, 1, t1);
		}

		gm = t1;
		igm = gm + n;
		fx = igm + n;
		gx = fx + n;

		modp_mkgm2(gm, igm, logn, primes[u].g, p, p0i);

		if (u < slen) {
			for (v = 0, x = ft + u, y = gt + u;
				v < n; v ++, x += slen, y += slen)
			{
				fx[v] = *x;
				gx[v] = *y;
			}
			modp_iNTT2_ext(ft + u, slen, igm, logn, p, p0i);
			modp_iNTT2_ext(gt + u, slen, igm, logn, p, p0i);
		} else {
			uint32_t Rx;

			Rx = modp_Rx((unsigned)slen, p, p0i, R2);
			for (v = 0, x = ft, y = gt;
				v < n; v ++, x += slen, y += slen)
			{
				fx[v] = zint_mod_small_signed(x, slen,
					p, p0i, R2, Rx);
				gx[v] = zint_mod_small_signed(y, slen,
					p, p0i, R2, Rx);
			}
			modp_NTT2(fx, gm, logn, p, p0i);
			modp_NTT2(gx, gm, logn, p, p0i);
		}

		/*
		 * Get F' and G' modulo p and in NTT representation
		 * (they have degree n/2). These values were computed in
		 * a previous step, and stored in Ft and Gt.
		 */
		Fp = gx + n;
		Gp = Fp + hn;
		for (v = 0, x = Ft + u, y = Gt + u;
			v < hn; v ++, x += llen, y += llen)
		{
			Fp[v] = *x;
			Gp[v] = *y;
		}
		modp_NTT2(Fp, gm, logn - 1, p, p0i);
		modp_NTT2(Gp, gm, logn - 1, p, p0i);

		/*
		 * Compute our F and G modulo p.
		 *
		 * General case:
		 *
		 *   we divide degree by d = 2 or 3
		 *   f'(x^d) = N(f)(x^d) = f * adj(f)
		 *   g'(x^d) = N(g)(x^d) = g * adj(g)
		 *   f'*G' - g'*F' = q
		 *   F = F'(x^d) * adj(g)
		 *   G = G'(x^d) * adj(f)
		 *
		 * We compute things in the NTT. We group roots of phi
		 * such that all roots x in a group share the same x^d.
		 * If the roots in a group are x_1, x_2... x_d, then:
		 *
		 *   N(f)(x_1^d) = f(x_1)*f(x_2)*...*f(x_d)
		 *
		 * Thus, we have:
		 *
		 *   G(x_1) = f(x_2)*f(x_3)*...*f(x_d)*G'(x_1^d)
		 *   G(x_2) = f(x_1)*f(x_3)*...*f(x_d)*G'(x_1^d)
		 *   ...
		 *   G(x_d) = f(x_1)*f(x_2)*...*f(x_{d-1})*G'(x_1^d)
		 *
		 * In all cases, we can thus compute F and G in NTT
		 * representation by a few simple multiplications.
		 * Moreover, in our chosen NTT representation, roots
		 * from the same group are consecutive in RAM.
		 */
		for (v = 0, x = Ft + u, y = Gt + u; v < hn;
			v ++, x += (llen << 1), y += (llen << 1))
		{
			uint32_t ftA, ftB, gtA, gtB;
			uint32_t mFp, mGp;

			ftA = fx[(v << 1) + 0];
			ftB = fx[(v << 1) + 1];
			gtA = gx[(v << 1) + 0];
			gtB = gx[(v << 1) + 1];
			mFp = modp_montymul(Fp[v], R2, p, p0i);
			mGp = modp_montymul(Gp[v], R2, p, p0i);
			x[0] = modp_montymul(gtB, mFp, p, p0i);
			x[llen] = modp_montymul(gtA, mFp, p, p0i);
			y[0] = modp_montymul(ftB, mGp, p, p0i);
			y[llen] = modp_montymul(ftA, mGp, p, p0i);
		}
		modp_iNTT2_ext(Ft + u, llen, igm, logn, p, p0i);
		modp_iNTT2_ext(Gt + u, llen, igm, logn, p, p0i);
	}

	/*
	 * Rebuild F and G with the CRT.
	 */
	zint_rebuild_CRT(Ft, llen, llen, n, primes, 1, t1);
	zint_rebuild_CRT(Gt, llen, llen, n, primes, 1, t1);

	/*
	 * Record the sizes of Ft, Gt so we know the amount of bits we need to use
	 * during the calculation of Ft and Gt.
	 */
	*unreduced_bits = 0;
	uint32_t *ptr = Ft;
	for (u = 0; u < 2*n; u++) {
		size_t sz = number_of_bits(ptr, llen);
		ptr += llen;
		if (sz > *unreduced_bits)
			*unreduced_bits = sz;
	}

	/*
	 * At that point, Ft, Gt, ft and gt are consecutive in RAM (in that
	 * order).
	 */

	/*
	 * Apply Babai reduction to bring back F and G to size slen.
	 *
	 * We use the FFT to compute successive approximations of the
	 * reduction coefficient. We first isolate the top bits of
	 * the coefficients of f and g, and convert them to floating
	 * point; with the FFT, we compute adj(f), adj(g), and
	 * 1/(f*adj(f)+g*adj(g)).
	 *
	 * Then, we repeatedly apply the following:
	 *
	 *   - Get the top bits of the coefficients of F and G into
	 *     floating point, and use the FFT to compute:
	 *        (F*adj(f)+G*adj(g))/(f*adj(f)+g*adj(g))
	 *
	 *   - Convert back that value into normal representation, and
	 *     round it to the nearest integers, yielding a polynomial k.
	 *     Proper scaling is applied to f, g, F and G so that the
	 *     coefficients fit on 32 bits (signed).
	 *
	 *   - Subtract k*f from F and k*g from G.
	 *
	 * Under normal conditions, this process reduces the size of F
	 * and G by some bits at each iteration. For constant-time
	 * operation, we do not want to measure the actual length of
	 * F and G; instead, we do the following:
	 *
	 *   - f and g are converted to floating-point, with some scaling
	 *     if necessary to keep values in the representable range.
	 *
	 *   - For each iteration, we _assume_ a maximum size for F and G,
	 *     and use the values at that size. If we overreach, then
	 *     we get zeros, which is harmless: the resulting coefficients
	 *     of k will be 0 and the value won't be reduced.
	 *
	 *   - We conservatively assume that F and G will be reduced by
	 *     at least 25 bits at each iteration.
	 *
	 * Even when reaching the bottom of the reduction, reduction
	 * coefficient will remain low. If it goes out-of-range, then
	 * something wrong occurred and the whole NTRU solving fails.
	 */

	/*
	 * Memory layout:
	 *  - We need to compute and keep adj(f), adj(g), and
	 *    1/(f*adj(f)+g*adj(g)) (sizes N, N and N/2 fp numbers,
	 *    respectively).
	 *  - At each iteration we need two extra fp buffer (N fp values),
	 *    and produce a k (N 32-bit words). k will be shared with one
	 *    of the fp buffers.
	 *  - To compute k*f and k*g efficiently (with the NTT), we need
	 *    some extra room; we reuse the space of the temporary buffers.
	 *
	 * Arrays of 'fpr' are obtained from the temporary array itself.
	 * We ensure that the base is at a properly aligned offset (the
	 * source array tmp[] is supposed to be already aligned).
	 */

	rt3 = align_fpr(tmp, t1);
	rt4 = rt3 + n;
	rt5 = rt4 + n;
	rt1 = rt5 + (n >> 1);
	k = (int32_t *)align_u32(tmp, rt1);
	rt2 = align_fpr(tmp, k + n);
	if (rt2 < (rt1 + n)) {
		rt2 = rt1 + n;
	}
	t1 = (uint32_t *)k + n;

	/*
	 * Get f and g into rt3 and rt4 as floating-point approximations.
	 *
	 * We need to "scale down" the floating-point representation of
	 * coefficients when they are too big. We want to keep the value
	 * below 2^310 or so. Thus, when values are larger than 10 words,
	 * we consider only the top 10 words. Array lengths have been
	 * computed so that average maximum length will fall in the
	 * middle or the upper half of these top 10 words.
	 */

	rlen = (slen > 10) ? 10 : slen;
	poly_big_to_fp(rt3, ft + slen - rlen, rlen, slen, logn);
	poly_big_to_fp(rt4, gt + slen - rlen, rlen, slen, logn);

	/*
	 * Values in rt3 and rt4 are downscaled by 2^(scale_fg).
	 */
	scale_fg = 31 * (int)(slen - rlen);

	/*
	 * Estimated boundaries for the maximum size (in bits) of the
	 * coefficients of (f,g). We use the measured average, and
	 * allow for a deviation of at most six times the standard
	 * deviation.
	 */
	minbl_fg = BITLENGTH[depth].avg - 6 * BITLENGTH[depth].std;
	maxbl_fg = BITLENGTH[depth].avg + 6 * BITLENGTH[depth].std;

	/*
	 * Compute 1/(f*adj(f)+g*adj(g)) in rt5. We also keep adj(f)
	 * and adj(g) in rt3 and rt4, respectively.
	 */

	Zf(FFT)(rt3, logn);
	Zf(FFT)(rt4, logn);
	Zf(poly_invnorm2_fft)(rt5, rt3, rt4, logn);
	Zf(poly_adj_fft)(rt3, logn);
	Zf(poly_adj_fft)(rt4, logn);

	/*
	 * Reduce F and G repeatedly.
	 *
	 * The expected maximum bit length of coefficients of F and G
	 * is kept in maxbl_FG, with the corresponding word length in
	 * FGlen.
	 */
	FGlen = llen;
	maxbl_FG = 31 * (int)llen;

	/*
	 * Each reduction operation computes the reduction polynomial
	 * "k". We need that polynomial to have coefficients that fit
	 * on 32-bit signed integers, with some scaling; thus, we use
	 * a descending sequence of scaling values, down to zero.
	 *
	 * The size of the coefficients of k is (roughly) the difference
	 * between the size of the coefficients of (F,G) and the size
	 * of the coefficients of (f,g). Thus, the maximum size of the
	 * coefficients of k is, at the start, maxbl_FG - minbl_fg;
	 * this is our starting scale value for k.
	 *
	 * We need to estimate the size of (F,G) during the execution of
	 * the algorithm; we are allowed some overestimation but not too
	 * much (poly_big_to_fp() uses a 310-bit window). Generally
	 * speaking, after applying a reduction with k scaled to
	 * scale_k, the size of (F,G) will be size(f,g) + scale_k + dd,
	 * where 'dd' is a few bits to account for the fact that the
	 * reduction is never perfect (intuitively, dd is on the order
	 * of sqrt(N), so at most 5 bits; we here allow for 10 extra
	 * bits).
	 *
	 * The size of (f,g) is not known exactly, but maxbl_fg is an
	 * upper bound.
	 */
	scale_k = maxbl_FG - minbl_fg;

	for (;;) {
		int scale_FG, dc, new_maxbl_FG;
		uint32_t scl, sch;
		fpr pdc, pt;

		/*
		 * Convert current F and G into floating-point. We apply
		 * scaling if the current length is more than 10 words.
		 */
		rlen = (FGlen > 10) ? 10 : FGlen;
		poly_big_to_fp(rt1, Ft + FGlen - rlen, rlen, llen, logn);
		poly_big_to_fp(rt2, Gt + FGlen - rlen, rlen, llen, logn);

		/*
		 * Values in rt1 and rt2 are downscaled by 2^(scale_FG).
		 */
		scale_FG = 31 * (int)(FGlen - rlen);

		/*
		 * Compute (F*adj(f)+G*adj(g))/(f*adj(f)+g*adj(g)) in rt2.
		 */
		Zf(FFT)(rt1, logn);
		Zf(FFT)(rt2, logn);
		Zf(poly_mul_fft)(rt1, rt3, logn);
		Zf(poly_mul_fft)(rt2, rt4, logn);
		Zf(poly_add)(rt2, rt1, logn);
		Zf(poly_mul_autoadj_fft)(rt2, rt5, logn);
		Zf(iFFT)(rt2, logn);

		/*
		 * (f,g) are scaled by 'scale_fg', meaning that the
		 * numbers in rt3/rt4 should be multiplied by 2^(scale_fg)
		 * to have their true mathematical value.
		 *
		 * (F,G) are similarly scaled by 'scale_FG'. Therefore,
		 * the value we computed in rt2 is scaled by
		 * 'scale_FG-scale_fg'.
		 *
		 * We want that value to be scaled by 'scale_k', hence we
		 * apply a corrective scaling. After scaling, the values
		 * should fit in -2^31-1..+2^31-1.
		 */
		dc = scale_k - scale_FG + scale_fg;

		/*
		 * We will need to multiply values by 2^(-dc). The value
		 * 'dc' is not secret, so we can compute 2^(-dc) with a
		 * non-constant-time process.
		 * (We could use ldexp(), but we prefer to avoid any
		 * dependency on libm. When using FP emulation, we could
		 * use our fpr_ldexp(), which is constant-time.)
		 */
		if (dc < 0) {
			dc = -dc;
			pt = fpr_two;
		} else {
			pt = fpr_onehalf;
		}
		pdc = fpr_one;
		while (dc != 0) {
			if ((dc & 1) != 0) {
				pdc = fpr_mul(pdc, pt);
			}
			dc >>= 1;
			pt = fpr_sqr(pt);
		}


		for (u = 0; u < n; u ++) {
			fpr xv;

			xv = fpr_mul(rt2[u], pdc);

			/*
			 * Sometimes the values can be out-of-bounds if
			 * the algorithm fails; we must not call
			 * fpr_rint() (and cast to int32_t) if the value
			 * is not in-bounds. Note that the test does not
			 * break constant-time discipline, since any
			 * failure here implies that we discard the current
			 * secret key (f,g).
			 */
			if (!fpr_lt(fpr_mtwo31m1, xv)
				|| !fpr_lt(xv, fpr_ptwo31m1))
			{
				return 0;
			}
			k[u] = (int32_t)fpr_rint(xv);
		}

		/*
		 * Values in k[] are integers. They really are scaled
		 * down by maxbl_FG - minbl_fg bits.
		 *
		 * If we are at low depth, then we use the NTT to
		 * compute k*f and k*g.
		 */
		sch = (uint32_t)(scale_k / 31);
		scl = (uint32_t)(scale_k % 31);
		if (depth <= DEPTH_INT_FG) {
			poly_sub_scaled_ntt(Ft, FGlen, llen, ft, slen, slen, k, sch, scl, logn, t1);
			poly_sub_scaled_ntt(Gt, FGlen, llen, gt, slen, slen, k, sch, scl, logn, t1);
		} else {
			poly_sub_scaled(Ft, FGlen, llen, ft, slen, slen, k, sch, scl, logn);
			poly_sub_scaled(Gt, FGlen, llen, gt, slen, slen, k, sch, scl, logn);
		}

		/*
		 * We compute the new maximum size of (F,G), assuming that
		 * (f,g) has _maximal_ length (i.e. that reduction is
		 * "late" instead of "early". We also adjust FGlen
		 * accordingly.
		 */
		new_maxbl_FG = scale_k + maxbl_fg + 10;
		if (new_maxbl_FG < maxbl_FG) {
			maxbl_FG = new_maxbl_FG;
			if ((int)FGlen * 31 >= maxbl_FG + 31) {
				FGlen --;
			}
		}

		/*
		 * We suppose that scaling down achieves a reduction by
		 * at least 25 bits per iteration. We stop when we have
		 * done the loop with an unscaled k.
		 */
		if (scale_k <= 0) {
			break;
		}
		scale_k -= 25;
		if (scale_k < 0) {
			scale_k = 0;
		}
	}

	/*
	 * If (F,G) length was lowered below 'slen', then we must take
	 * care to re-extend the sign.
	 */
	if (FGlen < slen) {
		for (u = 0; u < n; u ++, Ft += llen, Gt += llen) {
			size_t v;
			uint32_t sw;

			sw = -(Ft[FGlen - 1] >> 30) >> 1;
			for (v = FGlen; v < slen; v ++) {
				Ft[v] = sw;
			}
			sw = -(Gt[FGlen - 1] >> 30) >> 1;
			for (v = FGlen; v < slen; v ++) {
				Gt[v] = sw;
			}
		}
	}

	/*
	 * Compress encoding of all values to 'slen' words (this is the
	 * expected output format).
	 */
	for (u = 0, x = tmp, y = tmp;
		u < (n << 1); u ++, x += slen, y += llen)
	{
		memmove(x, y, slen * sizeof *y);
	}
	return 1;
}

// =============================================================================
// | END                                                                       |
// =============================================================================

const size_t logn = 9, n = 1<<logn;

void sample_FG(inner_shake256_context *rng, int8_t *f, int8_t *g,
	size_t *bits_FG, uint8_t *tmp)
{
	/*
	 * Algorithm is the following:
	 *
	 *  - Generate f and g with the Gaussian distribution.
	 *
	 *  - If either Res(f,phi) or Res(g,phi) is even, try again.
	 *
	 *  - Solve the NTRU equation fG - gF = 1; if the solving fails,
	 *    try again. Usual failure condition is when Res(f,phi)
	 *    and Res(g,phi) are not prime to each other.
	 */
	prng p;
	/*
	 * In the binary case, coefficients of f and g are generated
	 * independently of each other, with a discrete Gaussian
	 * distribution of standard deviation 1.425. Then,
	 * the two vectors have expected norm 2n * 1.425.
	 *
	 * We require that Res(f,phi) and Res(g,phi) are both odd (the
	 * NTRU equation solver requires it).
	 */
	while (1) {
		// Normal sampling. We use a fast PRNG seeded from our SHAKE context ('rng').
sample:
		Zf(prng_init)(&p, rng);

		poly_small_mkgauss(&p, f, logn);
		poly_small_mkgauss(&p, g, logn);
		
		if (!solve_NTRU_deepest(logn, f, g, (uint32_t *)tmp))
			continue;

		unsigned depth = logn;

		size_t len = MAX_BL_SMALL[logn];
		uint32_t *ptr = (uint32_t *)tmp;
		size_t szF = number_of_bits(ptr, len);
		size_t szG = number_of_bits(ptr + len, len);
		bits_FG[logn] = szF > szG ? szF : szG;

		while (depth -- > 0) {
			// use the modified version:
			if (!_solve_NTRU_intermediate(logn, f, g,
					&bits_FG[depth], depth, (uint32_t *)tmp)) {
				goto sample;
			}
		}

		break;
	}
}

void sample_FG_sizes(inner_shake256_context *rng, uint8_t *tmp)
{
	int8_t f[n], g[n];
	long long sum_b[logn + 1], sum_bsq[logn + 1];
	size_t bits_FG[logn + 1];

	memset(sum_b, 0, sizeof sum_b);
	memset(sum_bsq, 0, sizeof sum_bsq);

	size_t nsamples = 100;
	for (size_t i = 0; i < nsamples; i++) {
		if (i % (nsamples/100) == 0) printf("."), fflush(stdout);
		sample_FG(rng, f, g, bits_FG, tmp);
		for (size_t depth = 0; depth <= logn; depth++) {
			long long x = bits_FG[depth];
			sum_b[depth] += x;
			sum_bsq[depth] += x * x;
		}
	}

	printf("\n");
	// calculate the average and standard deviation
	for (size_t depth = 0; depth <= logn; depth++) {
		double avg = ((double)sum_b[depth]) / nsamples;
		double stddev = sqrt(((double)sum_bsq[depth]) / nsamples - avg*avg);
		size_t nr_ints = (int)(avg + 6.0 * stddev + 30) / 31;
		printf("Depth %zu: %.2f (%.2f) --> (%zu ints)\n", depth, avg, stddev, nr_ints);
	}

	printf("static const size_t MAX_BL_LARGE[9] = {\n\t");
	for (size_t depth = 0; depth < logn; depth++) {
		double avg = ((double)sum_b[depth]) / nsamples;
		double stddev = sqrt(((double)sum_bsq[depth]) / nsamples - avg*avg);
		size_t nr_ints = (int)(avg + 6.0 * stddev + 30) / 31;
		printf("%zu", nr_ints);
		if (depth == logn-1) printf("\n");
		else printf(", ");
	};
	printf("};\n");

}


// =============================================================================
uint8_t tmp[68 * 512];
int main() {
	srand(time(NULL));

	inner_shake256_context sc;
	uint8_t seed[48];
	for (int i = 0; i < 48; i++)
		seed[i] = (uint8_t) rand();
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	sample_FG_sizes(&sc, tmp);
	return 0;
}
