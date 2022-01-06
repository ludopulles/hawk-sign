#include <stdint.h>

#include "inner.h"

// ================================================================================

static const size_t MAX_BL_SMALL64[] = {
//  1, 1, 2, 2, 4, 7, 14, 27, 53, 106, 209 // (FALCON)
	1, 1, 2, 2, 4, 7, 13, 22, 43, 90, 209
};

static const size_t MAX_BL_LARGE64[] = {
//  2, 2, 5, 7, 12, 21, 40, 78, 157, 308 // (FALCON)
	2, 2, 5, 7, 12, 19, 30, 55, 120, 308
};

/*
 * Average and standard deviation for the maximum size (in bits) of
 * coefficients of (f,g), depending on depth. These values are used
 * to compute bounds for Babai's reduction.
 */
static const struct {
	int avg;
	int std;
} BITLENGTH64[] = {
	{    4,  0 },
	{   11,  1 },
	{   24,  1 },
	{   50,  1 },
	{  102,  1 },
	{  202,  2 },
	{  300,  4 },
	{  524,  5 }, // {  794,  5 },
	{ 1277,  8 }, // { 1577,  8 },
	{ 3138, 13 },
	{ 6308, 25 }
};

// ================================================================================

/*
 * Input: f,g of degree N = 2^logn; 'depth' is used only to get their
 * individual length.
 *
 * Output: f',g' of degree N/2, with the length for 'depth+1'.
 *
 * Values are in RNS; input and/or output may also be in NTT.
 */
static void
make_fg_step64(uint32_t *data, unsigned logn, unsigned depth,
	int in_ntt, int out_ntt)
{
	size_t n, hn, u;
	size_t slen, tlen;
	uint32_t *fd, *gd, *fs, *gs, *gm, *igm, *t1;
	const small_prime *primes;

	n = (size_t)1 << logn;
	hn = n >> 1;
	slen = MAX_BL_SMALL64[depth];
	tlen = MAX_BL_SMALL64[depth + 1];
	primes = PRIMES;

	/*
	 * Prepare room for the result.
	 */
	fd = data;
	gd = fd + hn * tlen;
	fs = gd + hn * tlen;
	gs = fs + n * slen;
	gm = gs + n * slen;
	igm = gm + n;
	t1 = igm + n;
	memmove(fs, data, 2 * n * slen * sizeof *data);

	/*
	 * First slen words: we use the input values directly, and apply
	 * inverse NTT as we go.
	 */
	for (u = 0; u < slen; u ++) {
		uint32_t p, p0i, R2;
		size_t v;
		uint32_t *x;

		p = primes[u].p;
		p0i = modp_ninv31(p);
		R2 = modp_R2(p, p0i);
		modp_mkgm2(gm, igm, logn, primes[u].g, p, p0i);

		for (v = 0, x = fs + u; v < n; v ++, x += slen) {
			t1[v] = *x;
		}
		if (!in_ntt) {
			modp_NTT2(t1, gm, logn, p, p0i);
		}
		for (v = 0, x = fd + u; v < hn; v ++, x += tlen) {
			uint32_t w0, w1;

			w0 = t1[(v << 1) + 0];
			w1 = t1[(v << 1) + 1];
			*x = modp_montymul(
				modp_montymul(w0, w1, p, p0i), R2, p, p0i);
		}
		if (in_ntt) {
			modp_iNTT2_ext(fs + u, slen, igm, logn, p, p0i);
		}

		for (v = 0, x = gs + u; v < n; v ++, x += slen) {
			t1[v] = *x;
		}
		if (!in_ntt) {
			modp_NTT2(t1, gm, logn, p, p0i);
		}
		for (v = 0, x = gd + u; v < hn; v ++, x += tlen) {
			uint32_t w0, w1;

			w0 = t1[(v << 1) + 0];
			w1 = t1[(v << 1) + 1];
			*x = modp_montymul(
				modp_montymul(w0, w1, p, p0i), R2, p, p0i);
		}
		if (in_ntt) {
			modp_iNTT2_ext(gs + u, slen, igm, logn, p, p0i);
		}

		if (!out_ntt) {
			modp_iNTT2_ext(fd + u, tlen, igm, logn - 1, p, p0i);
			modp_iNTT2_ext(gd + u, tlen, igm, logn - 1, p, p0i);
		}
	}

	/*
	 * Since the fs and gs words have been de-NTTized, we can use the
	 * CRT to rebuild the values.
	 */
	zint_rebuild_CRT(fs, slen, slen, n, primes, 1, gm);
	zint_rebuild_CRT(gs, slen, slen, n, primes, 1, gm);

	/*
	 * Remaining words: use modular reductions to extract the values.
	 */
	for (u = slen; u < tlen; u ++) {
		uint32_t p, p0i, R2, Rx;
		size_t v;
		uint32_t *x;

		p = primes[u].p;
		p0i = modp_ninv31(p);
		R2 = modp_R2(p, p0i);
		Rx = modp_Rx((unsigned)slen, p, p0i, R2);
		modp_mkgm2(gm, igm, logn, primes[u].g, p, p0i);
		for (v = 0, x = fs; v < n; v ++, x += slen) {
			t1[v] = zint_mod_small_signed(x, slen, p, p0i, R2, Rx);
		}
		modp_NTT2(t1, gm, logn, p, p0i);
		for (v = 0, x = fd + u; v < hn; v ++, x += tlen) {
			uint32_t w0, w1;

			w0 = t1[(v << 1) + 0];
			w1 = t1[(v << 1) + 1];
			*x = modp_montymul(
				modp_montymul(w0, w1, p, p0i), R2, p, p0i);
		}
		for (v = 0, x = gs; v < n; v ++, x += slen) {
			t1[v] = zint_mod_small_signed(x, slen, p, p0i, R2, Rx);
		}
		modp_NTT2(t1, gm, logn, p, p0i);
		for (v = 0, x = gd + u; v < hn; v ++, x += tlen) {
			uint32_t w0, w1;

			w0 = t1[(v << 1) + 0];
			w1 = t1[(v << 1) + 1];
			*x = modp_montymul(
				modp_montymul(w0, w1, p, p0i), R2, p, p0i);
		}

		if (!out_ntt) {
			modp_iNTT2_ext(fd + u, tlen, igm, logn - 1, p, p0i);
			modp_iNTT2_ext(gd + u, tlen, igm, logn - 1, p, p0i);
		}
	}
}


/*
 * Compute f and g at a specific depth, in RNS notation.
 *
 * Returned values are stored in the data[] array, at slen words per integer.
 *
 * Conditions:
 *   0 <= depth <= logn
 *
 * Space use in data[]: enough room for any two successive values (f', g',
 * f and g).
 */
static void
make_fg64(uint32_t *data, const int8_t *f, const int8_t *g,
	unsigned logn, unsigned depth, int out_ntt)
{
	size_t n, u;
	uint32_t *ft, *gt, p0;
	unsigned d;
	const small_prime *primes;

	n = MKN(logn);
	ft = data;
	gt = ft + n;
	primes = PRIMES;
	p0 = primes[0].p;
	for (u = 0; u < n; u ++) {
		ft[u] = modp_set(f[u], p0);
		gt[u] = modp_set(g[u], p0);
	}

	if (depth == 0 && out_ntt) {
		uint32_t *gm, *igm;
		uint32_t p0i;

		p0i = modp_ninv31(p0);
		gm = gt + n;
		igm = gm + n;
		modp_mkgm2(gm, igm, logn, primes[0].g, p0, p0i);
		modp_NTT2(ft, gm, logn, p0, p0i);
		modp_NTT2(gt, gm, logn, p0, p0i);
		return;
	}

	for (d = 0; d < depth; d ++) {
		make_fg_step64(data, logn - d, d,
			d != 0, (d + 1) < depth || out_ntt);
	}
}

/*
 * Solving the NTRU equation for q = 1, deepest level: compute the resultants
 * of f and g with X^N+1, and use binary GCD. The F and G values are returned
 * in tmp[].
 *
 * Returned value: 1 on success, 0 on error.
 */
static int
lilipu_solve_NTRU_deepest64(unsigned logn_top,
	const int8_t *f, const int8_t *g, uint32_t *tmp)
{
	size_t len;
	uint32_t *Fp, *Gp, *fp, *gp, *t1;
	const small_prime *primes;

	len = MAX_BL_SMALL64[logn_top];
	primes = PRIMES;

	Fp = tmp;
	Gp = Fp + len;
	fp = Gp + len;
	gp = fp + len;
	t1 = gp + len;

	make_fg64(fp, f, g, logn_top, logn_top, 0);

	/*
	 * We use the CRT to rebuild the resultants as big integers.
	 * There are two such big integers. The resultants are always
	 * nonnegative.
	 */
	zint_rebuild_CRT(fp, len, len, 2, primes, 0, t1);

/*
	printf("Initial normDown(f), normDown(g):\n");
	for (size_t u = 0; u < len; u++) printf("%u,", fp[u]);
	printf("\n");
	for (size_t u = 0; u < len; u++) printf("%u,", gp[u]);
	printf("\n");
*/

	/*
	 * Apply the binary GCD. The zint_bezout() function works only
	 * if both inputs are odd.
	 *
	 * We can test on the result and return 0 because that would
	 * imply failure of the NTRU solving equation, and the (f,g)
	 * values will be abandoned in that case.
	 */
	return zint_bezout(Gp, Fp, fp, gp, len, t1);
}

static int
solve_NTRU_intermediate64(unsigned logn_top,
	const int8_t *f, const int8_t *g, unsigned depth, uint32_t *tmp)
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
	slen = MAX_BL_SMALL64[depth];
	dlen = MAX_BL_SMALL64[depth + 1];
	llen = MAX_BL_LARGE64[depth];
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
	make_fg64(ft, f, g, logn_top, depth, 1);

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

/*	printf("Initial ft, gt:\n");
	for (size_t v = 0; v < n; v++) {
		for (u = 0; u < slen; u++) printf("%u,", ft[v*slen + u]);
		printf(" |\n");
	}
	for (size_t v = 0; v < n; v++) {
		for (u = 0; u < slen; u++) printf("%u,", gt[v*slen + u]);
		printf(" |\n");
	}

	printf("Initial Ft, Gt:\n");
	for (size_t v = 0; v < n; v++) {
		for (u = 0; u < llen; u++) printf("%u,", Ft[v*llen + u]);
		printf(" |\n");
	}
	for (size_t v = 0; v < n; v++) {
		for (u = 0; u < llen; u++) printf("%u,", Gt[v*llen + u]);
		printf(" |\n");
	} */

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
	minbl_fg = BITLENGTH64[depth].avg - 6 * BITLENGTH64[depth].std;
	maxbl_fg = BITLENGTH64[depth].avg + 6 * BITLENGTH64[depth].std;

	/*
	 * Compute 1/(f*adj(f)+g*adj(g)) in rt5. We also keep adj(f)
	 * and adj(g) in rt3 and rt4, respectively.
	 */

/* 	printf("Scale_fg: %d = 31*(%d - %d)\n", scale_fg, (int)slen, (int)rlen);
	for (u = 0; u < n; u++) printf("%e ", rt3[u].v);
	printf("\n");
	for (u = 0; u < n; u++) printf("%e ", rt4[u].v);
	printf("\n"); */

	Zf(FFT)(rt3, logn);
	Zf(FFT)(rt4, logn);
	Zf(poly_invnorm2_fft)(rt5, rt3, rt4, logn);
	Zf(poly_adj_fft)(rt3, logn);
	Zf(poly_adj_fft)(rt4, logn);

/*	for (u = 0; u < n; u++) printf("%e ", rt3[u].v);
	printf("\n");
	for (u = 0; u < n; u++) printf("%e ", rt4[u].v);
	printf("\n");
	for (u = 0; u < n; u++) printf("%e ", rt5[u].v);
	printf("\n"); */

	// for (u = 0; u < n; u++) assert(rt5[u].v == rt5[u].v);

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
		scale_FG = 31 * (int)(FGlen - rlen);
		poly_big_to_fp(rt1, Ft + FGlen - rlen, rlen, llen, logn);
		poly_big_to_fp(rt2, Gt + FGlen - rlen, rlen, llen, logn);

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
/* printf("Scales: %d %d %d (=%d - %d)\n", scale_fg, scale_FG, scale_k, maxbl_FG, minbl_fg);
printf("Values found: \n");
for (u = 0; u < n; u ++) {
	printf("%e ", fpr_mul(rt2[u], pdc).v);
}
printf("\n"); */
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
			poly_sub_scaled_ntt(Ft, FGlen, llen, ft, slen, slen,
				k, sch, scl, logn, t1);
			poly_sub_scaled_ntt(Gt, FGlen, llen, gt, slen, slen,
				k, sch, scl, logn, t1);
		} else {
			poly_sub_scaled(Ft, FGlen, llen, ft, slen, slen,
				k, sch, scl, logn);
			poly_sub_scaled(Gt, FGlen, llen, gt, slen, slen,
				k, sch, scl, logn);
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
