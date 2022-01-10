#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#include "codec.c"
#include "common.c"
#include "fft.c"
#include "fpr.c"


int
Zf(compute_public)(uint16_t *h,
	const int8_t *f, const int8_t *g, unsigned logn, uint8_t *tmp)
{ return 1; } // needed for compilation....

#include "keygen.c"
// #include "nist.c"
#include "rng.c"
#include "shake.c"
#include "sign.c"

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

const size_t UPPER_BOUND_SIZE[11] = {
	// 1, 1, 4, 4, 16, 16, 32, 64, 128, 128, 256
	1, 1, 2, 2, 4, 6, 12, 22, 42, 83, -1
};

/*
 * Input: f,g of degree N = 2^logn; 'depth' is used only to get their
 * individual length.
 *
 * Output: f',g' of degree N/2, with the length for 'depth+1'.
 *
 * Values are in RNS; input and/or output may also be in NTT.
 */
static void
lilipu_make_fg_step(uint32_t *data, unsigned logn, unsigned depth,
	int in_ntt, int out_ntt)
{
	size_t n, hn, u;
	size_t slen, tlen;
	uint32_t *fd, *gd, *fs, *gs, *gm, *igm, *t1;
	const small_prime *primes;

	n = (size_t)1 << logn;
	hn = n >> 1;
	slen = UPPER_BOUND_SIZE[depth];
	tlen = UPPER_BOUND_SIZE[depth + 1];
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
lilipu_make_fg(uint32_t *data, const int8_t *f, const int8_t *g,
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
		uint32_t p, p0i;

		p = primes[0].p;
		p0i = modp_ninv31(p);
		gm = gt + n;
		igm = gm + MKN(logn);
		modp_mkgm2(gm, igm, logn, primes[0].g, p, p0i);
		modp_NTT2(ft, gm, logn, p, p0i);
		modp_NTT2(gt, gm, logn, p, p0i);
		return;
	}

	for (d = 0; d < depth; d ++) {
		lilipu_make_fg_step(data, logn - d, d,
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
lilipu_solve_NTRU_deepest(unsigned logn_top,
	const int8_t *f, const int8_t *g, uint32_t *tmp)
{
	size_t len;
	uint32_t *Fp, *Gp, *fp, *gp, *t1;
	const small_prime *primes;

	len = UPPER_BOUND_SIZE[logn_top]; // should be large enough
	primes = PRIMES;

	Fp = tmp;
	Gp = Fp + len;
	fp = Gp + len;
	gp = fp + len;
	t1 = gp + len;

	lilipu_make_fg(fp, f, g, logn_top, logn_top, 0);

	/*
	 * We use the CRT to rebuild the resultants as big integers.
	 * There are two such big integers. The resultants are always
	 * nonnegative.
	 */
	zint_rebuild_CRT(fp, len, len, 2, primes, 0, t1);

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


// const fpr sigma_kg = fpr_div(fpr_of(1425), fpr_of(1000));
const fpr isigma_kg = { v: 1.0 / 1.425 }; // fpr_div(fpr_of(1000), fpr_of(1425));

/*
 * Generate a random polynomial with a Gaussian distribution. This function
 * also makes sure that the resultant of the polynomial with phi is odd.
 */
static void
lilipu_poly_small_mkgauss(samplerZ samp, void *samp_ctx, int8_t *f, unsigned logn, int lim)
{
	size_t n, u;
	int s;
	unsigned mod2;

	n = MKN(logn);
	mod2 = 0;

	for (u = n; u -- > 1; ) {
		do {
			s = samp(samp_ctx, fpr_zero, isigma_kg);
			/*
			 * We need the coefficient to fit within -127..+127;
			 * realistically, this is always the case except for
			 * the very low degrees (N = 2 or 4), for which there
			 * is no real security anyway.
			 */
		} while (s <= -lim || s >= lim);
		mod2 ^= (unsigned)(s & 1);
		f[u] = (int8_t)s;
	}

	do {
		s = samp(samp_ctx, fpr_zero, isigma_kg);
		/*
		 * We need the sum of all coefficients to be 1; otherwise,
		 * the resultant of the polynomial with X^N+1 will be even,
		 * and the binary GCD will fail.
		 */
	} while (s <= -lim || s >= lim || mod2 == (unsigned)(s & 1));
	f[0] = (int8_t)s;
}

size_t number_of_bits(uint32_t *fp, size_t sz)
{
	size_t sgn = (fp[sz-1] >> 30) & 1; // 0 = positive, 1 = negative
	while (sz > 0 && (fp[sz-1] == 0 || fp[sz-1] == 2147483647U)) sz--;

	if (sz == 0) return 1 + sgn; // 0 or -1 as value

	size_t res = (sz-1) * 31;
	for (size_t b = 30; b >= 0; b--)
		if (((fp[sz-1] >> b) & 1) != sgn) {
			if (sz == 1 && fp[0] == 2147483646U) assert(b == 0);
			if (sz == 1 && fp[0] == 3) assert(b == 1);
			return res + (b+1);
		}
	fprintf(stderr, "Unexpected value: %zu\n", fp[sz-1]);
	assert(0);
}

const size_t logn = 9, n = 1<<logn;

void sample_fg(inner_shake256_context *rng, int8_t *f, int8_t *g, uint8_t *tmp)
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
	int lim = 1 << (Zf(max_fg_bits)[logn] - 1);
	/*
	 * In the binary case, coefficients of f and g are generated
	 * independently of each other, with a discrete Gaussian
	 * distribution of standard deviation 1/isigma_kg. Then,
	 * the two vectors have expected norm 2n/isigma_kg.
	 *
	 * We require that Res(f,phi) and Res(g,phi) are both odd (the
	 * NTRU equation solver requires it).
	 */
	do {
		// Normal sampling. We use a fast PRNG seeded from our SHAKE context ('rng').
		sampler_context spc;
		samplerZ samp;
		void *samp_ctx;
		spc.sigma_min = fpr_sigma_min[logn];
		falcon_inner_prng_init(&spc.p, rng);
		samp = Zf(sampler);
		samp_ctx = &spc;

		lilipu_poly_small_mkgauss(samp, samp_ctx, f, logn, lim);
		lilipu_poly_small_mkgauss(samp, samp_ctx, g, logn, lim);
	} while (!lilipu_solve_NTRU_deepest(logn, f, g, (uint32_t *)tmp));
}

void sample_fg_sizes(inner_shake256_context *rng, uint8_t *tmp)
{
	int8_t f[n], g[n];

	long long sum_b[logn + 1];
	long long sum_bsq[logn + 1];

	memset(sum_b, 0, sizeof sum_b);
	memset(sum_bsq, 0, sizeof sum_bsq);

	size_t nsamples = 100000, len;
	uint32_t *fp, *gp, *t1;
	for (size_t i = 0; i < nsamples; i++) {
		sample_fg(rng, f, g, tmp);
		// now do some statistics
		for (int depth = 0; depth <= logn; depth++) {
			len = UPPER_BOUND_SIZE[depth]; // should be large enough

			fp = (uint32_t *)tmp;
			gp = fp + (len << (logn - depth));
			t1 = gp + (len << (logn - depth));

			memset(tmp, 0, (len << (logn - depth)) * sizeof fp);
			lilipu_make_fg(fp, f, g, logn, depth, 0);
			// Rebuild fp, gp as polynomials of big integers
			zint_rebuild_CRT(fp, len, len, 2 << (logn - depth), PRIMES, 1, t1);

			// determine sizes of fp, gp
			uint32_t *ptr = fp;
			size_t longest = 0;
			for (size_t u = 0; u < (2 << (logn - depth)); u++, ptr += len) {
				// for (size_t v = 0; v < len; v++) printf("%zu ", ptr[v]);
				// printf("\n");
				size_t sz = number_of_bits(ptr, len);
				if (sz > longest) longest = sz;
			}
			assert(longest >= 2);
			// printf("\nDepth %d: %d\n", depth, (int) longest);
			sum_b[depth] += longest;
			sum_bsq[depth] += (long long)longest * (long long)longest;
		}
	}

	// calculate the average and standard deviation
	for (int depth = 0; depth <= logn; depth++) {
		double avg = ((double)sum_b[depth]) / nsamples;
		double stddev = sqrt(((double)sum_bsq[depth]) / nsamples - avg*avg);
		size_t nr_ints = (int)(avg + 6.0 * stddev + 30) / 31;
		printf("Depth %d: %.2f (%.2f) --> (%zu ints)\n", depth, avg, stddev, nr_ints);
	}
}


// =============================================================================
uint8_t tmp[16 * 1024 * 1024];
int main() {
	srand(time(NULL));

	inner_shake256_context sc;
	uint8_t seed[48];
	for (int i = 0; i < 48; i++)
		seed[i] = (uint8_t) rand();
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	sample_fg_sizes(&sc, tmp);
	return 0;
}
