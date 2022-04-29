#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>


#define __CONFIG_fg
// since MAX_BL_SMALL is const, we use a #define __CONFIG_fg to mark in 
// keygen.c that we need to use large upper bounds instead.
#include "keygen.c"

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

/* unused:
void print_big(uint32_t *x, size_t xlen)
{
	uint32_t *y = malloc(xlen * sizeof(uint32_t));
	memcpy(y, x, xlen * sizeof(uint32_t));

	uint8_t *out = malloc(xlen * 31);
	size_t nout = 0, ylen = xlen, sgn = y[ylen-1]>>30;
	if (sgn) {
		for (size_t u = 0; u < ylen; u++)
			y[u] = (1U<<31) - 1 - y[u];
		y[0]++;
		printf("-");
	}

	while (ylen > 0) {
		uint64_t C = 0;
		for (size_t p = ylen; p -- > 0;  C %= 10) {
			C = (C << 31) | y[p];
			y[p] = C / 10;
		}
		// now C is the remainder
		out[nout++] = (int8_t) C;
		while (ylen > 0 && y[ylen - 1] == 0) --ylen;
	}
	while (nout -- > 0)
		printf("%d", out[nout]);
	free(y);
	free(out);
}
*/

const size_t logn = 9, n = MKN(logn);

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
	 *
	 * In the binary case, coefficients of f and g are generated
	 * independently of each other, with a discrete Gaussian
	 * distribution of standard deviation 1.425 Then,
	 * the two vectors have expected norm 2n * 1.425.
	 *
	 * We require that Res(f,phi) and Res(g,phi) are both odd (the
	 * NTRU equation solver requires it).
	 */
	do {
		// Normal sampling. We use a fast PRNG seeded from our SHAKE context ('rng').
		prng p;
		Zf(prng_init)(&p, rng);

		poly_small_mkgauss(&p, f, logn);
		poly_small_mkgauss(&p, g, logn);
	} while (!solve_NTRU_deepest(logn, f, g, (uint32_t *)tmp));
}

void sample_fg_sizes(inner_shake256_context *rng, uint8_t *tmp)
{
	int8_t f[n], g[n];
	long long sum_b[logn + 1], sum_bsq[logn + 1];

	memset(sum_b, 0, sizeof sum_b);
	memset(sum_bsq, 0, sizeof sum_bsq);

	size_t nsamples = 100, len;
	uint32_t *fp, *gp, *t1;
	for (size_t i = 0; i < nsamples; i++) {
		sample_fg(rng, f, g, tmp);
		// now do some statistics
		for (size_t depth = 0; depth <= logn; depth++) {
			len = MAX_BL_SMALL[depth]; // should be large enough

			fp = (uint32_t *)tmp;
			gp = fp + (len << (logn - depth));
			t1 = gp + (len << (logn - depth));

			memset(tmp, 0, (len << (logn - depth)) * sizeof fp);
			make_fg(fp, f, g, logn, depth, 0);
			// Rebuild fp, gp as polynomials of big integers
			zint_rebuild_CRT(fp, len, len, 2 << (logn - depth), PRIMES, 1, t1);

			// determine sizes of fp, gp
			uint32_t *ptr = fp;
			size_t longest = 0;
			for (size_t u = 0; u < (2U << (logn - depth)); u++, ptr += len) {
				size_t sz = number_of_bits(ptr, len);
				if (sz > longest) longest = sz;
			}
			assert(longest >= 2);
			sum_b[depth] += longest;
			sum_bsq[depth] += (long long)longest * (long long)longest;
		}
	}

	// calculate the average and standard deviation
	for (size_t depth = 0; depth <= logn; depth++) {
		double avg = ((double)sum_b[depth]) / nsamples;
		double stddev = sqrt(((double)sum_bsq[depth]) / nsamples - avg*avg);
		size_t nr_ints = (int)(avg + 6.0 * stddev + 30) / 31;
		printf("Depth %zu: %.2f (%.2f) --> (%zu ints)\n", depth, avg, stddev, nr_ints);
	}

	printf("static const size_t MAX_BL_SMALL[10] = {\n\t");
	for (size_t depth = 0; depth <= logn; depth++) {
		double avg = ((double)sum_b[depth]) / nsamples;
		double stddev = sqrt(((double)sum_bsq[depth]) / nsamples - avg*avg);
		int nr_ints = (int)(avg + 6.0 * stddev + 30) / 31;
		printf("%d", nr_ints);
		if (depth == logn) printf("\n");
		else printf(", ");
	};
	printf("};\n");

	printf("static const struct {\n\tint avg, std;\n} BITLENGTH[10] = {\n");
	for (size_t depth = 0; depth <= logn; depth++) {
		double avg = ((double)sum_b[depth]) / nsamples;
		double stddev = sqrt(((double)sum_bsq[depth]) / nsamples - avg*avg);
		printf("\t{ %d, %d }", (int) llroundl(avg), (int) ceil(stddev));
		if (depth == logn) printf("\n");
		else printf(",\n");
	};
	printf("};\n");
}


// =============================================================================
uint8_t tmp[26 * 512];

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
