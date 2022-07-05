/*
 * Speed benchmark code for Hawk implementation.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*
 * This code uses only the external API.
 */

#include "../hawk.h"

#ifdef __AVX2__
/*
 * Use floating points when the AVX2 instruction set is available, since a CPU
 * supporting AVX2 probably has a FPU.
 */
#define HAWK_AVX
#endif

static void *
xmalloc(size_t len)
{
	void *buf;

	if (len == 0) {
		return NULL;
	}
	buf = malloc(len);
	if (buf == NULL) {
		fprintf(stderr, "memory allocation error\n");
		exit(EXIT_FAILURE);
	}
	return buf;
}

static void
xfree(void *buf)
{
	if (buf != NULL) {
		free(buf);
	}
}

/*
 * Benchmark function takes an opaque context and an iteration count;
 * it returns 0 on success, a negative error code on error.
 */
typedef int (*bench_fun)(void *ctx, unsigned long num);

/*
 * Returned value is the time per iteration in nanoseconds. If the
 * benchmark function reports an error, 0.0 is returned.
 */
static double
do_bench(bench_fun bf, void *ctx, double threshold)
{
	unsigned long num;
	int r;

	/*
	 * Alsways do a few blank runs to "train" the caches and branch
	 * prediction.
	 */
	r = bf(ctx, 5);
	if (r != 0) {
		fprintf(stderr, "ERR: %d\n", r);
		return 0.0;
	}

	num = 1;
	for (;;) {
		clock_t begin, end;
		double tt;

		begin = clock();
		r = bf(ctx, num);
		end = clock();
		if (r != 0) {
			fprintf(stderr, "ERR: %d\n", r);
			return 0.0;
		}
		tt = (double)(end - begin) / (double)CLOCKS_PER_SEC;
		if (tt >= threshold) {
			return tt * 1000000000.0 / (double)num;
		}

		/*
		 * If the function ran for less than 0.1 seconds then
		 * we simply double the iteration number; otherwise, we
		 * use the run time to try to get a "correct" number of
		 * iterations quickly.
		 */
		if (tt < 0.1) {
			num <<= 1;
		} else {
			unsigned long num2;

			num2 = (unsigned long)((double)num
				* (threshold * 1.1) / tt);
			if (num2 <= num) {
				num2 = num + 1;
			}
			num = num2;
		}
	}
}

typedef struct {
	unsigned logn;
	shake256_context rng;
	uint8_t *tmp, *pk, *sk, *esk, *uncompsig, *sig;
	size_t tmp_len, uncompsig_len, sig_len;
} bench_context;

#define CC(x)   do { \
		int ccr = (x); \
		if (ccr != 0) { \
			return ccr; \
		} \
	} while (0)

static int
bench_keygen(void *ctx, unsigned long num)
{
	bench_context *bc;

	bc = ctx;
	while (num -- > 0) {
		CC(hawk_keygen_make(&bc->rng, bc->logn,
			bc->sk, HAWK_SECKEY_SIZE(bc->logn),
			bc->pk, HAWK_PUBKEY_SIZE[bc->logn],
			bc->tmp, bc->tmp_len));
	}
	return 0;
}

static int
bench_expand_seckey(void *ctx, unsigned long num)
{
	bench_context *bc;

	bc = ctx;
	while (num -- > 0) {
		CC(hawk_expand_seckey(
			bc->esk, HAWK_EXPANDEDKEY_SIZE(bc->logn),
			bc->sk, HAWK_SECKEY_SIZE(bc->logn)));
	}
	return 0;
}

static int
bench_uncompressed_sign(void *ctx, unsigned long num)
{
	bench_context *bc;

	bc = ctx;
	while (num -- > 0) {
		bc->uncompsig_len = HAWK_UNCOMPRESSED_SIG_COMPACT_MAXSIZE(bc->logn);
		CC(hawk_uncompressed_sign(&bc->rng,
			bc->uncompsig, &bc->uncompsig_len, HAWK_SIG_COMPACT,
			bc->sk, HAWK_SECKEY_SIZE(bc->logn),
			"data", 4, bc->tmp, bc->tmp_len));
	}
	return 0;
}

static int
bench_sign_dyn(void *ctx, unsigned long num)
{
	bench_context *bc;

	bc = ctx;
	while (num -- > 0) {
		bc->sig_len = HAWK_SIG_COMPACT_MAXSIZE(bc->logn);
		CC(hawk_sign_dyn(&bc->rng,
			bc->sig, &bc->sig_len, HAWK_SIG_COMPACT,
			bc->sk, HAWK_SECKEY_SIZE(bc->logn),
			"data", 4, bc->tmp, bc->tmp_len));
	}
	return 0;
}

static int
bench_sign(void *ctx, unsigned long num)
{
	bench_context *bc;

	bc = ctx;
	while (num -- > 0) {
		bc->sig_len = HAWK_SIG_COMPACT_MAXSIZE(bc->logn);
		CC(hawk_sign(&bc->rng,
			bc->sig, &bc->sig_len, HAWK_SIG_COMPACT,
			bc->esk,
			"data", 4, bc->tmp, bc->tmp_len));
	}
	return 0;
}

static int
bench_uncompressed_verify(void *ctx, unsigned long num)
{
	bench_context *bc;

	bc = ctx;
	while (num -- > 0) {
		CC(hawk_uncompressed_verify(
			bc->uncompsig, bc->uncompsig_len, HAWK_SIG_COMPACT,
			bc->pk, HAWK_PUBKEY_SIZE[bc->logn],
			"data", 4, bc->tmp, bc->tmp_len));
	}
	return 0;
}

static int
bench_verify(void *ctx, unsigned long num)
{
	bench_context *bc;

	bc = ctx;
	while (num -- > 0) {
		CC(hawk_verify(
			bc->sig, bc->sig_len, HAWK_SIG_COMPACT,
			bc->pk, HAWK_PUBKEY_SIZE[bc->logn],
			"data", 4, bc->tmp, bc->tmp_len));
	}
	return 0;
}

static void
test_speed_hawk(unsigned logn, double threshold)
{
	bench_context bc;

	printf("%4u:", 1u << logn);
	fflush(stdout);

	bc.logn = logn;

	if (shake256_init_prng_from_system(&bc.rng) != 0) {
		fprintf(stderr, "random seeding failed\n");
		exit(EXIT_FAILURE);
	}
	// shake256_init_prng_from_seed(&bc.rng, "test", 4);
	shake256_flip(&bc.rng);

	bc.sk = xmalloc(HAWK_SECKEY_SIZE(logn));
	bc.pk = xmalloc(HAWK_PUBKEY_SIZE[logn]);
	bc.esk = xmalloc(HAWK_EXPANDEDKEY_SIZE(logn));
	bc.uncompsig = xmalloc(HAWK_UNCOMPRESSED_SIG_COMPACT_MAXSIZE(logn));
	bc.uncompsig_len = 0;
	bc.sig = xmalloc(HAWK_SIG_COMPACT_MAXSIZE(logn));
	bc.sig_len = 0;

	// keygen
	bc.tmp = xmalloc(bc.tmp_len = HAWK_TMPSIZE_KEYGEN(logn));
	printf(" %8.2f", do_bench(&bench_keygen, &bc, threshold) / 1000000.0);
	fflush(stdout);

	printf(" %8.2f", do_bench(&bench_expand_seckey, &bc, threshold) / 1000.0);
	fflush(stdout);

	// sign
#ifdef HAWK_AVX
	bc.tmp = realloc(bc.tmp, bc.tmp_len = HAWK_TMPSIZE_UNCOMPRESSED_SIGN(logn));
#else
	bc.tmp = realloc(bc.tmp, bc.tmp_len = HAWK_TMPSIZE_UNCOMPRESSED_SIGN_NTT(logn));
#endif
	printf(" %8.2f", do_bench(&bench_uncompressed_sign, &bc, threshold) / 1000.0);
	fflush(stdout);
#ifdef HAWK_AVX
	bc.tmp = realloc(bc.tmp, bc.tmp_len = HAWK_TMPSIZE_SIGNDYN(logn));
#else
	bc.tmp = realloc(bc.tmp, bc.tmp_len = HAWK_TMPSIZE_SIGNDYN_NTT(logn));
#endif
	printf(" %8.2f", do_bench(&bench_sign_dyn, &bc, threshold) / 1000.0);
	fflush(stdout);
	bc.tmp = realloc(bc.tmp, bc.tmp_len = HAWK_TMPSIZE_SIGN(logn));
	printf(" %8.2f", do_bench(&bench_sign, &bc, threshold) / 1000.0);
	fflush(stdout);

	// verify
#ifdef HAWK_AVX
	bc.tmp = realloc(bc.tmp, bc.tmp_len = HAWK_TMPSIZE_UNCOMPRESSED_VERIFY(logn));
#else
	bc.tmp = realloc(bc.tmp, bc.tmp_len = HAWK_TMPSIZE_UNCOMPRESSED_VERIFY_NTT(logn));
#endif
	printf(" %8.2f", do_bench(&bench_uncompressed_verify, &bc, threshold) / 1000.0);
	fflush(stdout);
#ifdef HAWK_AVX
	bc.tmp = realloc(bc.tmp, bc.tmp_len = HAWK_TMPSIZE_VERIFY(logn));
#else
	bc.tmp = realloc(bc.tmp, bc.tmp_len = HAWK_TMPSIZE_VERIFY_NTT(logn));
#endif
	printf(" %8.2f", do_bench(&bench_verify, &bc, threshold) / 1000.0);
	fflush(stdout);

	printf("\n");
	fflush(stdout);

	xfree(bc.tmp);
	xfree(bc.pk);
	xfree(bc.sk);
	xfree(bc.esk);
	xfree(bc.uncompsig);
	xfree(bc.sig);
}

int
main(int argc, char *argv[])
{
	double threshold;

	if (argc < 2) {
		threshold = 2.0;
	} else if (argc == 2) {
		threshold = atof(argv[1]);
	} else {
		threshold = -1.0;
	}
	if (threshold <= 0.0 || threshold > 60.0) {
		fprintf(stderr,
"usage: speed [ threshold ]\n"
"'threshold' is the minimum time for a bench run, in seconds (must be\n"
"positive and less than 60).\n");
		exit(EXIT_FAILURE);
	}
	printf("time threshold = %.4f s\n", threshold);
	printf("kg = keygen, ek = expand private key\n");
	printf("su = uncompressed sign, sd = sign (without ek), sek = sign (with ek), vu = uncompressed verify, vv = verify\n");
	printf("keygen in milliseconds, other values in microseconds\n");
	printf("\n");
	printf("degree  kg(ms)   ek(us)   su(us)   sd(us)   sek(us)  vu(us)   vv(us)\n");
	fflush(stdout);

	for (unsigned logn = 1; logn <= 10; logn++) {
		test_speed_hawk(logn, logn <= 8 ? 0.01 : threshold);
	}

	return 0;
}
