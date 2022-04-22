// c99 -W -Wall -O3 -march=native test_sizes.c -o test_sizes build/* -lm
#include <assert.h>
#include <stdio.h>
#include <stdint.h>
// x86_64 specific:
#include <sys/time.h>

#include "inner.h"

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
#define MAX_LOGN (9)
#define MAX_N MKN(MAX_LOGN)

long long pk_sum[MAX_LOGN+1] = {}, pk_sumsq[MAX_LOGN+1] = {}, pk_min[MAX_LOGN+1] = {}, pk_max[MAX_LOGN+1] = {};
long long sk_sum[MAX_LOGN+1] = {}, sk_sumsq[MAX_LOGN+1] = {}, sk_min[MAX_LOGN+1] = {}, sk_max[MAX_LOGN+1] = {};
long long q00_sum[MAX_LOGN+1] = {}, q00_sumsq[MAX_LOGN+1] = {};
long long q10_sum[MAX_LOGN+1] = {}, q10_sumsq[MAX_LOGN+1] = {};
long long q11_sum[MAX_LOGN+1] = {}, q11_sumsq[MAX_LOGN+1] = {};

void measure_keygen(size_t n_repetitions, size_t logn) {
	uint8_t b[48 << MAX_LOGN];
	int8_t f[MAX_N], g[MAX_N], F[MAX_N], G[MAX_N];
	fpr q00[MAX_N], q10[MAX_N], q11[MAX_N];
	int16_t q00n[MAX_N], q10n[MAX_N], q11n[MAX_N];
	unsigned char seed[48];
	inner_shake256_context sc;

	struct timeval t0, t1;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);

	long long sq_fg = 0, sq_FG = 0;

	double avg = 0.0, std = 0.0;

	for (int i = 0; i < n_repetitions; i++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b);

		for (size_t u = 0; u < MKN(logn); u++) {
			sq_fg += f[u]*f[u] + g[u]*g[u];
			sq_FG += F[u]*F[u] + G[u]*G[u];
		}

		Zf(fft_to_int16)(q00n, q00, logn);
		Zf(fft_to_int16)(q10n, q10, logn);
		Zf(fft_to_int16)(q11n, q11, logn);

		int pk_sz = Zf(encode_pubkey)(NULL, 0, q00n, q10n, logn);
		int sk_sz = Zf(encode_seckey)(NULL, 0, f, g, F, logn);

		pk_sum[logn] += pk_sz;
		pk_sumsq[logn] += pk_sz * pk_sz;
		if (pk_sz < pk_min[logn]) pk_min[logn] = pk_sz;
		if (pk_sz > pk_max[logn]) pk_max[logn] = pk_sz;

		sk_sum[logn] += sk_sz;
		sk_sumsq[logn] += sk_sz * sk_sz;
		if (sk_sz < sk_min[logn]) sk_min[logn] = sk_sz;
		if (sk_sz > sk_max[logn]) sk_max[logn] = sk_sz;

		double t = ( *(double*)&q00[0] ) / (0.25*Zf(l2bound)[logn]);
		avg += t;
		std += t*t;

		for  (size_t u = 0; u < MKN(logn); u++) {
			long long x = fpr_rint(q00[u]);
			long long y = fpr_rint(q10[u]);
			long long z = fpr_rint(q11[u]);
			if (u != 0) {
				q00_sum[logn] += x;
				q00_sumsq[logn] += x*x;
				q11_sum[logn] += z;
				q11_sumsq[logn] += z*z;
			}
			q10_sum[logn] += y;
			q10_sumsq[logn] += y*y;
		}
	}

	gettimeofday(&t1, NULL);
	// double kg_duration = (double)time_diff(&t0, &t1) / n_repetitions; // (in us)
	// printf("Average time per keygen: %.3f ms\n", kg_duration / 1000.0);

	printf("logn = %d: sig_{f,g} ~ %.8f, sig_{F} ~ %.8f, sig_{q00} ~ %.8f, sig_{q10} ~ %.8f, sig_{q11} ~ %.8f\n",
		(int) logn,
		sqrt(sq_fg / (2.0 * n_repetitions * MKN(logn))),
		sqrt(sq_FG / (2.0 * n_repetitions * MKN(logn))),
		sqrt(q00_sumsq[logn] / (double) n_repetitions / (MKN(logn) - 1)),
		sqrt(q10_sumsq[logn] / (double) n_repetitions / MKN(logn)),
		sqrt(q11_sumsq[logn] / (double) n_repetitions / (MKN(logn) - 1))
	);
	printf("averages Q: %.8f %.8f %.8f\n",
		q00_sum[logn] / (double) n_repetitions / (MKN(logn) - 1),
		q10_sum[logn] / (double) n_repetitions / MKN(logn),
		q11_sum[logn] / (double) n_repetitions / (MKN(logn) - 1)
	);
	fflush(stdout);

	avg /= n_repetitions;
	std /= n_repetitions;
	printf("Q00 ~ %.8f +/- %.8f\n", avg,
		sqrt( std - avg*avg ));
}

int main() {
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	size_t n_repetitions = 1000;
	for (size_t logn = 1; logn <= 9; logn++) {
		pk_min[logn] = sk_min[logn] = 1000 * 1000;
		measure_keygen(n_repetitions, logn);
	}

	printf("logn | Average +/- stddev | min  | max \n");
	double avg, std;

	printf("-----+- Secret key -------+------+-----\n");
	// Secret key
	for (size_t logn = 1; logn <= 9; logn++) {
		avg = (double) sk_sum[logn] / n_repetitions;
		std = sqrt( (double) sk_sumsq[logn] / n_repetitions - avg*avg );
		printf("%4d | %7.2f +/- %6.2f | %4lld | %4lld\n", (int) logn, avg, std, sk_min[logn], sk_max[logn]);
	}

	printf("-----+- Public key -------+------+-----\n");
	// Public key
	for (size_t logn = 1; logn <= 9; logn++) {
		avg = (double) pk_sum[logn] / n_repetitions;
		std = sqrt( (double) pk_sumsq[logn] / n_repetitions - avg*avg );
		printf("%4d | %7.2f +/- %6.2f | %4lld | %4lld\n", (int) logn, avg, std, pk_min[logn], pk_max[logn]);
	}
	return 0;
}
