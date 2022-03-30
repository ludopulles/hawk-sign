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

void fpr_to_int16(int16_t *d, fpr *p, size_t logn) {
	for (size_t u = 0; u < MKN(logn); u++) {
		d[u] = fpr_rint(p[u]);
	}
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
#define MAX_LOGN (9)
#define MAX_N MKN(MAX_LOGN)

const int n_repetitions = 100;

long long pk_sum[MAX_LOGN+1] = {}, pk_sumsq[MAX_LOGN+1] = {}, pk_min[MAX_LOGN+1] = {}, pk_max[MAX_LOGN+1] = {};
long long sk_sum[MAX_LOGN+1] = {}, sk_sumsq[MAX_LOGN+1] = {}, sk_min[MAX_LOGN+1] = {}, sk_max[MAX_LOGN+1] = {};
long long q00_sum[MAX_LOGN+1] = {}, q00_sumsq[MAX_LOGN+1] = {}, q00_min[MAX_LOGN+1] = {}, q00_max[MAX_LOGN+1] = {};
long long q10_sum[MAX_LOGN+1] = {}, q10_sumsq[MAX_LOGN+1] = {}, q10_min[MAX_LOGN+1] = {}, q10_max[MAX_LOGN+1] = {};

void test_compress(size_t logn) {
	uint8_t b[48 << MAX_LOGN], encode_buf[10000];
	int8_t f[MAX_N], g[MAX_N], F[MAX_N], G[MAX_N];
	int8_t _f[MAX_N], _g[MAX_N], _F[MAX_N];
	fpr q00[MAX_N], q10[MAX_N], q11[MAX_N];
	int16_t q00n[MAX_N], q10n[MAX_N], _q00n[MAX_N], _q10n[MAX_N];
	unsigned char seed[48];
	inner_shake256_context sc;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	for (int i = 0; i < n_repetitions; i++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b);

		Zf(iFFT)(q00, logn);
		Zf(iFFT)(q10, logn);

		fpr_to_int16(q00n, q00, logn);
		fpr_to_int16(q10n, q10, logn);

		int pk_sz = Zf(encode_pubkey)((void *)&encode_buf, 10000, q00n, q10n, logn);
		assert(pk_sz != 0);
		int _pksz = Zf(decode_pubkey)(q00n, q10n, (void *)&encode_buf, 10000, logn);
		assert(_pksz == pk_sz);

		int sk_sz = Zf(encode_seckey)((void *)&encode_buf, 10000, f, g, F, logn);
		assert(sk_sz != 0);
		int _sksz = Zf(decode_seckey)(f, g, F, (void *)&encode_buf, 10000, logn);
		assert(_sksz == sk_sz);

	}
}

void measure_keygen(size_t logn) {
	uint8_t b[48 << MAX_LOGN];
	int8_t f[MAX_N], g[MAX_N], F[MAX_N], G[MAX_N];
	fpr q00[MAX_N], q10[MAX_N], q11[MAX_N];
	int16_t q00n[MAX_N], q10n[MAX_N];
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

	for (int i = 0; i < n_repetitions; i++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b);

		for (size_t u = 0; u < MKN(logn); u++) {
			sq_fg += f[u]*f[u] + g[u]*g[u];
			sq_FG += F[u]*F[u];
		}

		Zf(iFFT)(q00, logn);
		Zf(iFFT)(q10, logn);

		fpr_to_int16(q00n, q00, logn);
		fpr_to_int16(q10n, q10, logn);

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

		for  (size_t u = 0; u < MKN(logn); u++) {
			long long x = fpr_rint(q00[u]);
			long long y = fpr_rint(q10[u]);
			if (u != 0) {
				q00_sum[logn] += x;
				q00_sumsq[logn] += x*x;
			}
			q10_sum[logn] += y;
			q10_sumsq[logn] += y*y;
		}
	}

	gettimeofday(&t1, NULL);
	// double kg_duration = (double)time_diff(&t0, &t1) / n_repetitions; // (in us)
	// printf("Average time per keygen: %.3f ms\n", kg_duration / 1000.0);

	printf("Sigma_{seckey} ~ %.8f, %.8f\n",
		sqrt(sq_fg / (2.0 * n_repetitions * MKN(logn))),
		sqrt(sq_FG / n_repetitions / MKN(logn)));

	printf("sigma_{q00}(%d) = %.8f, sigma_{q10} = %.8f\n",
		(int) logn,
		sqrt(q00_sumsq[logn] / (double) n_repetitions / (MKN(logn) - 1)),
		sqrt(q10_sumsq[logn] / (double) n_repetitions / MKN(logn))
	);
}

int main() {
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	for (size_t logn = 1; logn <= 9; logn++) {
		test_compress(logn);
	}

	for (size_t logn = 1; logn <= 9; logn++) {
		pk_min[logn] = sk_min[logn] = 1000 * 1000;
		measure_keygen(logn);
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
