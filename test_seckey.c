#include <math.h>
#include <stdio.h>
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

const size_t logn = 9, n = MKN(logn);

void measure_keygen(size_t logn) {
	uint8_t b[48 << MAX_LOGN];
	int8_t f[MAX_N], g[MAX_N], F[MAX_N], G[MAX_N];
	fpr q00[MAX_N], q10[MAX_N], q11[MAX_N];
	unsigned char seed[48];
	inner_shake256_context sc;

	struct timeval t0, t1;
	const int n_repetitions = 100;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);

	long long sum = 0, sumsq = 0, maxv = 0;
	for (int i = 0; i < n_repetitions; i++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b);

		int len = Zf(encode_seckey)(NULL, 0, f, g, F, logn);
		if (len > maxv)
			maxv = len;
		sum += len;
		sumsq += len*len;
	}

	gettimeofday(&t1, NULL);
	double kg_duration = (double)time_diff(&t0, &t1) / n_repetitions; // (in us)
	printf("Average time per keygen: %.3f ms\n", kg_duration / 1000.0);

	double avg = (double) sum / n_repetitions;
	double std = sqrt( (double) sumsq / n_repetitions - avg*avg );
	printf("Size of encode_seckey (logn = %zu) = %.2f +/- %.2f, max = %lld\n", logn, avg, std, maxv);
}

int main() {
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	for (size_t logn = 3; logn <= MAX_LOGN; logn++) {
		measure_keygen(logn);
	}
	return 0;
}
