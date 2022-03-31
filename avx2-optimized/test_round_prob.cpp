#include <cassert>
#include <climits>
#include <cstdio>
#include <vector>
#include <mutex>
#include <thread>

// x86_64 specific:
#include<sys/time.h>

extern "C" {
	#ifndef restrict
		#define restrict
	#endif

	#include "inner.h"
}

typedef long long ll;

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

ll time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
constexpr size_t logn = 9, n = MKN(logn);

void output_poly(int16_t *x) {
	for (size_t u = 0; u < n; u++) {
		printf("%d ", x[u]);
	}
	printf("\n");
}

const int num_samples = 512;

void measure_sign_speed()
{
	uint8_t b[48 << logn];
	int8_t f[n], g[n], F[n], G[n], h[n];
	int16_t s0[n], s1[n];
	ll Q00[n];
	fpr q00[n], q10[n], q11[n];
	unsigned char seed[48];
	inner_shake256_context sc;
	fpr fs0[n], fs1[n];

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	double avg = 0.0;
	double exp_sq = 0.0;

	// One key generation
	for (int rep = 0; rep < num_samples; rep++) {
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b);

		inner_shake256_extract(&sc, h, sizeof h);
		Zf(complete_sign)(&sc, s0, s1, f, g, F, G, h, logn, b);

		// calculate error (s0 q00 + s1 q10) / q00
		for (size_t u = 0; u < n; u++) {
			fs0[u] = fpr_half(fpr_of(2*s0[u] + h[u]));
			fs1[u] = fpr_of(s1[u]);
		}
		Zf(FFT)(fs0, logn);
		Zf(FFT)(fs1, logn);
		Zf(poly_mul_fft)(fs0, q00, logn);
		Zf(poly_mul_fft)(fs1, q10, logn);
		Zf(poly_add)(fs0, fs1, logn);
		Zf(poly_div_autoadj_fft)(fs0, q00, logn);
		Zf(iFFT)(fs0, logn);
		// measure errors

		for (size_t u = 0; u < n; u++) {
			double err = fs0[u].v;
			avg += err;
			exp_sq += err*err;
		}

		Zf(iFFT)(q00, logn);
		for (size_t u = 0; u < n; u++) {
			Q00[u] += fpr_rint(q00[u]);
		}
	}

	printf("Averages:");
	for (size_t u = 0; u < n; u++) {
		printf(" %.1f", (double)Q00[u] / num_samples);
	}
	printf("\n");
	printf("Final ratio dstribution: %.8f +/- %.8f\n", avg / num_samples / n, sqrt(exp_sq / num_samples / n));

}

int main() {
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	measure_sign_speed();
	return 0;
}
