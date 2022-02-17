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

void random_hash(int8_t *h, unsigned logn) {
	assert(RAND_MAX == INT_MAX); // rand() should generate 31 random bits
	int x = rand();
	size_t RAND_BITS = 31, rand_bits = RAND_BITS;
	for (size_t u = MKN(logn); u -- > 0; ) {
		if (rand_bits == 0) {
			x = rand();
			rand_bits = RAND_BITS;
		}
		h[u] = (x & 1);
		x >>= 1;
		rand_bits--;
	}
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

void measure_sign_speed(fpr isigma_kg, fpr isigma_sig, uint32_t bound)
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
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, isigma_kg, logn, b);

		random_hash(h, logn);
		Zf(complete_sign)(&sc, s0, s1, f, g, F, G, h, isigma_sig, bound, logn, b);

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

constexpr fpr sigma_kg  = { v: 1.425 };
constexpr fpr sigma_sig = { v: 1.292 };
constexpr fpr verif_margin = { v: 1.1 };

int main() {
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	uint32_t bound = fpr_floor(fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_of(n))));
	measure_sign_speed(fpr_inv(sigma_kg), fpr_inv(sigma_sig), bound);
	return 0;
}
