#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// x86_64 specific:
#include<sys/time.h>

// Let C-code work with the C++ compiler:
#ifndef restrict
#define restrict   __restrict
#endif

#include "codec.c"
#include "common.c"
#include "fft.c"
#include "fpr.c"
#include "rng.c"
#include "shake.c"
#include "lilipu_keygen.c"

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
constexpr size_t logn = 9, n = MKN(logn);

void measure_keygen(fpr isigma_kg) {
	uint8_t b[28 * n]; // 14 kB temporary memory, 17.5 kB total
	int8_t f[n], g[n], F[n], G[n];
	fpr q00[n], q10[n], q11[n];
	unsigned char seed[48];
	inner_shake256_context sc;

	struct timeval t0, t1;
	const int n_repetitions = 10000;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);

	for (int i = 0; i < n_repetitions; i++) {
		// Generate key pair.
		lilipu_keygen(&sc, f, g, F, G, q00, q10, q11, logn, b, isigma_kg);
	}

	gettimeofday(&t1, NULL);
	double kg_duration = (double)time_diff(&t0, &t1) / n_repetitions; // (in us)
	double kg_ps = 1000000.0 / kg_duration;
	printf("Lilipu keygen/s = %.1f\n", kg_ps);
	printf("Average time (ms): %.3f\n", kg_duration / 1000.0);
	// printf("# of fails: %d (%.2f%%)\n", solve_NTRU_fails, 100.0 * solve_NTRU_fails / n_repetitions);
}

void measure_signatures(fpr isigma_kg, fpr isigma_sig, fpr verif_bound) {
	uint8_t b[42 << logn];
	int8_t f[n], g[n], F[n], G[n], h[n];
	int16_t s0[n], reconstructed_s0[n], s1[n];
	fpr q00[n], q10[n], q11[n];
	unsigned char seed[48];
	inner_shake256_context sc;

	const int n_repetitions = 500;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	for (int rep = 0; rep < n_repetitions; rep++) {
		// Generate key pair.
		lilipu_keygen(&sc, f, g, F, G, q00, q10, q11, logn, b, isigma_kg);

		// make a signature of a random message...
		randombytes((unsigned char *)h, sizeof h);

		// Compute the signature.
		lilipu_complete_sign(&sc, s0, s1, f, g, F, G, h, logn, isigma_sig, b);

		int sum0 = 0, sum1 = 0;
		for (size_t u = 0; u < n; u++) {
			// printf("%d,", s0[u]);
			sum0 += s0[u];
		}
		// printf("\n");
		for (size_t u = 0; u < n; u++) {
			// printf("%d,", s1[u]);
			sum1 += s1[u];
		}
		// printf("\n");
		printf("(%d,%d) ", sum0, sum1);

		if (!lilipu_verify(h, s0, s1, q00, q10, q11, logn, verif_bound, b)) {
			fprintf(stderr, "Invalid signature generated!\n");
		} else {
			if (!lilipu_verify(h, reconstructed_s0, s1, q00, q10, q11, logn, verif_bound, b)) {
				fprintf(stderr, "Babai was not succesful!\n");
			} else {
				int s0_eq = 1;
				for (size_t u = 0; u < n; u++)
					s0_eq &= (s0[u] == reconstructed_s0[u]);
				if (!s0_eq)
					fprintf(stderr, "Reconstructed s0 was different\n");
			}
		}
	}
	printf("All signatures were verified\n");
}

int8_t valid_sigma(fpr sigma_sig) {
	return !fpr_lt(sigma_sig, fpr_sigma_min[logn])
		&& fpr_lt(sigma_sig, fpr_div(fpr_of(18205), fpr_of(10000)));
}

int main() {
	unsigned seed = time(NULL);

	const fpr sigma_kg  = fpr_div(fpr_of(1425), fpr_of(1000));
	const fpr sigma_sig = fpr_div(fpr_of(1292), fpr_of(1000));
	// verif_margin = 1 + √(64 * ln(2) / 1024)   (see scheme.sage)
	const fpr verif_margin = fpr_add(fpr_one, fpr_sqrt(fpr_mul(fpr_log2, fpr_div(fpr_of(2), fpr_of(n)))));
	// verif_bound = (verif_margin * 2 * sigma_sig)^2 * (2*d) * d
	// Here, the vector (x0, x1) \in Z^{2d} is sampled from a Discrete Gaussian with sigma equal to 2*sigma_sig
	// and lattice coset (h%2) + 2Z^{2d}, so it has a SQUARED norm of around ~(2sigma_sig)^2 * 2d.
	// Using trace(s Q s^H) = trace(x x^H) = ||x||^2 [K:\QQ] = ||x||^2 d, we arrive at the verif_bound.
	const fpr verif_bound = fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_sqr(fpr_of(n))));

#ifdef __AVX2__
	const fpr sigma_FALCON = fpr_sqrt(fpr_div(fpr_of(117*117*Q), fpr_of(100*100*2*n))); // 1.17 √(q/2n)
	printf("Sigmas: %.3f (lilipu) vs %.2f (falcon)\n", sigma_kg.v, sigma_FALCON.v);
	printf("Verif margin: %.2f, bound: %.2f\n", verif_margin.v, verif_bound.v);
#endif

	printf("Seed: %u\n", seed);
	srand(seed);
	assert(valid_sigma(sigma_kg) && valid_sigma(sigma_sig));
	fpr isigma_kg = fpr_inv(sigma_kg), isigma_sig = fpr_inv(sigma_sig);

	// measure_keygen(isigma_kg);
	measure_signatures(isigma_kg, isigma_sig, verif_bound);
	return 0;
}
