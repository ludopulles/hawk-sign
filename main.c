#include <assert.h>
#include <limits.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// x86_64 specific:
#include<sys/time.h>

#include "inner.h"

// Simple randomness generator:
long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
const size_t logn = 9, n = MKN(logn);

void benchmark() {
	union {
		uint8_t b[50 * 512];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	int8_t f[512], g[512], F[512], G[512];
	uint8_t h[512 / 4];
	fpr q00[512], q10[512], q11[512];
	fpr exp_sk[EXPANDED_SECKEY_SIZE(logn)]; // if logn is not known at compile-time, take fixed value
	int16_t sig[512];
	unsigned char seed[48];
	inner_shake256_context sc;

	struct timeval t0, t1;
	const int n_repetitions = 10000;

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, tmp.b);
	Zf(expand_seckey)(exp_sk, f, g, F, logn);

	gettimeofday(&t1, NULL);
	printf("Key generation took %lld microseconds\n", time_diff(&t0, &t1));

	// =========================================================================
	// | Benchmark the signing of random messages                              |
	// =========================================================================
	gettimeofday(&t0, NULL);

	for (int rep = 0; rep < n_repetitions; rep++) {
		// make a signature of a random message...
		inner_shake256_extract(&sc, h, sizeof h);

		// Compute the signature.
		Zf(sign)(&sc, sig, exp_sk, h, logn, tmp.b);
	}

	gettimeofday(&t1, NULL);
	double sign_ps = 1000000LL * n_repetitions / (double)time_diff(&t0, &t1);
	printf("Lilipu sign/s = %.1f\n", sign_ps);
}

void test_valid_signature() {
	union {
		uint8_t b[50 * 512];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	int8_t f[512], g[512], F[512], G[512];
	uint8_t h[512/4];
	fpr q00[512], q10[512], q11[512];
	int16_t sig[512];
	unsigned char seed[48];
	inner_shake256_context sc;
	fpr exp_sk[EXPANDED_SECKEY_SIZE(logn)]; // if logn is not known at compile-time, take fixed value

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, tmp.b);
	Zf(expand_seckey)(exp_sk, f, g, F, logn);

	for (int rep = 0; rep < 1000; rep++) {
		// make a signature of a random message...
		inner_shake256_extract(&sc, h, sizeof h);

		// Compute the signature.
		Zf(sign)(&sc, sig, exp_sk, h, logn, tmp.b);

		assert(Zf(verify_simple_rounding_fft)(h, sig, q00, q10, q11, logn, tmp.b));
		for (size_t u = 0; u < n; u ++)
			sig[u] = 0;
		assert(!Zf(verify_simple_rounding_fft)(h, sig, q00, q10, q11, logn, tmp.b));
	}

	printf("Valid signatures were signed.\n");
	printf("No simple forgeries were possible.\n");
}

int main() {
	// const fpr sigma_kg  = fpr_div(fpr_of(1425), fpr_of(1000));
	// const fpr sigma_sig = fpr_div(fpr_of(1292), fpr_of(1000));
	// verif_margin ~ 1 + √(64 * ln(2) / N)   (see scheme.sage)
	// const fpr verif_margin = fpr_add(fpr_one, fpr_sqrt(fpr_mul(fpr_log2, fpr_div(fpr_of(64), fpr_of(n)))));
	// verif_bound = (verif_margin * 2 * sigma_sig)^2 * (2*d) * d
	// Here, the vector (x0, x1) \in Z^{2d} is sampled from a Discrete Gaussian with sigma equal to 2*sigma_sig
	// and lattice coset (h%2) + 2Z^{2d}, so it has a SQUARED norm of around ~(2sigma_sig)^2 * 2d.
	// Using trace(s Q s^H) = trace(x x^H) = ||x||^2 [K:\QQ] = ||x||^2 d, we arrive at the verif_bound.

	// verif_margin = sigma_kg/sigma_sig:
	// fail probability < 0.00003
	// $ ./python
	//     >>> r = 1.425 / 1.292
	//     >>> (1 + eps) * r**(2*d) * math.exp(-d*(r*r-1))
	//     2.7387274647723652e-05

	// sigmas used in FALCON:
	// for (int i = 1; i <= 10; i++) printf("%d: %.2f\n", i, 1.17 * sqrt(Q / (2 << i)));

	benchmark();
	test_valid_signature();
	return 0;
}
