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

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

void to_sage(const char *varname, const int8_t *f, unsigned logn) {
	printf("%s = %d", varname, f[0]);
	for (size_t u = 1; u < MKN(logn); u ++) {
		if (f[u] > 0) printf("+%d*z^%zu", f[u], u);
		if (f[u] < 0) printf("-%d*z^%zu", -f[u], u);
	}
	printf("\n");
}

void to_sage16(const char *varname, const int16_t *f, unsigned logn) {
	printf("%s = %d", varname, f[0]);
	for (size_t u = 1; u < MKN(logn); u ++) {
		if (f[u] > 0) printf("+%d*z^%zu", f[u], u);
		if (f[u] < 0) printf("-%d*z^%zu", -f[u], u);
	}
	printf("\n");
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
const size_t logn = 9, n = MKN(logn);

void benchmark(uint32_t bound) {
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
	randombytes(seed, sizeof seed);
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
		randombytes(h, n / 4);

		// Compute the signature.
		// Zf(sign)(&sc, sig, f, g, F, G, h0, h1, bound, logn, tmp.b);
		Zf(fft_sign)(&sc, sig, exp_sk, h, bound, logn, tmp.b);
	}

	gettimeofday(&t1, NULL);
	double sign_ps = 1000000LL * n_repetitions / (double)time_diff(&t0, &t1);
	printf("Lilipu sign/s = %.1f\n", sign_ps);
}

void test_valid_signature(uint32_t bound) {
	union {
		uint8_t b[50 * 512];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	int8_t f[512], g[512], F[512], G[512], h0[512], h1[512];
	fpr q00[512], q10[512], q11[512];
	int16_t sig[512], s0[512];
	unsigned char seed[48];
	inner_shake256_context sc;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, tmp.b);

	for (int rep = 0; rep < 1000; rep++) {
		// make a signature of a random message...
		random_hash(h0, logn);
		random_hash(h1, logn);

		// Compute the signature.
		Zf(sign)(&sc, sig, f, g, F, G, h0, h1, bound, logn, tmp.b);

		assert(Zf(verify_simple_rounding)(h0, h1, s0, sig, q00, q10, q11, bound, logn, tmp.b));
		for (size_t u = 0; u < n; u ++)
			sig[u] = 0;
		assert(!Zf(verify_simple_rounding)(h0, h1, s0, sig, q00, q10, q11, bound, logn, tmp.b));
	}

	printf("Valid signatures were signed.\n");
	printf("No simple forgeries were possible.\n");
}

int8_t valid_sigma(fpr sigma_sig) {
	return !fpr_lt(sigma_sig, fpr_sigma_min[logn])
		&& fpr_lt(sigma_sig, fpr_div(fpr_of(18205), fpr_of(10000)));
}

int main() {
	unsigned seed = time(NULL);

	const fpr sigma_kg  = fpr_div(fpr_of(1425), fpr_of(1000));
	const fpr sigma_sig = fpr_div(fpr_of(1292), fpr_of(1000));
	// verif_margin ~ 1 + √(64 * ln(2) / N)   (see scheme.sage)
	const fpr verif_margin = fpr_add(fpr_one, fpr_sqrt(fpr_mul(fpr_log2, fpr_div(fpr_of(64), fpr_of(n)))));
	// verif_bound = (verif_margin * 2 * sigma_sig)^2 * (2*d) * d
	// Here, the vector (x0, x1) \in Z^{2d} is sampled from a Discrete Gaussian with sigma equal to 2*sigma_sig
	// and lattice coset (h%2) + 2Z^{2d}, so it has a SQUARED norm of around ~(2sigma_sig)^2 * 2d.
	// Using trace(s Q s^H) = trace(x x^H) = ||x||^2 [K:\QQ] = ||x||^2 d, we arrive at the verif_bound.
	uint32_t bound = fpr_floor(fpr_mul(fpr_double(fpr_of(n)),
		fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig)))));

	// verif_margin = sigma_kg/sigma_sig:
	// fail probability < 0.00003
	// $ ./python
	//     >>> r = 1.425 / 1.292
	//     >>> (1 + eps) * r**(2*d) * math.exp(-d*(r*r-1))
	//     2.7387274647723652e-05

#ifdef __AVX2__
	const fpr sigma_FALCON = fpr_sqrt(fpr_div(fpr_of(117*117*12289), fpr_of(100*100*2*n))); // 1.17 √(q/2n)
	printf("Sigma: %.3f vs %.3f of falcon\n", sigma_kg.v, sigma_FALCON.v);
	printf("Verif margin: %.2f, bound: %u\n", verif_margin.v, bound);
#else
#endif

	// sigmas used in FALCON:
	// for (int i = 1; i <= 10; i++) printf("%d: %.2f\n", i, 1.17 * sqrt(Q / (2 << i)));

	printf("Seed: %u\n", seed);
	srand(seed);
	// assert(valid_sigma(sigma_kg) && valid_sigma(sigma_sig));

	benchmark(bound);
	test_valid_signature(bound);
	return 0;
}
