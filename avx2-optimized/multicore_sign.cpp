#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <vector>
// C++ parallelization:
#include <atomic>
#include <thread>

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

struct WorkerResult {
	int num_signed;
	int num_invalid;
	int num_babai_fail;
	int num_sig_differ;

	WorkerResult() =default;
};

WorkerResult measure_signatures(fpr isigma_kg, fpr isigma_sig, fpr verif_bound) {
	uint8_t b[42 << logn];
	int8_t f[n], g[n], F[n], G[n], h[n];
	int16_t s0[n], reconstructed_s0[n], s1[n];
	fpr q00[n], q10[n], q11[n];
	unsigned char seed[48];
	inner_shake256_context sc;
	const int n_repetitions = 10000;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	WorkerResult result;
	result.num_signed = n_repetitions;

	for (int rep = 0; rep < n_repetitions; rep++) {
		// Generate key pair.
		lilipu_keygen(&sc, f, g, F, G, q00, q10, q11, logn, b, isigma_kg);

		// make a signature of a random message...
		randombytes((unsigned char *)h, sizeof h);

		// Compute the signature.
		lilipu_complete_sign(&sc, s0, s1, f, g, F, G, h, logn, isigma_sig, b);

		if (!lilipu_complete_verify(h, s0, s1, q00, q10, q11, logn, verif_bound, b))
			result.num_invalid++;

		if (!lilipu_verify(h, reconstructed_s0, s1, q00, q10, q11, logn, verif_bound, b))
			result.num_babai_fail++;
		else
			assert(lilipu_complete_verify(h, reconstructed_s0, s1, q00, q10, q11, logn, verif_bound, b));

		for (size_t u = 0; u < n; u++) {
			if (s0[u] != reconstructed_s0[u]) {
				result.num_sig_differ++;
				break;
			}
		}
	}
	return result;
}

std::atomic<int> tot_signed(0), tot_invalid(0), tot_babai_fail(0), tot_sig_differ(0);

constexpr fpr sigma_kg  = { v: 1.425 };
constexpr fpr sigma_sig = { v: 1.292 };
constexpr fpr verif_margin = { v: sigma_kg.v / sigma_sig.v };

fpr getverif_bound() {
	return fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_sqr(fpr_of(n))));
}

void work() {
	WorkerResult result = measure_signatures(fpr_inv(sigma_kg), fpr_inv(sigma_sig), getverif_bound());

	tot_signed += result.num_signed;
	tot_invalid += result.num_invalid;
	tot_babai_fail += result.num_babai_fail;
	tot_sig_differ += result.num_sig_differ;;
}


int8_t valid_sigma(fpr sigma) {
	return !fpr_lt(sigma, fpr_sigma_min[logn]) && fpr_lt(sigma, fpr_div(fpr_of(18205), fpr_of(10000)));
}

int main() {
	unsigned seed = time(NULL);

	// verif_bound = (verif_margin 2 \sigma_sig)^2 (2 n^2)
	// Note 2n * n, where 2n comes from the rank of lattice R^2 over ZZ,
	// and n comes from the ratio between coefficient embedding and canonical embedding.
	const fpr verif_bound = fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_sqr(fpr_of(n))));

#ifdef __AVX2__
	const fpr sigma_FALCON = fpr_sqrt(fpr_div(fpr_of(117*117*Q), fpr_of(100*100*2*n))); // 1.17 √(q/2n)
	printf("Sigmas: %.3f (lilipu) vs %.2f (falcon)\n", sigma_kg.v, sigma_FALCON.v);
	printf("Verif margin: %.2f, bound: %.2f\n", verif_margin.v, verif_bound.v);
#endif

	printf("Seed: %u\n", seed);
	srand(seed);
	assert(valid_sigma(sigma_kg) && valid_sigma(sigma_sig));

	int nthreads = 4;
	if (nthreads == 1) {
		work();
	} else {
		std::vector<std::thread*> pool(nthreads);
		for (int i = 0; i < nthreads; i++) {
			pool[i] = new std::thread(work);
		}

		for (int i = 0; i < nthreads; i++) {
			pool[i]->join();
		}
	}

	printf("# Signatures signed:      %d\n", static_cast<int>(tot_signed));
	printf("# Signatures invalid:     %d\n", static_cast<int>(tot_invalid));
	printf("# Babai roundings failed: %d\n", static_cast<int>(tot_babai_fail));
	printf("# Babai != original s0:   %d\n", static_cast<int>(tot_sig_differ));

	return 0;
}
