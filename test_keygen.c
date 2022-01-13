#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// x86_64 specific:
#include<sys/time.h>

#include "keygen.c"

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

/* see inner.h for keygen */
void
keygen_count_fails(inner_shake256_context *rng,
	int8_t *restrict f, int8_t *restrict g, // secret key
	int8_t *restrict F, int8_t *restrict G, // secret key
	fpr *restrict q00, fpr *restrict q10, fpr *restrict q11, // public key
	unsigned logn, uint8_t *restrict tmp, fpr isigma_kg, int *num_fails)
{
	size_t n = MKN(logn);
	for (;;) {
		fpr *rt1, *rt2;
		sampler_context spc;
		void *samp_ctx;
		spc.sigma_min = fpr_sigma_min[logn];
		Zf(prng_init)(&spc.p, rng);
		samp_ctx = &spc;
		poly_small_mkgauss(samp_ctx, f, logn, isigma_kg, 128);
		poly_small_mkgauss(samp_ctx, g, logn, isigma_kg, 128);
		if (!solve_NTRU(logn, F, G, f, g, 128, (uint32_t *)tmp)) {
			(*num_fails)++;
			continue;
		}

		rt1 = (fpr *)tmp;
		rt2 = rt1 + n;
		poly_small_to_fp(q00, f, logn);
		poly_small_to_fp(rt1, g, logn);
		poly_small_to_fp(q11, F, logn);
		poly_small_to_fp(rt2, G, logn);
		Zf(FFT)(q00, logn); // f
		Zf(FFT)(rt1, logn); // g
		Zf(FFT)(q11, logn); // F
		Zf(FFT)(rt2, logn); // G
		Zf(poly_add_muladj_fft)(q10, q11, rt2, q00, rt1, logn);
		Zf(poly_mulselfadj_fft)(q00, logn); // f*adj(f)
		Zf(poly_mulselfadj_fft)(rt1, logn); // g*adj(g)
		Zf(poly_add)(q00, rt1, logn);
		Zf(poly_mulselfadj_fft)(q11, logn); // F*adj(F)
		Zf(poly_mulselfadj_fft)(rt2, logn); // G*adj(G)
		Zf(poly_add)(q11, rt2, logn);
		break;
	}
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
const size_t logn = 9, n = MKN(logn);

void measure_keygen(fpr isigma_kg) {
	uint8_t b[28 << logn]; // 14 kB temporary memory, 17.5 kB total
	int8_t f[n], g[n], F[n], G[n];
	fpr q00[n], q10[n], q11[n];
	unsigned char seed[48];
	inner_shake256_context sc;

	struct timeval t0, t1;
	const int n_repetitions = 1000;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);

	int fails = 0;
	for (int i = 0; i < n_repetitions; i++) {
		// Generate key pair.
		keygen_count_fails(&sc, f, g, F, G, q00, q10, q11, logn, b, isigma_kg, &fails);
	}

	gettimeofday(&t1, NULL);
	double kg_duration = (double)time_diff(&t0, &t1) / n_repetitions; // (in us)
	printf("Average time per keygen: %.3f ms\n", kg_duration / 1000.0);
	// This requires catching failed attempts
	printf("Probability failure: %.2f%%\n", 100.0 * fails / n_repetitions);
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

	int histogram[10000] = {};

	for (int rep = 0; rep < n_repetitions; rep++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b, isigma_kg);

		// make a signature of a random message...
		randombytes((unsigned char *)h, sizeof h);

		// Compute the signature.
		Zf(complete_sign)(&sc, s0, s1, f, g, F, G, h, logn, isigma_sig, b);

		for (size_t u = 0; u < n; u++)
			histogram[5000 + s1[u]]++;

		if (!Zf(complete_verify)(h, s0, s1, q00, q10, q11, logn, verif_bound, b)) {
			fprintf(stderr, "Invalid signature generated!\n");
		} else {
			if (!Zf(verify)(h, reconstructed_s0, s1, q00, q10, q11, logn, verif_bound, b)) {
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
	for (int i = 0; i < 10000; i++) {
		int freq = histogram[i];
		if (freq != 0) {
			printf("(%d,%d),", i - 5000, freq);
		}
	}
	printf("\n");
}

int8_t valid_sigma(fpr sigma_sig) {
	return !fpr_lt(sigma_sig, fpr_sigma_min[logn])
		&& fpr_lt(sigma_sig, fpr_div(fpr_of(18205), fpr_of(10000)));
}

int main() {
	unsigned seed = time(NULL);
	const fpr sigma_kg  = fpr_div(fpr_of(1425), fpr_of(1000));
	const fpr sigma_sig = fpr_div(fpr_of(1292), fpr_of(1000));
	const fpr verif_margin = fpr_div(sigma_kg, sigma_sig);
	printf("Seed: %u\n", seed);
	srand(seed);
	assert(valid_sigma(sigma_kg) && valid_sigma(sigma_sig));

	fpr isigma_kg = fpr_inv(sigma_kg), isigma_sig = fpr_inv(sigma_sig);
	fpr verif_bound = fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_sqr(fpr_of(n))));

	measure_keygen(isigma_kg);
	measure_signatures(isigma_kg, isigma_sig, verif_bound);
	return 0;
}
