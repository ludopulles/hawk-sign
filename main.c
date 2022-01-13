#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// x86_64 specific:
#include<sys/time.h>

#include "inner.h"

/*
 * Compute degree N from logarithm 'logn'.
 */
#define MKN(logn)   ((size_t)1 << (logn))


// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
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

void benchmark_lilipu(fpr isigma_kg, fpr isigma_sig) {
	union {
		uint8_t b[28 * 512];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	int8_t f[512], g[512], F[512], G[512], h[512];
	fpr q00[512], q10[512], q11[512];
	int16_t sig[512];
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

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, tmp.b, isigma_kg);

	gettimeofday(&t1, NULL);
	printf("Key generation took %lld microseconds\n", time_diff(&t0, &t1));

	// =========================================================================
	// | Benchmark the signing of random messages                              |
	// =========================================================================
	gettimeofday(&t0, NULL);

	for (int rep = 0; rep < n_repetitions; rep++) {
		// make a signature of a random message...
		randombytes((unsigned char *)h, sizeof h);

		// Compute the signature.
		Zf(sign)(&sc, sig, f, g, h, logn, isigma_sig, tmp.b);
	}

	gettimeofday(&t1, NULL);
	double sign_ps = 1000000LL * n_repetitions / (double)time_diff(&t0, &t1);
	printf("Lilipu sign/s = %.1f\n", sign_ps);
}

void test_lilipu_valid_signature(fpr isigma_kg, fpr isigma_sig, fpr verif_bound) {
	union {
		uint8_t b[42 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	int8_t f[1024], g[1024], F[1024], G[1024], h[1024];
	fpr q00[1024], q10[1024], q11[1024];
	int16_t sig[1024], s0[1024];
	unsigned char seed[48];
	inner_shake256_context sc;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, tmp.b, isigma_kg);

	for (int rep = 0; rep < 1000; rep++) {
		// make a signature of a random message...
		randombytes((unsigned char *)h, sizeof h);
		// to_sage16("h", (int16_t *)h, logn);

		// Compute the signature.
		Zf(sign)(&sc, sig, f, g, h, logn, isigma_sig, tmp.b);
		// to_sage16("s1", sig, logn);

		assert(Zf(verify)(h, s0, sig, q00, q10, q11, logn, verif_bound, tmp.b));
		// randombytes((unsigned char *)sig, sizeof sig);
		for (size_t u = 0; u < n; u ++)
			sig[u] = 0;
		assert(!Zf(verify)(h, s0, sig, q00, q10, q11, logn, verif_bound, tmp.b));
	}

	printf("Valid signatures were signed.\n");
	printf("No simple forgeries were possible.\n");
}


// try to find a forgery, put answer in s1
int try_forge(int16_t *s0, int16_t *s1, const fpr *q00, const fpr *q10, const fpr *q11,
	const int8_t *hm, const fpr verif_bound, uint8_t *tmp)
{
	size_t u;
	int16_t s0p[1024], s1p[1024];
	fpr *t0, *t1, *t2, *t3, trace;

	t0 = (fpr *)tmp;
	t1 = t0 + n;
	t2 = t1 + n;
	t3 = t2 + n;

	for (u = 0; u < n; u ++)
		s1[u] = 0;

	if (Zf(verify)(hm, s0, s1, q00, q10, q11, logn, verif_bound, tmp))
		return 1; // we are done

	unsigned char seed[48];
	inner_shake256_context sc;
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	sampler_context spc;
	void *samp_ctx;
	spc.sigma_min = fpr_sigma_min[logn];
	Zf(prng_init)(&spc.p, &sc);
	samp_ctx = &spc;

	fpr sigma = fpr_div(fpr_of(28660196105754623LL), fpr_of(10000000000000000LL));
	while (1) {
		for (u = 0; u < n; u ++)
			s0p[u] = 2*Zf(sampler)(samp_ctx, fpr_half(fpr_of(hm[u] & 1)), sigma) - (hm[u] & 1);
		for (u = 0; u < n; u ++)
			s1p[u] = 2*Zf(sampler)(samp_ctx, fpr_half(fpr_of(hm[u] & 1)), sigma) - (hm[u] & 1);

		for (u = 0; u < n; u++) t0[u] = fpr_of(s0p[u]);
		Zf(FFT)(t0, logn);
		memcpy(t3, t0, n * sizeof *t0);

		for (u = 0; u < n; u++) t1[u] = fpr_of(s1p[u]);
		Zf(FFT)(t1, logn);
		memcpy(t2, t1, n * sizeof *t0);
		// Currently in memory: s0, s1, s1, s0 (in FFT representation)

		// Compute s0 q00 s0* + s0 q01 s1* + s1 q10 s0* + s1 q11 s1*
		Zf(poly_mulselfadj_fft)(t2, logn);
		Zf(poly_mulselfadj_fft)(t3, logn);
		Zf(poly_mul_autoadj_fft)(t2, q11, logn); // t2 = s1 q11 s1*
		Zf(poly_mul_autoadj_fft)(t3, q00, logn); // t3 = s0 q00 s0*
		Zf(poly_muladj_fft)(t1, t0, logn); // t1 = s1 s0*
		Zf(poly_mul_fft)(t1, q10, logn); // t1 = s1 q10 s0*

		Zf(poly_addselfadj_fft)(t1, logn); // t1 = s1 q10 s0* + s0 q01 s1*
		Zf(poly_add_autoadj_fft)(t1, t2, logn);
		Zf(poly_add_autoadj_fft)(t1, t3, logn);

		trace = fpr_zero;
		for (u = 0; u < n/2; u ++)
			trace = fpr_add(trace, t1[u]);

		// note: only n/2 embeddings are stored,
		// the others are simply the conjugate embeddings.
		// TODO: this can be optimized in the verif_bound, cancelling with 2 in (2d).
		trace = fpr_double(trace);

		// Signature is valid if and only if `v` is short enough and s2%2 == 0 (<=> s1mod2 = 0).
		if (fpr_lt(trace, verif_bound)) {
			memcpy(s1, s1p, sizeof s1p);
			return 1; // we are done
		}
#ifdef AVX2
		printf("Nope, trace = %.2f\tvs %.2f.\n", trace.v, verif_bound.v);
		printf("s1p = ");
		for (u = 0; u < n; u++) printf("%d,", s1p[u]);
		printf("\n");
#endif
	}
}

void test_forge_signature(fpr isigma_kg, fpr isigma_sig, fpr verif_bound) {
	union {
		uint8_t b[42 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	int8_t f[1024], g[1024], F[1024], G[1024], h[1024];
	fpr q00[1024], q10[1024], q11[1024];
	int16_t sig[1024], s0[1024];
	unsigned char seed[48];
	inner_shake256_context sc;
	int res;


	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, tmp.b, isigma_kg);

	// generate random message
	randombytes((unsigned char *)h, sizeof h);
	// to_sage16("h", (int16_t *)h, logn);

	Zf(sign)(&sc, sig, f, g, h, logn, isigma_sig, tmp.b);
	printf("s1 = ");
	for (size_t u = 0; u < n; u++) printf("%d,", sig[u]);
	printf("\n");
	// to_sage16("s1", sig, logn);

	// try to forge a signature for 'h', and put result in s0
	res = try_forge(s0, sig, q00, q10, q11, h, verif_bound, tmp.b);
	assert(res);

	assert(!Zf(verify)(h, s0, sig, q00, q10, q11, logn, verif_bound, tmp.b));
	printf("No forgery found.\n");
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
	const fpr verif_margin = fpr_add(fpr_one, fpr_sqrt(fpr_mul(fpr_log2, fpr_div(fpr_of(64), fpr_of(n)))));
	// verif_bound = (verif_margin * 2 * sigma_sig)^2 * (2*d) * d
	// Here, the vector (x0, x1) \in Z^{2d} is sampled from a Discrete Gaussian with sigma equal to 2*sigma_sig
	// and lattice coset (h%2) + 2Z^{2d}, so it has a SQUARED norm of around ~(2sigma_sig)^2 * 2d.
	// Using trace(s Q s^H) = trace(x x^H) = ||x||^2 [K:\QQ] = ||x||^2 d, we arrive at the verif_bound.
	const fpr verif_bound = fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_sqr(fpr_of(n))));

#ifdef __AVX2__
	const fpr sigma_FALCON = fpr_sqrt(fpr_div(fpr_of(117*117*12289), fpr_of(100*100*2*n))); // 1.17 √(q/2n)
	printf("Sigmas: %.3f (lilipu) vs %.2f (falcon)\n", sigma_kg.v, sigma_FALCON.v);
	printf("Verif margin: %.2f, bound: %.2f\n", verif_margin.v, verif_bound.v);
#else
#endif

	// sigmas used in FALCON:
	// for (int i = 1; i <= 10; i++) printf("%d: %.2f\n", i, 1.17 * sqrt(Q / (2 << i)));

	printf("Seed: %u\n", seed);
	srand(seed);
	assert(valid_sigma(sigma_kg) && valid_sigma(sigma_sig));

	// test_forge_signature(fpr_inv(sigma_kg), fpr_inv(sigma_sig), verif_bound);
	benchmark_lilipu(fpr_inv(sigma_kg), fpr_inv(sigma_sig));

	test_lilipu_valid_signature(fpr_inv(sigma_kg), fpr_inv(sigma_sig), verif_bound);
	return 0;
}
