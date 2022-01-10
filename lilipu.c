#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// x86_64 specific:
#include<sys/time.h>

// #include "api.h"
// #include "inner.h"
#include "codec.c"
#include "common.c"
#include "fft.c"
#include "fpr.c"
// #include "keygen.c"
// #include "nist.c"
#include "rng.c"
#include "shake.c"
#include "sign.c"
#include "lilipu_keygen.c"
// #include "vrfy.c"

// See nist.c for info on how to use TEMPALLOC
#define TEMPALLOC
// See nist.c
#define NONCELEN   40
// See api.h
#define CRYPTO_SECRETKEYBYTES   2305
#define CRYPTO_PUBLICKEYBYTES   1793
#define CRYPTO_BYTES            1330
#define CRYPTO_ALGNAME          "Lilipu-1024"

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
// | TODO: make this into an actually NIST-compatible function like `nist.c`   |
// =============================================================================
int
crypto_lilipu_sign(unsigned char *sm, unsigned long long *smlen,
	const unsigned char *m, unsigned long long mlen,
	const int8_t *restrict f, const int8_t *restrict g, fpr isigma_sig)
{
	TEMPALLOC union {
		uint8_t b[28 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	TEMPALLOC union {
		int16_t sig[1024];
		uint16_t hm[1024];
	} r;
	TEMPALLOC unsigned char seed[48], nonce[NONCELEN];
	TEMPALLOC unsigned char esig[CRYPTO_BYTES - 2 - sizeof nonce];
	TEMPALLOC inner_shake256_context sc;
	size_t sig_len;

	/*
	 * Create a random nonce (40 bytes).
	 */
	randombytes(nonce, sizeof nonce);

	/*
	 * Hash message nonce + message into a vector.
	 */
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, nonce, sizeof nonce);
	inner_shake256_inject(&sc, m, mlen);
	inner_shake256_flip(&sc);
	Zf(hash_to_point_vartime)(&sc, (uint16_t *)r.hm, 10);

	/*
	 * Initialize a RNG.
	 */
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);


	/*
	 * Compute the signature.
	 */
	lilipu_sign(r.sig, &sc, f, g, r.hm, 10, isigma_sig, tmp.b);


	/*
	 * Encode the signature and bundle it with the message. Format is:
	 *   signature length     2 bytes, big-endian
	 *   nonce                40 bytes
	 *   message              mlen bytes
	 *   signature            slen bytes
	 */
	esig[0] = 0x20 + 10;
	sig_len = Zf(comp_encode)(esig + 1, (sizeof esig) - 1, r.sig, 10);
	if (sig_len == 0) {
		return -1;
	}
	sig_len ++;
	memmove(sm + 2 + sizeof nonce, m, mlen);
	sm[0] = (unsigned char)(sig_len >> 8);
	sm[1] = (unsigned char)sig_len;
	memcpy(sm + 2, nonce, sizeof nonce);
	memcpy(sm + 2 + (sizeof nonce) + mlen, esig, sig_len);
	*smlen = 2 + (sizeof nonce) + mlen + sig_len;
	return 0;
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
const size_t logn = 9, n = MKN(logn);

void benchmark_lilipu(fpr isigma_kg, fpr isigma_sig) {
	TEMPALLOC union {
		// uint8_t b[42 * 1024];
		uint8_t b[100 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	TEMPALLOC int8_t f[1024], g[1024], F[1024], G[1024];
	TEMPALLOC fpr q00[1024], q10[1024], q11[1024];
	TEMPALLOC uint16_t h[1024];
	TEMPALLOC int16_t sig[1024];
	TEMPALLOC unsigned char seed[48];
	TEMPALLOC inner_shake256_context sc;

	struct timeval t0, t1;
	const int n_repetitions = 1000;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);

	// Generate key pair.
	lilipu_keygen(&sc, f, g, F, G, q00, q10, q11, logn, tmp.b, isigma_kg);

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
		lilipu_sign(sig, &sc, f, g, h, logn, isigma_sig, tmp.b);
	}

	gettimeofday(&t1, NULL);
	double sign_ps = 1000000LL * n_repetitions / (double)time_diff(&t0, &t1);
	printf("Lilipu sign/s = %.1f\n", sign_ps);
}

void benchmark_falcon() {
	TEMPALLOC union {
		uint8_t b[FALCON_KEYGEN_TEMP_10];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	TEMPALLOC union {
		uint8_t b[72 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp2;
	TEMPALLOC union {
		int16_t sig[1024];
		uint16_t hm[1024];
	} r;
	TEMPALLOC int8_t f[1024], g[1024], F[1024], G[1024];
	TEMPALLOC uint16_t h[1024];
	TEMPALLOC unsigned char seed[48];
	TEMPALLOC inner_shake256_context sc;

	struct timeval t0, t1;
	const int n_repetitions = 1000;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);

	// Generate key pair.
	falcon_inner_keygen(&sc, f, g, F, G, h, logn, tmp.b);

	gettimeofday(&t1, NULL);
	printf("Key generation took %lld microseconds\n", time_diff(&t0, &t1));

	// =========================================================================
	// | Start the signing                                                     |
	// =========================================================================
	gettimeofday(&t0, NULL);

	for (int rep = 0; rep < n_repetitions; rep++) {
		// make a signature of a random message...
		randombytes((unsigned char *)h, sizeof h);

		// Compute the signature.
		falcon_inner_sign_dyn(r.sig, &sc, f, g, F, G, r.hm, logn, tmp2.b);
	}

	gettimeofday(&t1, NULL);
	double sign_ps = 1000000LL * n_repetitions / (double)time_diff(&t0, &t1);
	printf("Falcon sign/s = %.1f\n", sign_ps);
}

void test_lilipu_valid_signature(fpr isigma_kg, fpr isigma_sig, fpr verif_bound) {
	TEMPALLOC union {
		uint8_t b[42 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	TEMPALLOC int8_t f[1024], g[1024], F[1024], G[1024];
	TEMPALLOC fpr q00[1024], q10[1024], q11[1024];
	TEMPALLOC int16_t sig[1024], s0[1024];
	TEMPALLOC uint16_t h[1024];
	TEMPALLOC unsigned char seed[48];
	TEMPALLOC inner_shake256_context sc;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	// Generate key pair.
	lilipu_keygen(&sc, f, g, F, G, q00, q10, q11, logn, tmp.b, isigma_kg);

	for (int rep = 0; rep < 1000; rep++) {
		// make a signature of a random message...
		randombytes((unsigned char *)h, sizeof h);
		// to_sage16("h", (int16_t *)h, logn);

		// Compute the signature.
		lilipu_sign(sig, &sc, f, g, h, logn, isigma_sig, tmp.b);
		// to_sage16("s1", sig, logn);

		assert(lilipu_verify(h, s0, sig, q00, q10, q11, logn, verif_bound, (fpr *)tmp.b));
		// randombytes((unsigned char *)sig, sizeof sig);
		for (size_t u = 0; u < n; u ++)
			sig[u] = 0;
		assert(!lilipu_verify(h, s0, sig, q00, q10, q11, logn, verif_bound, (fpr *)tmp.b));
	}

	printf("Valid signatures were signed.\n");
	printf("No simple forgeries were possible.\n");
}


// try to find a forgery, put answer in s1
int try_forge(int16_t *s0, int16_t *s1, const fpr *q00, const fpr *q10, const fpr *q11,
	const uint16_t *hm, const fpr verif_bound, fpr *tmp)
{
	size_t u;
	int16_t s0p[1024], s1p[1024];
	fpr *t0, *t1, *t2, *t3, trace;

	t0 = tmp;
	t1 = t0 + n;
	t2 = t1 + n;
	t3 = t2 + n;

	for (u = 0; u < n; u ++)
		s1[u] = 0;

	if (lilipu_verify(hm, s0, s1, q00, q10, q11, logn, verif_bound, tmp))
		return 1; // we are done

	unsigned char seed[48];
	inner_shake256_context sc;
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	sampler_context spc;
	samplerZ samp;
	void *samp_ctx;
	spc.sigma_min = fpr_sigma_min[logn];
	falcon_inner_prng_init(&spc.p, &sc);
	samp = Zf(sampler);
	samp_ctx = &spc;

	fpr sigma = fpr_div(fpr_of(28660196105754623LL), fpr_of(10000000000000000LL));
	while (1) {
		for (u = 0; u < n; u ++)
			s0p[u] = 2*samp(samp_ctx, fpr_half(fpr_of(hm[u] & 1)), sigma) - (hm[u] & 1);
		for (u = 0; u < n; u ++)
			s1p[u] = 2*mkgauss(&sc, logn);

		for (int rep = 10; rep-->0; ) {
			for (u = 0; u < n; u ++)
				s0p[u] += 2*mkgauss(&sc, logn);
			for (u = 0; u < n; u ++)
				s1p[u] += 2*mkgauss(&sc, logn);
		}

		for (u = 0; u < n; u++) t0[u] = fpr_of(s0p[u]);
		falcon_inner_FFT(t0, logn);
		memcpy(t3, t0, n * sizeof *t0);

		for (u = 0; u < n; u++) t1[u] = fpr_of(s1p[u]);
		falcon_inner_FFT(t1, logn);
		memcpy(t2, t1, n * sizeof *t0);
		// Currently in memory: s0, s1, s1, s0 (in FFT representation)

		// Compute s0 q00 s0* + s0 q01 s1* + s1 q10 s0* + s1 q11 s1*
		falcon_inner_poly_mulselfadj_fft(t2, logn);
		falcon_inner_poly_mulselfadj_fft(t3, logn);
		falcon_inner_poly_mul_autoadj_fft(t2, q11, logn); // t2 = s1 q11 s1*
		falcon_inner_poly_mul_autoadj_fft(t3, q00, logn); // t3 = s0 q00 s0*
		falcon_inner_poly_muladj_fft(t1, t0, logn); // t1 = s1 s0*
		falcon_inner_poly_mul_fft(t1, q10, logn); // t1 = s1 q10 s0*

		lilipu_inner_poly_addselfadj_fft(t1, logn); // t1 = s1 q10 s0* + s0 q01 s1*
		lilipu_inner_poly_add_autoadj_fft(t1, t2, logn);
		lilipu_inner_poly_add_autoadj_fft(t1, t3, logn);

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
	TEMPALLOC union {
		uint8_t b[42 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	TEMPALLOC int8_t f[1024], g[1024], F[1024], G[1024];
	TEMPALLOC fpr q00[1024], q10[1024], q11[1024];
	TEMPALLOC int16_t sig[1024], s0[1024];
	TEMPALLOC uint16_t h[1024];
	TEMPALLOC unsigned char seed[48];
	TEMPALLOC inner_shake256_context sc;
	int res;


	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	// Generate key pair.
	lilipu_keygen(&sc, f, g, F, G, q00, q10, q11, logn, tmp.b, isigma_kg);

	// generate random message
	randombytes((unsigned char *)h, sizeof h);
	// to_sage16("h", (int16_t *)h, logn);

	lilipu_sign(sig, &sc, f, g, h, logn, isigma_sig, tmp.b);
	printf("s1 = ");
	for (size_t u = 0; u < n; u++) printf("%d,", sig[u]);
	printf("\n");
	// to_sage16("s1", sig, logn);

	// try to forge a signature for 'h', and put result in s0
	res = try_forge(s0, sig, q00, q10, q11, h, verif_bound, (fpr *)tmp.b);
	assert(res);

	assert(!lilipu_verify(h, s0, sig, q00, q10, q11, logn, verif_bound, (fpr *)tmp.b));
	printf("No forgery found.\n");
}

void testmod2() {
	int8_t f[1024], fg[1024], check[1024];
	uint16_t g[1024], tmp[2*1024];

	for (size_t rep = 0; rep < 10; rep ++) {
		for (size_t u = 0; u < n; u ++)
			f[u] = rand();
		for (size_t u = 0; u < n; u ++)
			g[u] = rand();
		for (size_t u = 0; u < n; u ++)
			check[u] = 0;

		for (size_t u = 0; u < n; u ++)
			for (size_t v = 0; v < n; v ++)
				check[modp_add(u, v, n)] ^= f[u] & g[v];
		for (size_t u = 0; u < n; u ++)
			check[u] &= 1;

		lilipu_inner_mulmod2(fg, f, g, logn, tmp);

		for (size_t u = 0; u < n; u ++)
			assert(fg[u] == check[u]);
	}

/*
	struct timeval t0, t1;
	gettimeofday(&t0, NULL);
	for (size_t rep = 0; rep < 100; rep ++) {
		for (size_t u = 0; u < n; u ++) f[u] = rand();
		for (size_t u = 0; u < n; u ++) g[u] = rand();
		lilipu_inner_mulmod2_old(fg, f, g, logn, tmp);
	}
	gettimeofday(&t1, NULL);
	printf("d log d took: %d μs\n", time_diff(&t0, &t1));

	gettimeofday(&t0, NULL);
	for (size_t rep = 0; rep < 100; rep ++) {
		for (size_t u = 0; u < n; u ++) f[u] = rand();
		for (size_t u = 0; u < n; u ++) g[u] = rand();
		lilipu_inner_mulmod2(fg, f, g, logn);
	}
	gettimeofday(&t1, NULL);
	printf("d^2 took: %d μs\n", time_diff(&t0, &t1));
*/
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
	const fpr sigma_FALCON = fpr_sqrt(fpr_div(fpr_of(117*117*Q), fpr_of(100*100*2*n))); // 1.17 √(q/2n)
	printf("Sigmas: %.3f (lilipu) vs %.2f (falcon)\n", sigma_kg.v, sigma_FALCON.v);
	printf("Verif margin: %.2f, bound: %.2f\n", verif_margin.v, verif_bound.v);
#else
#endif

	// sigmas used in FALCON:
	// for (int i = 1; i <= 10; i++) printf("%d: %.2f\n", i, 1.17 * sqrt(Q / (2 << i)));

	printf("Seed: %u\n", seed);
	srand(seed);
	// prepare_random();
	testmod2();
	assert(valid_sigma(sigma_kg) && valid_sigma(sigma_sig));

	// test_forge_signature(fpr_inv(sigma_kg), fpr_inv(sigma_sig), verif_bound);
	benchmark_lilipu(fpr_inv(sigma_kg), fpr_inv(sigma_sig));
	benchmark_falcon();

	test_lilipu_valid_signature(fpr_inv(sigma_kg), fpr_inv(sigma_sig), verif_bound);
	return 0;
}
