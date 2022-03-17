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

/* Should be moved to vrfy.c */
/* Requires space for 6 polynomials of size 2^logn */
int
Zf(verify_nearest_plane_tree)(const fpr *restrict tree,
	const int8_t *restrict h0, const int8_t *restrict h1,
	const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp)
{
	/*
	 * This works better than simple rounding.
	 * Reconstruct s0, by running Babai's NP algorithm with target
	 *     h0/2 + (h1/2 - s1) q10 / q00.
	 */

	size_t n, u;
	fpr *t0, *t1, *t2, *t3, trace;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;
	t2 = t1 + n;
	t3 = t2 + n;

	for (u = 0; u < n; u ++) {
		t0[u] = fpr_half(fpr_of(h0[u] & 1));
	}
	for (u = 0; u < n; u ++) {
		t1[u] = fpr_half(fpr_of((h1[u] & 1) - 2 * s1[u]));
	}

	Zf(FFT)(t0, logn);
	Zf(FFT)(t1, logn);
	memcpy(t2, t1, n * sizeof *t1);
	Zf(poly_mul_fft)(t2, q10, logn);
	Zf(poly_div_fft)(t2, q00, logn);
	Zf(poly_add)(t2, t0, logn); // t2 = h0/2 + (h1/2 - s1) q10/q00

	// Run Babai with target t2 and Gram-matrix q00.
	Zf(ffNearestPlane_tree)(t3, tree, t2, logn, t3 + n);

	Zf(poly_sub)(t0, t3, logn);
	Zf(poly_mulconst)(t0, fpr_two, logn);
	// t0 = 2 * (t0 - t3) = h0 - 2 s0
	Zf(poly_mulconst)(t1, fpr_two, logn);
	// t1 = h1 - 2 s1

	// Currently in memory: s0, s1, s1, s0 (in FFT representation)
	memcpy(t2, t1, n * sizeof *t0);
	memcpy(t3, t0, n * sizeof *t0);

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
	for (u = 0; u < n/2; u ++) {
		trace = fpr_add(trace, t1[u]);
	}

	/*
	 * Note: only n/2 embeddings are stored, because they come in pairs.
	 */
	trace = fpr_double(trace);
	/*
	 * Renormalize, so we get the norm of (s0, s1) w.r.t Q.
	 */
	trace = fpr_div(trace, fpr_of(n));

	/*
	 * Signature is valid if and only if
	 *     `Tr(s* Q s) / n (=Tr(x^* x)/n = sum_i x_i^2) <= bound`.
	 * Note: check whether the norm is actually storable in a uint32_t.
	 */
	return fpr_lt(fpr_zero, trace) && fpr_lt(trace, fpr_ptwo31m1)
		&& (uint32_t)fpr_rint(trace) <= bound;

}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
constexpr size_t logn = 9, n = MKN(logn);

void output_poly(int16_t *x)
{
	for (size_t u = 0; u < n; u++) {
		printf("%d ", x[u]);
	}
	printf("\n");
}

const int num_samples = 200'000;

uint8_t pregen_h[num_samples][n/4] = {};
int8_t pregen_h0[num_samples][n] = {}, pregen_h1[num_samples][n] = {};
int16_t pregen_s0[num_samples][n], pregen_s1[num_samples][n], pregen_s[num_samples][n];

void measure_sign_speed(uint32_t bound)
{
	// uint8_t b[42 << logn];
	uint8_t b[50 << logn];
	int8_t f[n], g[n], F[n], G[n];
	fpr exp_sk[EXPANDED_SECKEY_SIZE(logn)]; // if logn is not known at compile-time, take fixed value
	int16_t s0[n], s1[n];
	fpr q00[n], q10[n], q11[n];
	unsigned char seed[48];
	inner_shake256_context sc;
	timeval t0, t1;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	// One key generation
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b);
	Zf(expand_seckey)(exp_sk, f, g, F, logn);

	// memset(pregen_h, 0, num_samples * n / 4);
	for (int rep = 0; rep < num_samples; rep++) {
		random_hash(pregen_h0[rep], logn);
		// random_hash(pregen_h1[rep], logn);
		for (unsigned i = 0; i < n; i++) {
			pregen_h[rep][i/8] |= pregen_h0[rep][i] << (i % 8);
			pregen_h[rep][(n + i) / 8] |= pregen_h1[rep][i] << (i % 8);
		}
	}

	// Also pregenerate all signatures
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		Zf(complete_sign)(&sc, pregen_s0[rep], pregen_s1[rep], f, g, F, G, pregen_h0[rep], pregen_h1[rep], bound, logn, b);
	}
	gettimeofday(&t1, NULL);
	ll cs_us = time_diff(&t0, &t1);

	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		Zf(sign)(&sc, s1, f, g, F, G, pregen_h0[rep], pregen_h1[rep], bound, logn, b);
	}
	gettimeofday(&t1, NULL);
	ll  s_us = time_diff(&t0, &t1);

	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		// Zf(guaranteed_sign)(&sc, s1, f, g, F, G, q00, pregen_h0[rep], pregen_h1[rep], bound, logn, b);
		Zf(expand_seckey)(exp_sk, f, g, F, logn);
		Zf(fft_sign)(&sc, pregen_s[rep], exp_sk, pregen_h[rep], bound, logn, b);
	}
	gettimeofday(&t1, NULL);
	ll gs_us = time_diff(&t0, &t1);

	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		Zf(fft_sign)(&sc, pregen_s[rep], exp_sk, pregen_h[rep], bound, logn, b);
	}
	gettimeofday(&t1, NULL);
	ll fs_us = time_diff(&t0, &t1);

	printf("  complete_sign: %.1f us/sign\n", (double)cs_us / num_samples);
	printf("           sign: %.1f us/sign\n", (double) s_us / num_samples);
	printf("guaranteed_sign: %.1f us/sign\n", (double)gs_us / num_samples);
	printf("       fft_sign: %.1f us/sign\n", (double)fs_us / num_samples);

	// Run tests:
	int nr_cor = 0;
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++)
		nr_cor += Zf(complete_verify)(pregen_h0[rep], pregen_h1[rep], pregen_s0[rep], pregen_s1[rep], q00, q10, q11, bound, logn, b);
	gettimeofday(&t1, NULL);
	printf("# correct = %d\n", nr_cor);
	ll cv_us = time_diff(&t0, &t1);

	nr_cor = 0;
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		nr_cor += Zf(verify_simple_rounding)(pregen_h0[rep], pregen_h1[rep], s0, pregen_s[rep], q00, q10, q11, bound, logn, b);
	}
	gettimeofday(&t1, NULL);
	printf("# correct = %d\n", nr_cor);
	ll vsr_us = time_diff(&t0, &t1);

	nr_cor = 0;
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		nr_cor += Zf(verify_simple_rounding_fft)(pregen_h[rep], pregen_s[rep], q00, q10, q11, bound, logn, b);
	}
	gettimeofday(&t1, NULL);
	printf("# correct = %d\n", nr_cor);
	ll vfsr_us = time_diff(&t0, &t1);

	nr_cor = 0;
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		nr_cor += Zf(verify_nearest_plane)(pregen_h0[rep], pregen_h1[rep], pregen_s0[rep], pregen_s[rep], q00, q10, q11, bound, logn, b);
	}
	gettimeofday(&t1, NULL);
	ll vnp_us = time_diff(&t0, &t1);

	fpr tree[LDL_TREESIZE(logn)];
	// Construct the LDL tree
	Zf(ffLDL_fft)(tree, q00, logn, (fpr *)b);

	printf("# correct = %d\n", nr_cor);
	nr_cor = 0;
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		nr_cor += Zf(verify_nearest_plane_tree)(tree, pregen_h0[rep], pregen_h1[rep], pregen_s[rep], q00, q10, q11, bound, logn, b);
	}
	gettimeofday(&t1, NULL);
	ll vnpt_us = time_diff(&t0, &t1);
	printf("# correct = %d\n", nr_cor);

	printf("       complete_verify: %.1f us/vrfy\n", (double) cv_us / num_samples);
	printf("verify_simple_rounding: %.1f us/vrfy\n", (double) vsr_us / num_samples);
	printf("verify_simple_roundfft: %.1f us/vrfy\n", (double)vfsr_us / num_samples);
	printf("  verify_nearest_plane: %.1f us/vrfy\n", (double) vnp_us / num_samples);
	printf("verify_nearest_plane_t: %.1f us/vrfy\n", (double) vnpt_us / num_samples);
}

struct WorkerResult {
	ll iters, round_fail, fft_fail;
	ll sig_fail, sigsize, sqsigsize, minsigsz, maxsigsz;

	WorkerResult() : iters(0), round_fail(0), fft_fail(0),
		sig_fail(0), sigsize(0), sqsigsize(0), minsigsz(INT_MAX), maxsigsz(0) {}

	void combine(const WorkerResult &res) {
		iters += res.iters;
		round_fail += res.round_fail;
		fft_fail += res.fft_fail;
		sig_fail += res.sig_fail;
		sigsize += res.sigsize;
		sqsigsize += res.sqsigsize;
		minsigsz = std::min(minsigsz, res.minsigsz);
		maxsigsz = std::max(maxsigsz, res.maxsigsz);
	}
};

WorkerResult measure_signatures(uint32_t bound)
{
	uint8_t b[50 << logn];
	int8_t f[n], g[n], F[n], G[n], h0[n], h1[n];
	uint8_t h[n / 4];
	int16_t recs0[n], s1[n];
	fpr q00[n], q10[n], q11[n];
	unsigned char seed[48];
	inner_shake256_context sc;
	fpr exp_sk[EXPANDED_SECKEY_SIZE(logn)];

	const int num_reps = 128 * 1024;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	WorkerResult result;
	result.iters = num_reps;

	prng rng;
	Zf(prng_init)(&rng, &sc);

	for (int rep = 0; rep < num_reps; rep++) {
		// Generate key pair.
		if ((rep & 1023) == 0) {
			Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b);
			Zf(expand_seckey)(exp_sk, f, g, F, logn);
		}

		// make a signature of a random message...
		// random_hash(h0, logn);
		// random_hash(h1, logn);

		/* memset(h, 0, sizeof h);
		for (unsigned i = 0; i < n; i++) {
			h[i/8] |= h0[i] << (i % 8);
			h[(n + i) / 8] |= h1[i] << (i % 8);
		} */
		Zf(prng_get_bytes)(&rng, h, sizeof h);

		// Compute the signature.
		// Zf(complete_sign)(&sc, recs0, s1, f, g, F, G, h0, h1, bound, logn, b);
		// Zf(sign)(&sc, s1, f, g, h0, h1, bound, logn, b);
		// Zf(guaranteed_sign)(&sc, s1, f, g, F, G, q00, h0, h1, bound, logn, b);
		Zf(fft_sign)(&sc, s1, exp_sk, h, bound, logn, b);

		// if (!Zf(verify_simple_rounding)(h0, h1, recs0, s1, q00, q10, q11, bound, logn, b))
			// result.round_fail++;
		if (!Zf(verify_simple_rounding_fft)(h, s1, q00, q10, q11, bound, logn, b))
			result.fft_fail++;

		/* size_t sig_sz = Zf(encode_sig)(NULL, 0, s1, logn, 5);
		if (sig_sz == 0) {
			result.sig_fail++;
		} else {
			result.sigsize += sig_sz;
			result.sqsigsize += sig_sz * sig_sz;
			result.minsigsz = std::min(result.minsigsz, (ll) sig_sz);
			result.maxsigsz = std::max(result.maxsigsz, (ll) sig_sz);
		} */
	}
	return result;
}

WorkerResult tot;
std::mutex mx;

constexpr fpr sigma_kg  = { v: 1.425 };
constexpr fpr sigma_sig = { v: 1.292 };
constexpr fpr verif_margin = { v: 1.1 };

uint32_t bound;

void work()
{
	WorkerResult result = measure_signatures(bound);

	/* acquire mutex lock */ {
		std::lock_guard<std::mutex> guard(mx);
		tot.combine(result);
	}
}

int main() {
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	bound = fpr_floor(fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_of(n))));

	// measure_sign_speed(bound);
	// return 0;

	const int nthreads = 4;
	std::thread* pool[nthreads-1];
	for (int i = 0; i < nthreads-1; i++) pool[i] = new std::thread(work);
	work();
	for (int i = 0; i < nthreads-1; i++) pool[i]->join(), delete pool[i];

	printf("\n");
	printf("# Signatures signed:      %lld\n", tot.iters);
	printf("# Simple roundings failed: %lld\n", tot.round_fail);
	printf("# FFT    roundings failed: %lld\n", tot.fft_fail);
	printf("# Sig. coding failed:     %lld\n", tot.sig_fail);

	ll nsigs = tot.iters - tot.sig_fail;
	double avg_sz = (double)tot.sigsize / nsigs;
	double std_sz = (double)tot.sqsigsize / nsigs - avg_sz * avg_sz;
	printf("# Average sig. size = %.2f (± %.2f) and \\in [%lld, %lld]\n", avg_sz, std_sz, tot.minsigsz, tot.maxsigsz);

	return 0;
}
