#include <cassert>
#include <climits>
#include <cstdio>
#include <vector>
#include <mutex>
#include <thread>

// x86_64 specific:
#include<sys/time.h>

extern "C" {
	#define restrict
	#include "../inner.h"
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

#define SECOND_HASH(h, logn) \
	((h) + ((logn) <= 3 ? 1u : 1u << ((logn) - 3)))

/**
 * If s != NULL, set the polynomial p equal to h - 2*s.
 * Otherwise, set p equal to h.
 * Returns the polynomial in FFT format.
 */
static void
hash_to_fft(fpr *p, const uint8_t *h, const int16_t *s, unsigned logn)
{
	size_t n, u, v;
	uint8_t hash;

	n = MKN(logn);
	if (logn <= 3) {
		hash = h[0];
		for (v = 0; v < n; v ++) {
			p[v] = fpr_of((hash & 1) - (s == NULL ? 0 : 2 * s[v]));
			hash >>= 1;
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				p[u] = fpr_of((hash & 1) - (s == NULL ? 0 : 2 * s[u]));
				hash >>= 1;
			}
		}
	}
	Zf(FFT)(p, logn);
}


/* Requires space for 6 polynomials of size 2^logn */
static int
verify_nearest_plane_tree(const fpr *restrict tree,
	const uint8_t *restrict h, int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	unsigned logn, uint8_t *restrict tmp)
{
	/*
	 * This works better than simple rounding.
	 * Reconstruct s0, by running Babai's NP algorithm with target
	 *     h0/2 + (h1/2 - s1) q10 / q00.
	 */

	size_t n;
	fpr *t0, *t1;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;

	hash_to_fft(t0, h, NULL, logn);
	hash_to_fft(t1, SECOND_HASH(h, logn), s1, logn);

	Zf(poly_mul_fft)(t1, q10, logn);
	Zf(poly_div_autoadj_fft)(t1, q00, logn);
	Zf(poly_add)(t0, t1, logn);
	Zf(poly_mulconst)(t0, fpr_onehalf, logn);

	// Run Babai with target t0 and Gram-matrix q00.
	Zf(ffNearestPlane_tree)(t1, tree, t0, logn, t1 + n);
	Zf(fft_to_int16)(s0, t1, logn);

	return Zf(uncompressed_verify)(h, s0, s1, q00, q10, q11, logn, tmp);
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
constexpr size_t logn = 10, n = MKN(logn);

const int num_samples = 1000;

uint8_t pregen_h[num_samples][n / 4] = {};
int16_t pregen_s0[num_samples][n], pregen_s1[num_samples][n], pregen_s[num_samples][n];

void measure_sign_speed()
{
	// uint8_t b[42 << logn];
	uint8_t b[50 << logn];
	int8_t f[n], g[n], F[n], G[n];
	fpr exp_sk[EXPANDED_SECKEY_SIZE(logn)]; // if logn is not known at compile-time, take fixed value
	int16_t iq00[n], iq10[n], s0[n], s1[n];
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
	Zf(keygen)(&sc, f, g, F, G, iq00, iq10, logn, b);
	Zf(complete_pubkey)(iq00, iq10, q00, q10, q11, logn);
	Zf(expand_seckey)(exp_sk, f, g, F, logn);

	// memset(pregen_h, 0, num_samples * n / 4);
	for (int rep = 0; rep < num_samples; rep++) {
		inner_shake256_extract(&sc, pregen_h[rep], n / 4);
	}

	// Also pregenerate all signatures
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		Zf(uncompressed_sign)(&sc, pregen_s0[rep], pregen_s1[rep], f, g, F, G, pregen_h[rep], logn, b);
	}
	gettimeofday(&t1, NULL);
	ll cs_us = time_diff(&t0, &t1);

	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		Zf(sign_dyn)(&sc, s1, f, g, F, G, pregen_h[rep], logn, b);
	}
	gettimeofday(&t1, NULL);
	ll sd_us = time_diff(&t0, &t1);

	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		Zf(sign)(&sc, pregen_s[rep], exp_sk, pregen_h[rep], logn, b);
	}
	gettimeofday(&t1, NULL);
	ll ss_us = time_diff(&t0, &t1);

	printf("uncomp_sign: %.1f us/sign\n", (double)cs_us / num_samples);
	printf("   sign_dyn: %.1f us/sign\n", (double)sd_us / num_samples);
	printf("       sign: %.1f us/sign\n", (double)ss_us / num_samples);

	// Run tests:
	int nr_cor = 0;
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++)
		nr_cor += Zf(uncompressed_verify)(pregen_h[rep], pregen_s0[rep], pregen_s1[rep], q00, q10, q11, logn, b);
	gettimeofday(&t1, NULL);
	printf("# correct = %d\n", nr_cor);
	ll cv_us = time_diff(&t0, &t1);

	nr_cor = 0;
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		nr_cor += Zf(recover_and_verify)(pregen_h[rep], s0, pregen_s[rep], q00, q10, q11, logn, b);
	}
	gettimeofday(&t1, NULL);
	printf("# correct = %d\n", nr_cor);
	ll vsr_us = time_diff(&t0, &t1);

	nr_cor = 0;
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		nr_cor += Zf(verify)(pregen_h[rep], pregen_s[rep], q00, q10, q11, logn, b);
	}
	gettimeofday(&t1, NULL);
	printf("# correct = %d\n", nr_cor);
	ll vfsr_us = time_diff(&t0, &t1);

	nr_cor = 0;
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++) {
		nr_cor += Zf(verify_nearest_plane)(pregen_h[rep], pregen_s[rep], q00, q10, q11, logn, b);
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
		nr_cor += verify_nearest_plane_tree(tree, pregen_h[rep], s0, pregen_s[rep], q00, q10, q11, logn, b);
	}
	gettimeofday(&t1, NULL);
	ll vnpt_us = time_diff(&t0, &t1);
	printf("# correct = %d\n", nr_cor);

	printf("   uncompressed_verify: %.1f us/vrfy\n", (double) cv_us / num_samples);
	printf("    recover_and_verify: %.1f us/vrfy\n", (double) vsr_us / num_samples);
	printf("                verify: %.1f us/vrfy\n", (double)vfsr_us / num_samples);
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

WorkerResult measure_signatures()
{
	uint8_t b[50 << logn];
	int8_t f[n], g[n], F[n], G[n];
	uint8_t h[n / 4];
	int16_t iq00[n], iq10[n], s1[n];
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
			Zf(keygen)(&sc, f, g, F, G, iq00, iq10, logn, b);
			Zf(complete_pubkey)(iq00, iq10, q00, q10, q11, logn);
			Zf(expand_seckey)(exp_sk, f, g, F, logn);
		}

		// make a signature of a random message...
		Zf(prng_get_bytes)(&rng, h, sizeof h);

		// Compute the signature.
		Zf(sign)(&sc, s1, exp_sk, h, logn, b);

		if (!Zf(verify)(h, s1, q00, q10, q11, logn, b))
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

void work()
{
	WorkerResult result = measure_signatures();

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

	measure_sign_speed();
	return 0;

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
