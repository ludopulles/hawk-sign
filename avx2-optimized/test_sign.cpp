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

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
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

struct WorkerResult {
	ll iters, babai_fail, sig_fail, sigsize, sqsigsize, maxsigsz;

	WorkerResult() : iters(0), babai_fail(0),
		sig_fail(0), sigsize(0), sqsigsize(0), maxsigsz(0) {}

	void combine(const WorkerResult &res) {
		iters += res.iters;
		babai_fail += res.babai_fail;
		sig_fail += res.sig_fail;
		sigsize += res.sigsize;
		sqsigsize += res.sqsigsize;
		maxsigsz = std::max(maxsigsz, res.maxsigsz);
	}
};

WorkerResult measure_signatures(fpr isigma_kg, fpr isigma_sig, uint32_t bound) {
	uint8_t b[42 << logn];
	int8_t f[n], g[n], F[n], G[n], h[n];
	int16_t s0[n], recs0[n], s1[n];
	fpr q00[n], q10[n], q11[n];
	unsigned char seed[48];
	inner_shake256_context sc;
	const int n_repetitions = 1024 * 16;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	WorkerResult result;
	result.iters = n_repetitions;

	for (int rep = 0; rep < n_repetitions; rep++) {
		if ((rep & 1023) == 0) {
			// Generate key pair.
			Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, isigma_kg, logn, b);
		}

		// make a signature of a random message...
		random_hash(h, logn);

		// Compute the signature.
		// Zf(complete_sign)(&sc, s0, s1, f, g, F, G, h, isigma_sig, bound, logn, b);
		Zf(sign)(&sc, s1, f, g, h, isigma_sig, bound, logn, b);

		if (!Zf(verify)(h, recs0, s1, q00, q10, q11, bound, logn, b))
			result.babai_fail++;
		size_t sig_sz = Zf(encode_sig)(NULL, 0, s1, logn, 5);
		if (sig_sz == 0) {
			result.sig_fail++;
		} else {
			result.sigsize += sig_sz;
			result.sqsigsize += sig_sz * sig_sz;
			result.maxsigsz = std::max(result.maxsigsz, (ll) sig_sz);
		}
	}
	return result;
}

WorkerResult tot;
std::mutex mx;

constexpr fpr sigma_kg  = { v: 1.425 };
constexpr fpr sigma_sig = { v: 1.292 };
constexpr fpr verif_margin = { v: 1.1 };

uint32_t getbound() {
	return fpr_floor(fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_of(n))));
}

void work() {
	WorkerResult result = measure_signatures(fpr_inv(sigma_kg), fpr_inv(sigma_sig), getbound());

	{
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

	const int nthreads = 4;
	std::thread* pool[nthreads-1];
	for (int i = 0; i < nthreads-1; i++) pool[i] = new std::thread(work);
	work();
	for (int i = 0; i < nthreads-1; i++) pool[i]->join(), delete pool[i];

	printf("\n");
	printf("# Signatures signed:      %lld\n", tot.iters);
	printf("# Babai roundings failed: %lld\n", tot.babai_fail);
	printf("# Sig. coding failed:     %lld\n", tot.sig_fail);

	ll nsigs = tot.iters - tot.sig_fail;
	double avg_sz = (double)tot.sigsize / nsigs;
	double std_sz = (double)tot.sqsigsize / nsigs - avg_sz * avg_sz;
	printf("# Average sig. size = %.2f (± %.2f) and <= %lld\n", avg_sz, std_sz, tot.maxsigsz);

	return 0;
}
