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
using namespace std;

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

struct WorkerResult {
	long long num_signed[10], sig_failed[10], sz[10], sz_sq[10];
	size_t maxsz[10], maxhuf;
	long long num_huf, num_huf_fail, sz_h, sz_hsq;

	WorkerResult() : maxhuf(0), num_huf(0), num_huf_fail(0),
		sz_h(0), sz_hsq(0)
	{
		memset(num_signed, 0, sizeof num_signed);
		memset(sig_failed, 0, sizeof sig_failed);
		memset(sz, 0, sizeof sz);
		memset(sz_sq, 0, sizeof sz_sq);
		memset(maxsz, 0, sizeof sz_sq);
	}

	void add(size_t s, int lb) {
		if (s == 0) {
			sig_failed[lb]++;
		} else {
			num_signed[lb]++;
			sz[lb] += s;
			sz_sq[lb] += s*s;
			maxsz[lb] = max(maxsz[lb], s);
		}
	}

	void combine(const WorkerResult &res) {
		num_huf += res.num_huf;
		num_huf_fail += res.num_huf_fail;
		sz_h += res.sz_h;
		sz_hsq += res.sz_hsq;
		maxhuf = max(maxhuf, res.maxhuf);
		for (unsigned i=0; i < 10; i++) {
			num_signed[i] += res.num_signed[i];
			sig_failed[i] += res.sig_failed[i];
			sz[i] += res.sz[i];
			sz_sq[i] += res.sz_sq[i];
			maxsz[i] = max(maxsz[i], res.maxsz[i]);
		}
	}
};

const int n_repetitions = 1000;

WorkerResult measure_signatures(fpr isigma_kg, fpr isigma_sig, uint32_t bound) {
	uint8_t b[42 << logn];
	int8_t f[n], g[n], F[n], G[n], h[n];
	int16_t s0[n], s1[n];
	fpr q00[n], q10[n], q11[n];
	unsigned char seed[48];
	inner_shake256_context sc;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	WorkerResult result;
	for (int rep = 0; rep < n_repetitions; rep++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, isigma_kg, logn, b);

		// make a signature of a random message...
		random_hash(h, logn);

		// Compute the signature.
		Zf(complete_sign)(&sc, s0, s1, f, g, F, G, h, isigma_sig, bound, logn, b);

		for (int lobits = 3; lobits < 10; lobits++) {
			size_t sz = Zf(encode_sig)(NULL, 0, s1, logn, lobits);
			result.add(sz, lobits);
		}
		size_t szh = Zf(encode_sig_huffman)(NULL, 0, s1, logn);
		if (szh != 0) {
			result.num_huf++;
			result.sz_h += szh;
			result.sz_hsq += szh*szh;
			result.maxhuf = max(result.maxhuf, szh);
		} else result.num_huf_fail++;
	}
	return result;
}

WorkerResult tot;
mutex mx;

constexpr fpr sigma_kg  = { v: 1.425 };
constexpr fpr sigma_sig = { v: 1.292 };
constexpr fpr verif_margin = { v: 1.1 };

uint32_t getbound() {
	return fpr_floor(fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_of(n))));
}

void work() {
	WorkerResult result = measure_signatures(fpr_inv(sigma_kg), fpr_inv(sigma_sig), getbound());

	{
		lock_guard<mutex> guard(mx);
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

	printf("\nRuns: %d\n", n_repetitions * nthreads);
	printf("lo_bits | max sig | avg sig. size (B) | #fails\n");
	for (int lb = 3; lb < 10; lb++) {
		double avg_sz = (double)tot.sz[lb] / tot.num_signed[lb];
		double std_sz = (double)tot.sz_sq[lb] / tot.num_signed[lb] - avg_sz * avg_sz;
		printf("%7u | %7zu | %.1f (std %5.1f) | %5lld\n",
			lb, tot.maxsz[lb], avg_sz, std_sz, tot.sig_failed[lb]);
	}
	
	double avg_h = (double) tot.sz_h / tot.num_huf;
	double std_h = (double) tot.sz_hsq / tot.num_huf - avg_h * avg_h;
	printf("Huffman | %7zu | %.1f (std %5.1f) | %5lld\n", tot.maxhuf, avg_h, std_h, tot.num_huf_fail);

	return 0;
}
