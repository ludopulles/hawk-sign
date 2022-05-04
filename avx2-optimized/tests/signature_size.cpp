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
using namespace std;

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}


/*
 * Create a dynamic Huffman table of size N which is (almost) optimal for a
 * normally distributed variable with mean 0 and standard deviation sigma.
 *
 * This method only returns an array that at index 'i' contains the encoding
 * length of 'x' when i = |x|.
 */
vector<int> generate_huffman_table(size_t n, double sigma)
{
	vector<double> freq(2 * n, -1.0);
	// Fill the frequencies in the top half.
	for (size_t u = 0; u < n; u++) {
		double u2 = (double) u * (double) u;
		freq[n + u] = exp(-u2 / (2.0 * sigma * sigma));
	}
	vector<int> p(2 * n, -1), lengths(2 * n);

	/* construct the tree */
	for (size_t u = n; --u >= 1; ) {
		size_t l = 0, r = 0; /* find 2 nodes with smallest frequencies */
		for (size_t v = 2 * n; --v > u; ) {
			if (freq[v] < 0) continue; /* v is already used */

			if (!l || freq[v] < freq[l]) r = l, l = v;
			else if (!r || freq[v] < freq[r]) r = v;
		}

		freq[u] = freq[l] + freq[r];
		freq[l] = freq[r] = -1; /* mark l and r as used */
		p[l] = p[r] = u;
	}

	lengths[1] = 0;
	for (size_t u = 2; u < 2 * n; u++) {
		lengths[u] = 1 + lengths[p[u]];
	}
	return vector<int>(lengths.begin() + n, lengths.end());
}

/*
 * Return the number of bytes needed for encoding the polynomial sig of length
 * 2^logn with a given pre-computed huffman table.
 */
size_t encoding_length(const vector<int> &table, int16_t *sig, unsigned logn)
{
	size_t res = 0; /* res is the number of bits currently required */
	for (size_t u = 0; u < MKN(logn); u++) {
		uint16_t w = (uint16_t) sig[u];
		// int abs_val = sig[u] < 0 ? (-sig[u] - 1) : sig[u];
		w ^= -(w >> 15);
		if (w >= table.size()) {
			return 0; /* value too big */
		} else {
			res += 1 + table[w]; /* sign bit and encoding of w = |sig[u]|' */
		}
	}

	// return ceil(res / 8)
	return (res + 7) / 8;
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
constexpr size_t logn = 9, n = MKN(logn);
constexpr int n_repetitions = 2500;

struct WorkerResult {
	long long sig_failed[10], sz[10], sz_sq[10];
	size_t maxsz[10], maxhuf;
	long long huf_failed, sz_h, sz_hsq;

	long long sum_s0, sumsq_s0, sum_s1, sumsq_s1;

	WorkerResult() : maxhuf(0), huf_failed(0), sz_h(0), sz_hsq(0),
		sum_s0(0), sumsq_s0(0), sum_s1(0), sumsq_s1(0)
	{
		memset(sig_failed, 0, sizeof sig_failed);
		memset(sz, 0, sizeof sz);
		memset(sz_sq, 0, sizeof sz_sq);
		memset(maxsz, 0, sizeof sz_sq);
	}

	void add(size_t s, int lb) {
		if (s == 0) {
			sig_failed[lb]++;
		} else {
			sz[lb] += s;
			sz_sq[lb] += s*s;
			maxsz[lb] = max(maxsz[lb], s);
		}
	}

	void combine(const WorkerResult &res) {
		huf_failed += res.huf_failed;
		sz_h += res.sz_h;
		sz_hsq += res.sz_hsq;
		maxhuf = max(maxhuf, res.maxhuf);

		for (size_t i = 0; i < 10; i++) {
			sig_failed[i] += res.sig_failed[i];
			sz[i] += res.sz[i];
			sz_sq[i] += res.sz_sq[i];
			maxsz[i] = max(maxsz[i], res.maxsz[i]);
		}

		sum_s0 += res.sum_s0;
		sumsq_s0 += res.sumsq_s0;
		sum_s1 += res.sum_s1;
		sumsq_s1 += res.sumsq_s1;
	}
};

vector<int> huffman_s0, huffman_s1;

WorkerResult measure_signatures() {
	uint8_t b[50 << logn], h[n / 4];
	int8_t f[n], g[n], F[n], G[n];
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
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b);

		// make a signature of a random message...
		inner_shake256_extract(&sc, h, sizeof h);

		// Compute the signature.
		Zf(complete_sign)(&sc, s0, s1, f, g, F, G, h, logn, b);

		for (size_t u = 0; u < n; u++) {
			result.sum_s0 += s0[u];
			result.sumsq_s0 += s0[u] * s0[u];
			result.sum_s1 += s1[u];
			result.sumsq_s1 += s1[u] * s1[u];
		}

		for (int lobits = 3; lobits < 10; lobits++) {
			size_t sz = Zf(encode_sig)(NULL, 0, s1, logn, lobits);
			result.add(sz, lobits);
		}
		// size_t szh = Zf(encode_sig_huffman)(NULL, 0, s1, logn);
		size_t szh = encoding_length(huffman_s1, s1, logn);
		if (szh != 0) {
			result.sz_h += szh;
			result.sz_hsq += szh*szh;
			result.maxhuf = max(result.maxhuf, szh);
		} else result.huf_failed++;
	}
	return result;
}

WorkerResult tot;
mutex mx;

constexpr fpr sigma_sig = { v: 1.292 };
constexpr fpr verif_margin = { v: 1.1 };

void work() {
	WorkerResult result = measure_signatures();

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

	huffman_s0 = generate_huffman_table(512, 405.304);
	huffman_s1 = generate_huffman_table(512,  62.095);

	const int nthreads = 4;
	std::thread* pool[nthreads-1];
	for (int i = 0; i < nthreads-1; i++) pool[i] = new std::thread(work);
	work();
	for (int i = 0; i < nthreads-1; i++) pool[i]->join(), delete pool[i];

	long long runs = n_repetitions * nthreads;

	printf("\nNumber of signatures signed = %lld\n\n", runs);
	double avg_s0 = double(tot.sum_s0) / (runs * n);
	double avg_s1 = double(tot.sum_s1) / (runs * n);
	double var_s0 = (long double)(tot.sumsq_s0) / (runs * n) - avg_s0 * avg_s0;
	double var_s1 = (long double)(tot.sumsq_s1) / (runs * n) - avg_s1 * avg_s1;
	printf("Average coefficient of s0 ~ %7.3f +/- %7.3Lf\n", avg_s0, sqrtl(var_s0));
	printf("Average coefficient of s1 ~ %7.3f +/- %7.3Lf\n", avg_s1, sqrtl(var_s1));

	printf("lo_bits | max sig | avg sig. size (B) | #fails\n");
	for (int lb = 3; lb < 10; lb++) {
		double avg_sz = (double)tot.sz[lb] / (runs - tot.sig_failed[lb]);
		double std_sz = (double)tot.sz_sq[lb] / (runs - tot.sig_failed[lb]) - avg_sz * avg_sz;
		printf("%7u | %7zu | %.1f (std %5.1f) | %5lld\n",
			lb, tot.maxsz[lb], avg_sz, std_sz, tot.sig_failed[lb]);
	}
	
	double avg_h = (double) tot.sz_h / (runs - tot.huf_failed);
	double std_h = (double) tot.sz_hsq / (runs - tot.huf_failed) - avg_h * avg_h;
	printf("Huffman | %7zu | %.1f (std %5.1f) | %5lld\n", tot.maxhuf, avg_h, std_h, tot.huf_failed);

	return 0;
}
