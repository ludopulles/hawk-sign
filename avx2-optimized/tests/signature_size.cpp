#include <cassert>
#include <climits>
#include <cstdio>
#include <vector>

#include <mutex>
#include <thread>

extern "C" {
	#define restrict
	#include "../inner.h"
}
using namespace std;

constexpr int MAXN = 1024;

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
		if (w >= table.size()) return 0; /* value too big */
		res += 1 + table[w]; /* sign bit and encoding of w = |sig[u]|' */
	}
	return (res + 7u) / 8;
}


/*
 * Returns number of bytes for signature.
 */
size_t cost_golomb(const int16_t *x, unsigned logn, size_t lo_bits)
{
	size_t n = MKN(logn), u, nbits = 0;
	for (u = 0; u < n; u ++) {
		uint16_t w = (uint16_t) x[u];
		w ^= -(w >> 15);
		w >>= lo_bits;
		nbits += (lo_bits + 1) + (w + 1);
	}
	return (nbits + 7u) / 8;
}

int lg2(int16_t x) { return x == 1 ? 0 : 1 + lg2(x >> 1); }
int expGolombCost(size_t absx) {
	// 0: 0, +-1: 10x, +-2: 1100x, +-3: 1101x, etc.
	if (absx == 0) return 1;
	return 2 * lg2(absx) + 3; // 1 -> 3; 2,3 -> 5; 4,5,6,7 -> 7; etc.
	return 2 + absx; // Alternative: print unary & sign
}

/*
 * Returns number of bytes for signature.
 */
size_t cost_expgolomb(const int16_t *x, unsigned logn, size_t lo_bits)
{
	size_t n = MKN(logn), u, nbits = 0;
	for (u = 0; u < n; u ++) {
		int16_t w;
		w = ((x[u] >> (lo_bits - 1)) + 1) >> 1;
		w = w < 0 ? -w : w;
		nbits += lo_bits + expGolombCost(w);
	}
	return (nbits + 7u) / 8;
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
constexpr int n_repetitions = 100;

// Taken from hawk.h:
/*
 * The number of bytes that are required for the salt that gets prepended to a
 * message before hashing.
 * 
 * From GPV08, having \lambda security bits and a transcript size of Q_s, having k = \lambda + \log_2(Q_s) salt bits gives a hash collision of message to happen with probability at most 2^-k.
 * So for n = 512, take k = 192 and for n = 1024, take k = 320.
 */
#define HAWK_SALT_SIZE(logn) ((logn) <= 9 ? 24u : 40u)

struct WorkerResult {
	// s1 <- encode_sig()
	long long sig_failed[10], sz[10], sz_sq[10], maxsz[10];
	// (s0, s1) <- encode_uncomp_sig()
	long long csig_failed[11], csz[11], csz_sq[11], maxcsz[11];
	// encode_huffman
	long long huf_failed, sz_h, sz_hsq, maxhuf;
	long long chuf_failed, csz_h, csz_hsq, maxchuf;
	long long sum_s0, sumsq_s0, sum_s1, sumsq_s1;

	WorkerResult() :
		huf_failed(0), sz_h(0), sz_hsq(0), maxhuf(0),
		chuf_failed(0), csz_h(0), csz_hsq(0), maxchuf(0),
		sum_s0(0), sumsq_s0(0), sum_s1(0), sumsq_s1(0)
	{
		memset(sig_failed, 0, sizeof sig_failed);
		memset(sz, 0, sizeof sz);
		memset(sz_sq, 0, sizeof sz_sq);
		memset(maxsz, 0, sizeof sz_sq);

		memset(csig_failed, 0, sizeof csig_failed);
		memset(csz, 0, sizeof csz);
		memset(csz_sq, 0, sizeof csz_sq);
		memset(maxcsz, 0, sizeof csz_sq);
	}

	void add(size_t s, int lb) {
		if (s == 0) {
			sig_failed[lb]++;
		} else {
			sz[lb] += s;
			sz_sq[lb] += s*s;
			maxsz[lb] = max(maxsz[lb], (long long) s);
		}
	}

	void cadd(size_t s, int lb) {
		if (s == 0) {
			csig_failed[lb]++;
		} else {
			csz[lb] += s;
			csz_sq[lb] += s*s;
			maxcsz[lb] = max(maxcsz[lb], (long long) s);
		}
	}

	void combine(const WorkerResult &res) {
		huf_failed += res.huf_failed;
		sz_h += res.sz_h;
		sz_hsq += res.sz_hsq;
		maxhuf = max(maxhuf, res.maxhuf);

		chuf_failed += res.chuf_failed;
		csz_h += res.csz_h;
		csz_hsq += res.csz_hsq;
		maxchuf = max(maxchuf, res.maxchuf);

		for (size_t i = 0; i < 10; i++) {
			sig_failed[i] += res.sig_failed[i];
			sz[i] += res.sz[i];
			sz_sq[i] += res.sz_sq[i];
			maxsz[i] = max(maxsz[i], res.maxsz[i]);
		}

		for (size_t i = 0; i < 11; i++) {
			csig_failed[i] += res.csig_failed[i];
			csz[i] += res.csz[i];
			csz_sq[i] += res.csz_sq[i];
			maxcsz[i] = max(maxcsz[i], res.maxcsz[i]);
		}

		sum_s0 += res.sum_s0;
		sumsq_s0 += res.sumsq_s0;
		sum_s1 += res.sum_s1;
		sumsq_s1 += res.sumsq_s1;
	}
};

vector<int> huffman_512_s0, huffman_512_s1;
vector<int> huffman_1024_s0, huffman_1024_s1;

WorkerResult measure_signatures(unsigned logn) {
	uint8_t b[48 * MAXN], h[MAXN / 4];
	int8_t f[MAXN], g[MAXN], F[MAXN], G[MAXN];
	int16_t s0[MAXN], s1[MAXN];
	fpr q00[MAXN], q10[MAXN], q11[MAXN];
	size_t n;
	unsigned char seed[48];
	inner_shake256_context sc;

	n = MKN(logn);

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
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
		while (!Zf(uncompressed_sign)(&sc, s0, s1, f, g, F, G, h, logn, b)) {}

		for (size_t u = 0; u < n; u++) {
			result.sum_s0 += s0[u];
			result.sumsq_s0 += s0[u] * s0[u];
			result.sum_s1 += s1[u];
			result.sumsq_s1 += s1[u] * s1[u];
		}

		for (int lobits = 3; lobits < 10; lobits++) {
			// size_t sz = Zf(encode_sig)(NULL, 0, s1, logn, lobits);
			size_t sz = cost_golomb(s1, logn, lobits);
			result.add(sz, lobits);
		}

		// size_t szh = Zf(encode_sig_huffman)(NULL, 0, s1, logn);
		size_t szh = encoding_length(logn == 10 ? huffman_1024_s1 : huffman_512_s1, s1, logn);
		if (szh != 0) {
			result.sz_h += szh;
			result.sz_hsq += szh*szh;
			result.maxhuf = max(result.maxhuf, (long long) szh);
		} else {
			result.huf_failed++;
		}

		for (int lobits = 5; lobits < 11; lobits++) {
			size_t sz = Zf(encode_uncomp_sig)(NULL, 0, s0, s1, logn, lobits, logn - 4);
			result.cadd(sz, lobits);
		}
		size_t cszh = encoding_length(logn == 10 ? huffman_1024_s0 : huffman_512_s0, s0, logn)
			+ encoding_length(logn == 10 ? huffman_1024_s1 : huffman_512_s1, s1, logn);
		if (cszh != 0) {
			result.csz_h += cszh;
			result.csz_hsq += cszh*cszh;
			result.maxchuf = max(result.maxchuf, (long long) cszh);
		} else {
			result.chuf_failed++;
		}
	}
	return result;
}

WorkerResult tot;
mutex mx;

void work(unsigned logn) {
	WorkerResult result = measure_signatures(logn);

	{
		lock_guard<mutex> guard(mx);
		tot.combine(result);
	}
}

int main() {
	huffman_512_s0 = generate_huffman_table(2048, 355);
	huffman_512_s1 = generate_huffman_table(512,  61);

	huffman_1024_s0 = generate_huffman_table(8192, 959);
	huffman_1024_s1 = generate_huffman_table(1024,  118);

	const int nthreads = 4;
	std::thread* pool[nthreads-1];

	for (unsigned logn = 9; logn <= 10; logn++) {
		size_t n = MKN(logn);
		tot = WorkerResult();

		for (int i = 0; i < nthreads-1; i++) pool[i] = new std::thread(work, logn);
		work(logn);
		for (int i = 0; i < nthreads-1; i++) pool[i]->join(), delete pool[i];

		long long runs = n_repetitions * nthreads;

		printf("\nlogn = %u, num_samples = %lld\n", logn, runs);
		double avg_s0 = double(tot.sum_s0) / (runs * n);
		double avg_s1 = double(tot.sum_s1) / (runs * n);
		double var_s0 = double(tot.sumsq_s0) / (runs * n) - avg_s0 * avg_s0;
		double var_s1 = double(tot.sumsq_s1) / (runs * n) - avg_s1 * avg_s1;
		printf("Average coefficient of s0 ~ %7.3f +/- %7.3f\n", avg_s0, sqrt(var_s0));
		printf("Average coefficient of s1 ~ %7.3f +/- %7.3f\n", avg_s1, sqrt(var_s1));

		long long salt_and_header = 1 + HAWK_SALT_SIZE(logn);

		printf("lo_bits | max sig | avg sig. size (B) | #fails\n");
		for (int lb = 3; lb < 10; lb++) {
			double avg_sz = (double)tot.sz[lb] / (runs - tot.sig_failed[lb]);
			double var_sz = (double)tot.sz_sq[lb] / (runs - tot.sig_failed[lb]) - avg_sz * avg_sz;
			printf("%7u | %7lld | %.1f (std %5.1f) | %5lld\n",
				lb,
				salt_and_header + tot.maxsz[lb],
				salt_and_header + avg_sz,
				sqrt(var_sz),
				tot.sig_failed[lb]
			);
		}

		double avg_h = (double) tot.sz_h / (runs - tot.huf_failed);
		double var_h = (double) tot.sz_hsq / (runs - tot.huf_failed) - avg_h * avg_h;
		printf("Huffman | %7lld | %.1f (std %5.1f) | %5lld\n",
			salt_and_header + tot.maxhuf,
			salt_and_header + avg_h,
			sqrt(var_h),
			tot.huf_failed
		);

		printf("Complete signature (s0, s1) encodings:\n");
		for (int lb = 7; lb < 11; lb++) {
			double avg_sz = (double)tot.csz[lb] / (runs - tot.csig_failed[lb]);
			double var_sz = (double)tot.csz_sq[lb] / (runs - tot.csig_failed[lb]) - avg_sz * avg_sz;
			printf("(%2u, %d) | %7lld | %.1f (std %5.1f) | %5lld\n",
				lb,
				logn - 4,
				salt_and_header + tot.maxcsz[lb],
				salt_and_header + avg_sz,
				sqrt(var_sz),
				tot.csig_failed[lb]
			);
		}

		avg_h = (double) tot.csz_h / (runs - tot.chuf_failed);
		var_h = (double) tot.csz_hsq / (runs - tot.chuf_failed) - avg_h * avg_h;
		printf("Huffman | %7lld | %.1f (std %5.1f) | %5lld\n",
			salt_and_header + tot.maxchuf,
			salt_and_header + avg_h,
			sqrt(var_h),
			tot.chuf_failed
		);

	}

	return 0;
}
