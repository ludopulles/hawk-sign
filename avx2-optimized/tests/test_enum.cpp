#include <assert.h>
#include <stdio.h>
// x86_64 specific:
#include <sys/time.h>

extern "C" {
	#define restrict
	#include "../inner.h"
}

// concurrency:
#include <thread>
#include <mutex>

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

unsigned sqnorm(int8_t *p, size_t logn) {
	unsigned res = 0;
	for (size_t u = 0, n = MKN(logn); u < n; u ++)
		res += (unsigned)p[u] * p[u];
	return res;
}

unsigned sqnorm16(int16_t *p, size_t logn) {
	unsigned res = 0;
	for (size_t u = 0, n = MKN(logn); u < n; u ++)
		res += (unsigned)p[u] * p[u];
	return res;
}

void print_int16(int16_t *p, size_t logn) {
	for (size_t u = 0, n = MKN(logn); u < n; u ++) {
		if (u) printf(" ");
		printf("%d", p[u]);
	}
	printf("\n");
}


// Huffman Encoding:
#define MAX_Q10 (4096)
#define MAXLEN_Q10 (56)
struct {
	uint16_t a[2*MAX_Q10][2]; // { left child, right child }
	uint16_t p[2*MAX_Q10]; // parent
} huffman_tree;

void create_huffman_tree() {
	float freq[2*MAX_Q10], sigma = 512;
	int len[2*MAX_Q10];
	for (int x = 0; x < MAX_Q10; x++)
		freq[MAX_Q10 + x] = exp((float)-x * x / (2.0 * sigma * sigma)),
		len[MAX_Q10 + x] = 0;

	// construct the tree
	for (uint16_t node = MAX_Q10; --node >= 1; ) {
		// find 2 nodes with smallest frequencies
		uint16_t l = 0, r = 0;
		for (uint16_t idx = node; ++idx < 2*MAX_Q10; ) {
			if (freq[idx] < 0) continue;
			if (!l || freq[idx] < freq[l])
				r = l, l = idx;
			else if (!r || freq[idx] < freq[r])
				r = idx;
		}
		// hide frequency
		freq[node] = freq[l] + freq[r];
		freq[l] = freq[r] = -1;
		len[node] = 1 + (len[l] > len[r] ? len[l] : len[r]);
		huffman_tree.p[l] = huffman_tree.p[r] = node;
		huffman_tree.a[node][0] = l;
		huffman_tree.a[node][1] = r;
	}
	assert(len[1] <= MAXLEN_Q10 && "Longest path is longer than expected!");
}

size_t Zf(huffman_encode)(void *out, size_t max_out_len, const int16_t *x, unsigned logn) {
	uint8_t *buf = (uint8_t *)out;
	size_t n = MKN(logn), u, v = 0;
	uint8_t acc = 0, acc_len = 0, steps[MAXLEN_Q10];

	for (u = 0; u < n; u ++)
		if (x[u] <= -MAX_Q10 || x[u] >= MAX_Q10) return 0;
	for (u = 0; u < n; u ++) {
		// Get sign and absolute value of next integer; push the sign bit.
		acc <<= 1;
		int16_t t = x[u];
		if (t < 0) t = -t, acc |= 1;

		if (++acc_len == 8) {
			if (buf != NULL) {
				if (v >= max_out_len) return 0;
				buf[v] = acc;
			}
			acc_len = acc = 0; // reset acc
			v++;
		}

		size_t nsteps = 0;
		for (int16_t idx = MAX_Q10 + t; idx > 1; ) {
			int16_t next_idx = huffman_tree.p[idx];
			steps[nsteps++] = huffman_tree.a[next_idx][1] == idx ? 1 : 0;
			idx = next_idx;
		}

		/* printf("%d: ", x[u]);
		for (size_t w = nsteps; w --> 0; )
			printf("%d", steps[nsteps]);
		printf("\n"); */

		while (nsteps --> 0) {
			acc = (acc << 1) | steps[nsteps];
			if (++acc_len == 8) {
				if (buf != NULL) {
					if (v >= max_out_len) return 0;
					buf[v] = acc;
				}
				acc_len = acc = 0; // reset acc
				v++;
			}
		}
	}

	// Flush remaining bits (if any).
	if (acc_len > 0) {
		if (buf != NULL) {
			if (v >= max_out_len) return 0;
			buf[v] = (uint8_t)(acc << (8 - acc_len));
		}
		v++;
	}
	return v;
}


/* see codec.c, Zf(comp_encode) */
size_t
Zf(encode_q00)(
	void *out, size_t max_out_len,
	const int16_t *x, unsigned logn)
{
	uint8_t *buf;
	size_t n, u, v;
	uint64_t acc;
	unsigned acc_len;

	n = (size_t)1 << logn;
	buf = (uint8_t *)out;

	/*
	 * Make sure that all values are within the -2047..+2047 range.
	 * Moreover, q00 = q00^* so q00 + q00^* must be in fact an integer.
	 */
	for (u = 1; u < n/2; u ++) {
		if (x[u] < -2047 || x[u] > +2047 || x[u] + x[n-u] != 0) {
			return 0;
		}
	}

	/*
	 * The first value of q00 is of the order sigma_kg^2 2d,
	 * but definitely << 8d when sigma_kg = 1.425, by the standard
	 * Laurent-Massart bound (see https://en.wikipedia.org/wiki/Chi-squared_distribution#Concentration).
	 */
	acc_len = 3 + logn;
	acc = (uint32_t) x[0];
	v = 0;

	int32_t max_val = 1 << acc_len;
	if (x[0] <= -max_val || x[0] >= max_val) {
		return 0;
	}

	/*
	 * Produce all full bytes.
	 */
	while (acc_len >= 8) {
		acc_len -= 8;
		if (buf != NULL) {
			if (v >= max_out_len) {
				return 0;
			}
			buf[v] = (uint8_t)(acc >> acc_len);
		}
		v ++;
	}

	for (u = 1; u < n/2; u ++) {
		int t;
		unsigned w;

		/*
		 * Get sign and absolute value of next integer; push the
		 * sign bit.
		 */
		acc <<= 1;
		t = x[u];
		if (t < 0) {
			t = -t;
			acc |= 1;
		}
		w = (unsigned)t;

		/*
		 * Push the low `lo_bits` bits of the absolute value.
		 */

		const int lo_bits = 5;
		acc <<= lo_bits;
		acc |= w & ((1U << lo_bits) - 1);
		w >>= lo_bits;

		/*
		 * We pushed exactly `lo_bits + 1` bits.
		 */
		acc_len += (lo_bits + 1);

		/*
		 * Push as many zeros as necessary, then a one. Since the
		 * absolute value is at most 2047, w can only range up to
		 * 15 at this point, thus we will add at most 16 bits
		 * here. With the 8 bits above and possibly up to 7 bits
		 * from previous iterations, we may go up to 31 bits, which
		 * will fit in the accumulator, which is an uint32_t.
		 */
		acc <<= (w + 1);
		acc |= 1;
		acc_len += w + 1;

		/*
		 * Produce all full bytes.
		 */
		while (acc_len >= 8) {
			acc_len -= 8;
			if (buf != NULL) {
				if (v >= max_out_len) {
					return 0;
				}
				buf[v] = (uint8_t)(acc >> acc_len);
			}
			v ++;
		}
	}

	/*
	 * Flush remaining bits (if any).
	 */
	if (acc_len > 0) {
		if (buf != NULL) {
			if (v >= max_out_len) {
				return 0;
			}
			buf[v] = (uint8_t)(acc << (8 - acc_len));
		}
		v ++;
	}

	return v;
}

/* see codec.c, Zf(comp_encode) */
size_t
Zf(encode_q10)(void *out, size_t max_out_len, const int16_t *x, unsigned logn, const int lim, const int lo_bits) {
	uint8_t *buf = (uint8_t *)out;
	size_t n = MKN(logn), u, v = 0;
	uint64_t acc = 0;
	unsigned acc_len = 0;

	for (u = 0; u < n; u ++)
		if (x[u] < -lim || x[u] > lim) return 0;
	for (u = 0; u < n; u ++) {
		int t;
		unsigned w;

		// Get sign and absolute value of next integer; push the sign bit.
		acc <<= 1;
		t = x[u];
		if (t < 0) t = -t, acc |= 1;
		w = (unsigned)t;

		// Push the low `lo_bits` bits of the absolute value.
		// const int lo_bits = 8;

		acc <<= lo_bits;
		acc |= w & ((1U << lo_bits) - 1);
		w >>= lo_bits;

		// We pushed exactly `lo_bits + 1` bits.
		acc_len += (lo_bits + 1);

		// TODO: perhaps this still works, but perhaps we need uint64_t...
		/* Push as many zeros as necessary, then a one. Since the
		 * absolute value is at most 4095, w can only range up to
		 * 7 at this point, thus we will add at most 8 bits
		 * here. With the 10 bits above and possibly up to 7 bits
		 * from previous iterations, we may go up to 25 bits, which
		 * will fit in the accumulator, which is an uint32_t. */
		acc <<= (w + 1);
		acc |= 1;
		acc_len += w + 1;

		// Produce all full bytes.
		while (acc_len >= 8) {
			acc_len -= 8;
			if (buf != NULL) {
				if (v >= max_out_len) return 0;
				buf[v] = (uint8_t)(acc >> acc_len);
			}
			v ++;
		}
	}

	// Flush remaining bits (if any).
	if (acc_len > 0) {
		if (buf != NULL) {
			if (v >= max_out_len) return 0;
			buf[v] = (uint8_t)(acc << (8 - acc_len));
		}
		v ++;
	}
	return v;
}



size_t bitreverse(size_t num, size_t logn) {
	return logn == 1 ? num : (bitreverse(num / 2, logn - 1) | ((num & 1) << (logn - 1)));
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
const size_t logn = 9, n = MKN(logn);

struct WorkerResult
{
	long long num_iters, kgEnc, kgHuf, enumEnc, enumHuf;
	double time_sum;

	WorkerResult() : num_iters(0), kgEnc(0), kgHuf(0), enumEnc(0), enumHuf(0), time_sum(0.0) {}

	void combine(const WorkerResult &res) {
		num_iters += res.num_iters;
		kgEnc += res.kgEnc;
		kgHuf += res.kgHuf;
		enumEnc += res.enumEnc;
		enumHuf += res.enumHuf;
		time_sum += res.time_sum;
	}
};

void
do_enumeration(const int8_t *f, const int8_t *g, int8_t *F, int8_t *G, const fpr *q00, fpr *q10, const int16_t *iq00, int16_t *iq10)
{
	constexpr int NP = 10, STEPS = 10*1000;
	int16_t newq10[n];
	int pos[NP], add[NP];
	for (int i = (NP-1)*STEPS; i --> 0; ) {
		// Try to modify q10 at NP points
		memcpy(newq10, iq10, sizeof newq10);
		int np = 2 + i/STEPS;
		for (int it = np; it --> 0; ) {
			pos[it] = rand() % n;
			add[it] = rand() % 2;

			for (size_t j = 0; j < n - pos[it]; j++)
				newq10[pos[it]+j  ] += add[it] ? iq00[j] : -iq00[j];
			for (size_t j = n - pos[it]; j < n; j++)
				newq10[pos[it]+j-n] += add[it] ? -iq00[j] : iq00[j];
		}

		if (sqnorm16(newq10, logn) < sqnorm16(iq10, logn)) {
			// printf("Improvement (%d): %d ==> %d\n", np, sqnorm16(iq10, logn), sqnorm16(newq10, logn));
			for (size_t it = np; it --> 0; ) {
				for (size_t j = 0; j < n - pos[it]; j++) {
					F[pos[it]+j  ] += add[it] ? f[j] : -f[j];
					G[pos[it]+j  ] += add[it] ? g[j] : -g[j];
					q10[pos[it]+j  ] = fpr_add(q10[pos[it]+j  ], add[it] ? q00[j] : fpr_neg(q00[j]));
					iq10[pos[it]+j  ] += add[it] ? iq00[j] : -iq00[j];
				}
				for (size_t j = n - pos[it]; j < n; j++) {
					F[pos[it]+j-n] += add[it] ? -f[j] : f[j];
					G[pos[it]+j-n] += add[it] ? -g[j] : g[j];
					q10[pos[it]+j-n] = fpr_add(q10[pos[it]+j-n], add[it] ? fpr_neg(q00[j]) : q00[j]);
					iq10[pos[it]+j-n] += add[it] ? -iq00[j] : iq00[j];
				}
			}
		}
	}

	// bruteforce optimising last position
	int lp = n-1, p = 0;
	while (p != lp) {
		for (add[0] = 0; add[0] < 2; add[0]++) {
improveq10:
			memcpy(newq10, iq10, sizeof newq10);
			for (size_t j = 0; j < n - p; j++)
				newq10[p+j  ] += add[0] ? iq00[j] : -iq00[j];
			for (size_t j = n - p; j < n; j++)
				newq10[p+j-n] += add[0] ? -iq00[j] : iq00[j];

			if (sqnorm16(newq10, logn) < sqnorm16(iq10, logn)) {
				printf("Improvement (1, %d): %d ==> %d\n", p, sqnorm16(iq10, logn), sqnorm16(newq10, logn));
				lp = p;

				for (size_t j = 0; j < n - p; j++) {
					F[p+j  ] += add[0] ? f[j] : -f[j];
					G[p+j  ] += add[0] ? g[j] : -g[j];
					q10[p+j  ] = fpr_add(q10[p+j  ], add[0] ? q00[j] : fpr_neg(q00[j]));
					iq10[p+j  ] += add[0] ? iq00[j] : -iq00[j];
				}
				for (size_t j = n - p; j < n; j++) {
					F[p+j-n] += add[0] ? -f[j] : f[j];
					G[p+j-n] += add[0] ? -g[j] : g[j];
					q10[p+j-n] = fpr_add(q10[p+j-n], add[0] ? fpr_neg(q00[j]) : q00[j]);
					iq10[p+j-n] += add[0] ? -iq00[j] : iq00[j];
				}
				// repeat this position
				goto improveq10;
			}
		}
		if (++p == n) p = 0;
	}
}

WorkerResult measure_keygen()
{
	union { // use union to ensure alignment
		uint8_t b[28*512]; // 14 kB temporary memory, 17.5 kB total
		uint64_t dummy_i64;
		fpr dummy_fpr;
	} tmp;
	int8_t f[n], g[n], F[n], G[n];
	int16_t iq00[n], iq10[n];
	fpr q00[n], q10[n];
	unsigned char seed[48];
	inner_shake256_context sc;

	struct timeval t0, t1;
	const int n_repetitions = 1;

	WorkerResult result;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);
	for (int _ = 0; _ < n_repetitions; _++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, iq00, iq10, logn, tmp.b);
		result.num_iters++;

		// Check size before enumeration (with Babai though)
		size_t eq10 = Zf(encode_q10)(NULL, 0, iq10, logn, 2047, 8);
		size_t hq10 = Zf(huffman_encode)(NULL, 0, iq10, logn);
		printf("Size before: %zu, %zu\n", eq10, hq10);
		result.kgEnc += eq10;
		result.kgHuf += hq10;

		// Improve F, G with enumeration
		Zf(int16_to_fft)(q00, iq00, logn);
		Zf(int16_to_fft)(q10, iq10, logn);
		do_enumeration(f, g, F, G, q00, q10, iq00, iq10);

		// Check size after enumeration
		size_t epq10 = Zf(encode_q10)(NULL, 0, iq10, logn, 2047, 8);
		size_t hpq10 = Zf(huffman_encode)(NULL, 0, iq10, logn);
		printf("Size after:  %zu, %zu\n", epq10, hpq10);
		result.enumEnc += epq10;
		result.enumHuf += hpq10;
	}

	gettimeofday(&t1, NULL);
	result.time_sum += time_diff(&t0, &t1);
	return result;
}

WorkerResult final_result;
std::mutex final_result_mutex;

void work() {
	WorkerResult result = measure_keygen();

	// acquire mutex:
	{
		const std::lock_guard<std::mutex> lock(final_result_mutex);
		final_result.combine(result);
	}
}


int main()
{
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	create_huffman_tree();

	const int nthreads = 1;
	if (nthreads == 1) {
		work();
	} else {
		std::thread* pool[nthreads];
		for (int i = 0; i < nthreads; i++) {
			pool[i] = new std::thread(work);
		}

		for (int i = 0; i < nthreads; i++) {
			pool[i]->join();
			delete pool[i];
		}
	}

	double kg_duration = final_result.time_sum / final_result.num_iters; // (in us)
	printf("Average time per keygen: %.3f ms\n", kg_duration / 1000.0);

	printf("# iterations = %lld\n", final_result.num_iters);
	printf("# kgEnc   = %lld (%.1f B)\n", final_result.kgEnc,   (double)final_result.kgEnc  / final_result.num_iters);
	printf("# kgHuf   = %lld (%.1f B)\n", final_result.kgHuf,  (double)final_result.kgHuf / final_result.num_iters);
	printf("# enumEnc = %lld (%.1f B)\n", final_result.enumEnc,  (double)final_result.enumEnc / final_result.num_iters);
	printf("# enumHuf = %lld (%.1f B)\n", final_result.enumHuf, (double)final_result.enumHuf/ final_result.num_iters);
	return 0;
}

