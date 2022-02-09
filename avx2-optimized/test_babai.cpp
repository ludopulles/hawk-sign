#include <assert.h>
#include <stdio.h>
// x86_64 specific:
#include <sys/time.h>

extern "C" {
	#ifndef restrict
		#define restrict
	#endif

	#include "inner.h"
}

#include "inner.h"

#include <vector>
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

static void
smallints_to_fpr(fpr *r, const int8_t *t, unsigned logn)
{
	size_t n = MKN(logn), u;
	for (u = 0; u < n; u ++) {
		r[u] = fpr_of(t[u]);
	}
}

void fpr_to_int16(int16_t *buf, fpr *p, size_t logn) {
	for (size_t u = 0, n = MKN(logn); u < n; u ++) {
		int val = fpr_rint(p[u]);
		assert(-(1 << 15) <= val && val < (1 << 15));
		buf[u] = (int16_t) val;
	}
}

unsigned sqnorm(int8_t *p, size_t logn) {
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
#define MAX_Q00 (512) // ~6*sigma
#define ENCODING_LEN_Q00 (96) // max path in the tree
#define MAX_Q10 (2*2048)
#define MAX_Q10_LARGE (MAX_Q10)
#define ENCODING_LEN_Q10 (56)

static struct { uint16_t a[MAX_Q00][2], p[2*MAX_Q00]; } tree_q00;
static struct { uint16_t a[MAX_Q10][2], p[2*MAX_Q10]; } tree_q10;
static struct { uint16_t a[MAX_Q10_LARGE][2], p[2*MAX_Q10_LARGE]; } tree_q10_large;

static
void create_huffman_tree() {
	float freq[2*MAX_Q10_LARGE];
	uint16_t u, l, r, v;

#define BUILD_TREE(T, N, sigma)                                          \
	/* calculate PDF of normal distribution */                           \
	for (u = 0; u < N; u++)                                              \
		freq[N + u] = exp((float)-u * u / (2.0 * sigma * sigma));        \
	/* construct the tree */                                             \
	for (u = N; --u >= 1; ) {                                            \
		l = r = 0; /* find 2 nodes with smallest frequencies */          \
		for (v = 2*N; --v > u; ) {                                       \
			if (freq[v] < 0) continue; /* v is already used */           \
			if (!l || freq[v] < freq[l]) r = l, l = v;                   \
			else if (!r || freq[v] < freq[r]) r = v;                     \
		}                                                                \
		freq[u] = freq[l] + freq[r];                                     \
		freq[l] = freq[r] = -1; /* mark l and r as used */               \
		T.p[l] = T.p[r] = u;                                             \
		T.a[u][0] = r;                                                   \
		T.a[u][1] = l;                                                   \
	}

	BUILD_TREE(tree_q00, MAX_Q00, 45.75);
	BUILD_TREE(tree_q10, MAX_Q10, 512.0);
	BUILD_TREE(tree_q10_large, MAX_Q10_LARGE, 512.0);
}

static size_t
Zf(huffman_encode)(void *out, size_t max_out_len, const int16_t *x,
	unsigned logn)
{
	uint8_t *buf = (uint8_t *)out;
	size_t n = MKN(logn), u, v = 0;
	uint8_t acc = 0, acc_len = 0, steps[ENCODING_LEN_Q10];

#define ADDBIT(x) {                                                      \
	acc = (acc << 1) | (x);                                              \
	if (++acc_len == 8) {                                                \
		if (buf != NULL) {                                               \
			if (max_out_len <= v) return 0;                              \
			buf[v] = acc;                                                \
		}                                                                \
		acc_len = acc = 0; /* reset acc */                               \
		v++;                                                             \
	}                                                                    \
}

	for (u = 0; u < n; u ++)
		if (x[u] <= -MAX_Q10 || x[u] >= MAX_Q10) return 0;

	/*
	 * Then output q10 using the second Huffman tree.
	 */
	for (u = 0; u < n; u ++) {
		uint16_t t, s;
		size_t nsteps;

		ADDBIT(x[u] >> 15); // push the sign bit
		t = (uint16_t)(x[u] < 0 ? (-x[u]-1) : x[u]); // absolute value
		nsteps = 0; // store the steps to go up the tree in the buffer
		for (t += MAX_Q10; t > 1; t = s) {
			s = tree_q10.p[t];
			steps[nsteps++] = (tree_q10.a[s][1] == t);
		}

		// print the bits in reverse order, i.e. from root to leaf
		while (nsteps --> 0) {
			ADDBIT(steps[nsteps]);
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

static size_t
Zf(huffman_encode_large)(void *out, size_t max_out_len, const int16_t *x,
	unsigned logn)
{
	uint8_t *buf = (uint8_t *)out;
	size_t n = MKN(logn), u, v = 0;
	uint8_t acc = 0, acc_len = 0, steps[128];

#define ADDBIT(x) {                                                      \
	acc = (acc << 1) | (x);                                              \
	if (++acc_len == 8) {                                                \
		if (buf != NULL) {                                               \
			if (max_out_len <= v) return 0;                              \
			buf[v] = acc;                                                \
		}                                                                \
		acc_len = acc = 0; /* reset acc */                               \
		v++;                                                             \
	}                                                                    \
}

	for (u = 0; u < n; u ++)
		if (x[u] <= -MAX_Q10_LARGE || x[u] >= MAX_Q10_LARGE) return 0;

	/*
	 * Then output q10 using the second Huffman tree.
	 */
	for (u = 0; u < n; u ++) {
		uint16_t t, s;
		size_t nsteps;

		ADDBIT(x[u] >> 15); // push the sign bit
		t = (uint16_t)(x[u] < 0 ? (-x[u]-1) : x[u]); // absolute value
		nsteps = 0; // store the steps to go up the tree in the buffer
		for (t += MAX_Q10_LARGE; t > 1; t = s) {
			s = tree_q10_large.p[t];
			steps[nsteps++] = (tree_q10_large.a[s][1] == t);
		}

		// print the bits in reverse order, i.e. from root to leaf
		while (nsteps --> 0) {
			ADDBIT(steps[nsteps]);
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


/** see codec.c, Zf(comp_encode) */
size_t
Zf(encode_q10)(void *out, size_t max_out_len, const int16_t *x,
	unsigned logn, const int lim, const int lo_bits)
{
	uint8_t *buf = (uint8_t *)out;
	size_t n = MKN(logn), u, v = 0;
	uint64_t acc = 0;
	unsigned acc_len = 0;

	for (u = 0; u < n; u ++)
		if (x[u] <= -lim || x[u] >= lim) return 0;
	for (u = 0; u < n; u ++) {
		int t;
		unsigned w;

		// Get sign and absolute value of next integer; push the sign bit.
		acc <<= 1;
		t = x[u];
		if (t < 0) t = -t-1, acc |= 1;
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




static size_t
Zf(huffman_encode_q00)(void *out, size_t max_out_len, const int16_t *x,
	unsigned logn)
{
	uint8_t *buf = (uint8_t *)out;
	size_t n = MKN(logn), u, v = 0;
	uint8_t acc = 0, acc_len = 0, steps[ENCODING_LEN_Q00];

// #define ADDBIT(x) [...]
	for (u = 1; u <= n/2; u ++)
		if (x[u] <= -MAX_Q00 || x[u] >= MAX_Q00 || x[u] + x[n-u] != 0) return 0;

	if (buf != NULL) {
		if (2 >= max_out_len) return 0;
		buf[0] = (uint8_t)(x[0]);
		buf[1] = (uint8_t)((uint16_t)x[0] >> 8);
	}
	v = 2;

	for (u = 1; u < n/2; u ++) {
		uint16_t t, s;
		ADDBIT(x[u] >> 15); // push the sign bit
		t = (uint16_t)(x[u] < 0 ? (-x[u]-1) : x[u]); // absolute value
		size_t nsteps = 0; // store the steps to go up the tree in the buffer
		for (t += MAX_Q00; t > 1; t = s) {
			s = tree_q00.p[t];
			steps[nsteps++] = (tree_q00.a[s][1] == t);
		}

		while (nsteps --> 0) ADDBIT(steps[nsteps]);
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


size_t
Zf(encode_q00)(void *out, size_t max_out_len, const int16_t *x,
	unsigned logn, const int lim, const int lo_bits)
{
	uint8_t *buf = (uint8_t *)out;
	size_t n = MKN(logn), u, v = 0;
	uint64_t acc = 0;
	unsigned acc_len = 0;

	for (u = 1; u <= n/2; u ++)
		if (x[u] <= -lim || x[u] >= lim || x[u] + x[n-u] != 0) return 0;

	if (buf != NULL) {
		if (2 >= max_out_len) return 0;
		buf[0] = (uint8_t)(x[0]);
		buf[1] = (uint8_t)((uint16_t)x[0] >> 8);
	}
	v = 2;

	for (u = 1; u < n/2; u ++) {
		int t;
		unsigned w;

		// Get sign and absolute value of next integer; push the sign bit.
		acc <<= 1;
		t = x[u];
		if (t < 0) t = -t-1, acc |= 1;
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

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
const size_t logn = 9, n = MKN(logn);

struct WorkerResult
{
	long long num_iters, num_fails;
	long long enc, huf, encB, hufB;
	long long sqenc, sqhuf, sqencB, sqhufB;
	long long enc00, huf00, sqenc00, sqhuf00;
	double time_sum;

	WorkerResult() : num_iters(0), num_fails(0),
		enc(0), huf(0), encB(0), hufB(0),
		sqenc(0), sqhuf(0), sqencB(0), sqhufB(0),
		enc00(0), huf00(0), sqenc00(0), sqhuf00(0),
		time_sum(0.0) {}

	void combine(const WorkerResult &res) {
		num_iters += res.num_iters;
		num_fails += res.num_fails;

		enc += res.enc;
		huf += res.huf;
		encB += res.encB;
		hufB += res.hufB;

		sqenc += res.sqenc;
		sqhuf += res.sqhuf;
		sqencB += res.sqencB;
		sqhufB += res.sqhufB;

		enc00 += res.enc00;
		huf00 += res.huf00;
		sqenc00 += res.sqenc00;
		sqhuf00 += res.sqhuf00;

		time_sum += res.time_sum;
	}
};

WorkerResult measure_keygen(fpr isigma_kg)
{
	union {
		uint8_t b[64*512];
		uint64_t dummy_i64;
		fpr dummy_fpr;
	} tmp;
	int8_t f[n], g[n], F[n], G[n];
	fpr *_f = (fpr*)tmp.b, *_g = _f + n, *_F = _g + n, *_G = _F + n;
	fpr q00[n], q10[n], q11[n];
	int16_t q00i[n], q10i[n];
	unsigned char seed[48];
	inner_shake256_context sc;
	struct timeval t0, t1;

	WorkerResult result;
	result.num_iters = 2500;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);
	for (int _ = 0; _ < result.num_iters; _++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, isigma_kg, logn, tmp.b);

		Zf(iFFT)(q00, logn);
		Zf(iFFT)(q10, logn);
		fpr_to_int16(q00i, q00, logn);
		fpr_to_int16(q10i, q10, logn);

		size_t enc00 = Zf(encode_q00)(NULL, 0, q00i, logn, 512, 5);
		size_t huf00 = Zf(huffman_encode_q00)(NULL, 0, q00i, logn);
		size_t enc = Zf(encode_q10)(NULL, 0, q10i, logn, 4096, 8);
		size_t huf = Zf(huffman_encode_large)(NULL, 0, q10i, logn);

		if (!enc00 || !huf00 || !enc || !huf) { printf("."); fflush(stdout); result.num_fails++; _--; continue; }

		smallints_to_fpr(_f, f, logn);
		smallints_to_fpr(_g, g, logn);
		smallints_to_fpr(_F, F, logn);
		smallints_to_fpr(_G, G, logn);
		Zf(FFT)(_f, logn); // f
		Zf(FFT)(_g, logn); // g
		Zf(FFT)(_F, logn); // F
		Zf(FFT)(_G, logn); // G

		Zf(ffBabai_reduce)(_f, _g, _F, _G, F, G, logn, _G + n);

		// q10 = F*adj(f) + G*adj(g)
		Zf(poly_add_muladj_fft)(q10, _F, _G, _f, _g, logn);
		Zf(iFFT)(q10, logn);
		fpr_to_int16(q10i, q10, logn);

		size_t encB = Zf(encode_q10)(NULL, 0, q10i, logn, 2048, 8);
		size_t hufB = Zf(huffman_encode)(NULL, 0, q10i, logn);

		if (!encB || !hufB) { result.num_fails++; _--; continue; }

		result.enc += enc;
		result.sqenc += enc*enc;
		result.huf += huf;
		result.sqhuf += huf*huf;
		result.encB += encB;
		result.sqencB += encB*encB;
		result.hufB += hufB;
		result.sqhufB += hufB*hufB;

		result.enc00 += enc00;
		result.sqenc00 += enc00*enc00;
		result.huf00 += huf00;
		result.sqhuf00 += huf00*huf00;
	}

	gettimeofday(&t1, NULL);
	result.time_sum += time_diff(&t0, &t1);
	return result;
}

int8_t valid_sigma(fpr sigma_sig)
{
	return !fpr_lt(sigma_sig, fpr_sigma_min[logn])
		&& fpr_lt(sigma_sig, fpr_div(fpr_of(18205), fpr_of(10000)));
}

constexpr fpr sigma_kg  = { v: 1.425 };

WorkerResult final_result;
std::mutex final_result_mutex;

void work() {
	WorkerResult result = measure_keygen(fpr_inv(sigma_kg));

	// acquire mutex:
	{
		const std::lock_guard<std::mutex> lock(final_result_mutex);
		final_result.combine(result);
	}
}

const int nthreads = 4;

int main()
{
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	assert(valid_sigma(sigma_kg));

	create_huffman_tree();

	std::vector<std::thread*> pool(nthreads - 1);
	for (int i = 0; i < nthreads - 1; i++)
		pool[i] = new std::thread(work);
	work();
	for (std::thread *t : pool)
		t->join(), delete t;

	double kg_duration = final_result.time_sum / final_result.num_iters; // (in us)
	printf("Average time per keygen: %.3f ms\n", kg_duration / 1000.0);

	printf("# iterations  = %lld\n", final_result.num_iters);
	printf("# fails       = %lld\n", final_result.num_fails);


	double avg, sqavg;
#define PRINTF(fstr, sum, sqsum)                                         \
	avg = (double) (sum) / final_result.num_iters;                       \
	sqavg = (double) (sqsum) / final_result.num_iters;                   \
	printf(fstr, avg, sqrt(sqavg - avg*avg))

	PRINTF("# enc         = %.1f±%.1f B\n", final_result.enc,  final_result.sqenc);
	PRINTF("# huf         = %.1f±%.1f B\n", final_result.huf,  final_result.sqhuf);
	PRINTF("# enc (babai) = %.1f±%.1f B\n", final_result.encB, final_result.sqencB);
	PRINTF("# huf (babai) = %.1f±%.1f B\n", final_result.hufB, final_result.sqhufB);
	PRINTF("# enc q00     = %.1f±%.1f B\n", final_result.enc00, final_result.sqenc00);
	PRINTF("# huf q00     = %.1f±%.1f B\n", final_result.huf00, final_result.sqhuf00);

	PRINTF("# pubkey      = %.1f±%.1f B\n", final_result.encB + final_result.enc00, final_result.sqencB + final_result.sqenc00);
	PRINTF("# pubkey (huf)= %.1f±%.1f B\n", final_result.hufB + final_result.huf00, final_result.sqhufB + final_result.sqhuf00);
	return 0;
}

