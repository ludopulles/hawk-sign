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

static void poly_small_to_fp(fpr *x, const int8_t *f, unsigned logn) {
	for (size_t u = 0, n = MKN(logn); u < n; u ++)
		x[u] = fpr_of(f[u]);
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


/** see codec.c, Zf(comp_encode) */
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

/** see codec.c, Zf(comp_encode) */
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

void nearest_plane(const int8_t *f, const int8_t *g, int8_t *F, int8_t *G)
{
	int8_t F_orig[n], G_orig[n];
	fpr tree[(logn - 1)*n/2], tmp[11*n]; // TODO lower 11

	float Btilde[2*n][n], invnorm[n], dot;
	// B = Btilde * mu
	size_t REVi[n];


	memcpy(F_orig, F, sizeof F_orig);
	memcpy(G_orig, G, sizeof G_orig);

	// ffBabai(f, g, F_orig, G_orig, logn, tree, tmp);
	Zf(ffStraightBabai)(f, g, F_orig, G_orig, logn, tmp);

	int sqn = 0;
	for (size_t u = 0; u < n; u++)
		sqn += F[u]*F[u] + G[u]*G[u];
	printf("Square norm %d\n", sqn);

	printf("F, G (ffBabai) = \n");
	sqn = 0;
	for (size_t u = 0; u < n; u++)
		printf("%d ", F_orig[u]), sqn += F_orig[u]*F_orig[u];
	printf("\n");
	for (size_t u = 0; u < n; u++)
		printf("%d ", G_orig[u]), sqn += G_orig[u]*G_orig[u];
	printf("\n");
	printf("Square norm %d\n", sqn);



	for (size_t i = 0; i < n; i++)
		REVi[i] = bitreverse(i, logn);

	for (size_t i = 0; i < n; i++) {
		size_t rev_i = REVi[i];
		// multiply (f, g)^T by X^i
		// put in the ith column of Btilde
		for (size_t j = 0; i + j < n; j++) Btilde[i+j  ][rev_i] =  f[j];
		for (size_t j = n - i; j < n; j++) Btilde[i+j-n][rev_i] = -f[j];
		for (size_t j = 0; i + j < n; j++) Btilde[i+j+n][rev_i] =  g[j];
		for (size_t j = n - i; j < n; j++) Btilde[i+j  ][rev_i] = -g[j];
	}

	// now orthogonalize, O(n^3)
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < i; j++) {
			// mu_{ji} = <b_i, b_j> / <b_j, b_j>,
			dot = 0.0;
			for (size_t k = 0; k < n+n; k++)
				dot += Btilde[k][i] * Btilde[k][j];
			dot *= invnorm[j];
			// b_i -= mu_{ji} b_j
			for (size_t k = 0; k < n+n; k++)
				Btilde[k][i] -= dot * Btilde[k][j];
		}

		invnorm[i] = 0.0;
		for (size_t k = 0; k < n+n; k++)
			invnorm[i] += Btilde[k][i] * Btilde[k][i];
		invnorm[i] = 1.0 / invnorm[i];
	}

	int x[n];
	// Do Babai's algorithm with F, G as target
	for (size_t i = n; i --> 0; ) {
		dot = 0.0;
		for (size_t k = 0; k < n; k++)
			dot += F[k] * Btilde[k  ][i]
				 + G[k] * Btilde[k+n][i];
		int xtake = round(dot * invnorm[i]);

		x[REVi[i]] = xtake;
		// subtract (fx^i, gx^i) from (F, G)
		if (xtake == 0) continue;

		// Btilde[i] corresponds to (f,g)[rev_i]
		int rev_i = REVi[i];
		for (size_t k = 0; rev_i + k < n; k++)
			F[rev_i + k] -= xtake * f[k], G[rev_i + k] -= xtake * g[k];
		// Use X^n = -1!
		for (size_t k = n - rev_i; k < n; k++)
			F[rev_i+k-n] += xtake * f[k], G[rev_i+k-n] += xtake * g[k];
	}

	printf("x = ");
	for (size_t u = 0; u < n; u++) printf("%d ", x[u]);
	printf("\n");
	printf("F, G = \n");
	sqn = 0;
	for (size_t u = 0; u < n; u++)
		printf("%d ", F[u]), sqn += F[u] * F[u];
	printf("\n");
	for (size_t u = 0; u < n; u++)
		printf("%d ", G[u]), sqn += G[u] * G[u];
	printf("\n");
	printf("Square norm %d\n", sqn);
	exit(1);
}

struct WorkerResult {
	int16_t qvalues[4096];

	WorkerResult() {
		memset(qvalues, 0, sizeof qvalues);
	}

	void addResult(int16_t *q10i) {
		for (size_t u = 0; u < n; u++) {
			if (q10i[u] <= -2048 || q10i[u] >= 2048) {
				fprintf(stderr, "Skipping too large value!\n");
				fflush(stderr);
			} else {
				qvalues[ q10i[u] + 2048 ]++;
			}
		}
	}

	void combine(const WorkerResult *res) {
		for (size_t u = 0; u < 4096; u++)
			qvalues[u] += res->qvalues[u];
	}
};

WorkerResult measure_keygen(fpr isigma_kg)
{
	union {
		uint8_t b[28*512]; // 14 kB temporary memory, 17.5 kB total
		uint64_t dummy_i64;
		fpr dummy_fpr;
	} tmp;
	int8_t f[n], g[n], F[n], G[n];
	fpr q00[n], q10[n], q11[n];
	int16_t q00i[n], q10i[n];
	unsigned char seed[48];
	inner_shake256_context sc;

	fpr rt1[n], rt2[n], rt3[n], rt4[n];

	// uint8_t output_buf[1024];
	// struct timeval t0, t1;
	const int n_repetitions = 10;

	WorkerResult result;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	// gettimeofday(&t0, NULL);
	for (int _ = 0; _ < n_repetitions; _++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, tmp.b, isigma_kg);

		int8_t Fp[n], Gp[n];
		memcpy(Fp, F, sizeof F);
		memcpy(Gp, G, sizeof G);
		// now try to optimize Fp, Gp with babai.

		nearest_plane(f, g, Fp, Gp);

		/* 
		printf("Sqnorm saving: %5d -> %5d,\t%5d -> %5d\n",
			sqnorm(F, logn), sqnorm(Fp, logn),
			sqnorm(G, logn), sqnorm(Gp, logn));

		Zf(iFFT)(q00, logn);
		Zf(iFFT)(q10, logn);
		fpr_to_int16(q00i, q00, logn);
		fpr_to_int16(q10i, q10, logn);
		printf("First   q10 bits: %zu\n", Zf(encode_q10)(NULL, 0, q10i, logn, 4095, 9));
		int mv = 0;
		for (size_t u = 0; u < n; u++) {
			if (q10i[u] > mv) mv = q10i[u];
			if (-q10i[u] > mv) mv = -q10i[u];
		}
		printf("Max value: %d\n", mv); */

		// Calculate q10 (in FFT representation)
		poly_small_to_fp(rt1, f, logn);
		poly_small_to_fp(rt2, g, logn);
		poly_small_to_fp(rt3, Fp, logn);
		poly_small_to_fp(rt4, Gp, logn);
		Zf(FFT)(rt1, logn); // g
		Zf(FFT)(rt2, logn); // F
		Zf(FFT)(rt3, logn); // G
		Zf(FFT)(rt4, logn); // f

		// q10 = F*adj(f) + G*adj(g)
		Zf(poly_add_muladj_fft)(q10, rt3, rt4, rt1, rt2, logn);
		Zf(iFFT)(q10, logn);
		fpr_to_int16(q10i, q10, logn);

		result.addResult(q10i);

		printf("Comparison pubkey sizes: %zu (encode_q10) %zu (huffman_encode)\n",
			Zf(encode_q10)(NULL, 0, q10i, logn, 2047, 8),
			Zf(huffman_encode)(NULL, 0, q10i, logn));
		/* printf("Current pubkey sz: %zu + %zu", Zf(encode_q00)(NULL, 0, q00i, logn), Zf(encode_q10)(NULL, 0, q10i, logn, 4095, 9));
		for (size_t lo_bits = 7; lo_bits <= 10; lo_bits++)
			printf(" vs %zu", Zf(encode_q10)(NULL, 0, q10i, logn, 2047, lo_bits));
		printf("\n");
		print_int16(q10i, logn);
		mv = 0;
		for (size_t u = 0; u < n; u++) {
			if (q10i[u] > mv) mv = q10i[u];
			if (-q10i[u] > mv) mv = -q10i[u];
		}
		printf("Max value: %d\n", mv);

		// print compressed q10
		size_t len = Zf(encode_q10)(output_buf, sizeof output_buf, q10i, logn, 2047, 8);
		if (len == 0) {
			printf("Error occurred.\n");
			exit(1);
		}
		size_t curu = 0, curbit = 7, lo_bits = 8;

		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < lo_bits + 1; j++) { // (also count sign-bit)
				uint8_t bit = (output_buf[curu] >> curbit) & 1;
				if (curbit-- == 0) curbit = 7, curu++;
				printf("%u", bit);
			}
			while (1) {
				uint8_t bit = (output_buf[curu] >> curbit) & 1;
				if (curbit-- == 0) curbit = 7, curu++;
				printf("%u", bit);
				if (bit == 1) { printf(" "); break; }
			}
		}

		while (curu < len) {
			uint8_t bit = (output_buf[curu] >> curbit) & 1;
			if (curbit-- == 0) curbit = 7, curu++;
			printf("%u", bit);
		}

		printf("\n"); */

/*
		for (size_t u = 0; u < len; u++) {
			printf("%u%u%u%u%u%u%u%u",
				(output_buf[u]>>7)&1, (output_buf[u]>>6)&1, (output_buf[u]>>5)&1, (output_buf[u]>>4)&1,
				(output_buf[u]>>3)&1, (output_buf[u]>>2)&1, (output_buf[u]>>1)&1, (output_buf[u]>>0)&1);
		} */
	}

	/* gettimeofday(&t1, NULL);
	double kg_duration = (double)time_diff(&t0, &t1) / n_repetitions; // (in us)
	printf("Average time per keygen: %.3f ms\n", kg_duration / 1000.0); */

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
		final_result.combine(&result);
	}
}


int main(int argc, char **argv)
{
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	// seed = 2784228632;
	printf("Seed: %u\n", seed);
	srand(seed);

	assert(valid_sigma(sigma_kg));

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
		}
	}

	if (argc >= 2 && strcmp(argv[1], "values") == 0) {
		printf("values = [");
		for (size_t u = 0; u < 4096; u++) {
			printf("%d, ", final_result.qvalues[u]);
		}
	}
	return 0;
}

