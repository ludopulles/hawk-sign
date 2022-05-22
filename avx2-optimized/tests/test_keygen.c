/*
 * Estimates the probability that keygen fails and has to start over due to
 * various reasons.
 *
 * Moreover, it prints the size of the compressed public key, including the
 * size when using Huffman encoding compared to the compressed-gaussian
 * technique.
 */
#include <assert.h>
#include <stdio.h>
// x86_64 specific:
#include <sys/time.h>

#include "../keygen.c"

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

void poly_output(fpr *p, size_t logn) {
	for (size_t u = 0; u < MKN(logn); u++) {
		if (u) printf(" ");
		printf("%ld", fpr_rint(p[u]));
	}
	printf("\n");
}


#define MAX_Q00 (512) // ~6*sigma
#define MAXLEN_Q00 (52) // max path in the tree
// #define MAX_Q10 (2048) // sufficient enough
// #define MAXLEN_Q10 (21) // max path in the tree
#define MAX_Q10 (4096)
#define MAXLEN_Q10 (56)

struct {
	uint16_t a[MAX_Q00][2]; // { left child, right child }
	uint16_t p[2*MAX_Q00]; // parent
} huffman_tree_q00;
struct {
	uint16_t a[MAX_Q10][2]; // { left child, right child }
	uint16_t p[2*MAX_Q10]; // parent
} huffman_tree_q10;

void init_huffman_trees() {
	float freq[2*MAX_Q10];
	int len[2*MAX_Q10];

#define BUILD_TREE(T, N, L, sigma)                                       \
	/* calculate PDF of normal distribution */                           \
	memset(len, 0, sizeof len);                                          \
	for (int x = 0; x < N; x++)                                          \
		freq[N + x] = exp((float)-x * x / (2.0 * sigma * sigma));        \
                                                                         \
	/* construct the tree */                                             \
	for (uint16_t node = N; --node >= 1; ) {                             \
		/* find 2 nodes with smallest frequencies */                     \
		uint16_t l = 0, r = 0;                                           \
		for (uint16_t idx = 2*N-1; idx > node; --idx) {                  \
			if (freq[idx] < 0) continue;                                 \
			if (!l || freq[idx] < freq[l]) r = l, l = idx;               \
			else if (!r || freq[idx] < freq[r]) r = idx;                 \
		}                                                                \
		/* hide frequency */                                             \
		freq[node] = freq[l] + freq[r];                                  \
		freq[l] = freq[r] = -1;                                          \
		len[node] = 1 + (len[l] > len[r] ? len[l] : len[r]);             \
		T.p[l] = T.p[r] = node;                                          \
		T.a[node][0] = r;                                                \
		T.a[node][1] = l;                                                \
	}                                                                    \
	assert(len[1] <= L && "Longest codeword is too long!");              \
	if (len[1] != L) printf("Longest codeword has length %d\n", len[1]) // ;

	BUILD_TREE(huffman_tree_q00, MAX_Q00, MAXLEN_Q00, 64.0);
	BUILD_TREE(huffman_tree_q10, MAX_Q10, MAXLEN_Q10, 512.0);

#define DEBUG_TREE(T, N, L) \
	for (uint16_t x = 0; x < N; x++) { \
		size_t steps = 0; \
		memset(len, 0, sizeof len); \
		for (int16_t idx = x + N; idx > 1; ) { \
			int16_t next_idx = T.p[idx]; \
			len[steps++] = T.a[next_idx][1] == idx ? 1 : 0; \
			idx = next_idx; \
		} \
		assert(steps <= L); \
		printf("%d: ", x); \
		while (steps --> 0) printf("%d", len[steps]); \
		printf("\n"); \
	}

	// DEBUG_TREE(huffman_tree_q00, MAX_Q00, MAXLEN_Q00);
	// DEBUG_TREE(huffman_tree_q10, MAX_Q10, MAXLEN_Q10);
}

size_t Zf(huffman_encode_q00)(void *out, size_t max_out_len, const int16_t *x, unsigned logn) {
	uint8_t *buf = (uint8_t *)out;
	size_t n = MKN(logn), u, v = 0;
	uint8_t acc = 0, acc_len = 0, steps[MAXLEN_Q00];

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

	/*
	 * The first value of q00 is of the order sigma_kg^2 2d,
	 * but definitely << 8d when sigma_kg = 1.425, by the standard
	 * Laurent-Massart bound (see https://en.wikipedia.org/wiki/Chi-squared_distribution#Concentration).
	 * Here, be lazy and print the whole of x[0].
	 */

	// output x[0] directly, without using acc
	if (buf != NULL) {
		if (max_out_len < 2) return 0;
		buf[0] = ((uint16_t)x[0]) >> 8;
		buf[1] = (uint8_t)x[0];
	}
	v += 2;

	for (u = 1; u < n/2; u ++)
		if (x[u] <= -MAX_Q00 || x[u] >= MAX_Q00) return 0;
	for (u = 1; u < n/2; u ++) {
		uint16_t t;

		ADDBIT(x[u] >> 15); // push the sign bit
		t = (uint16_t)(x[u] < 0 ? (-x[u]) : x[u]); // absolute value

		// store the steps to go up the tree in the buffer
		size_t nsteps = 0;
		for (int16_t idx = MAX_Q00 + t; idx > 1; ) {
			int16_t next_idx = huffman_tree_q00.p[idx];
			steps[nsteps++] = huffman_tree_q00.a[next_idx][1] == idx ? 1 : 0;
			idx = next_idx;
		}

		// print the bits in reverse order, i.e. from top to bottom
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

size_t Zf(huffman_encode_q10)(void *out, size_t max_out_len, const int16_t *x, unsigned logn) {
	uint8_t *buf = (uint8_t *)out;
	size_t n = MKN(logn), u, v = 0;
	uint8_t acc = 0, acc_len = 0, steps[MAXLEN_Q10];

	for (u = 0; u < n; u ++)
		if (x[u] <= -MAX_Q10 || x[u] >= MAX_Q10) return 0;
	for (u = 0; u < n; u ++) {
		uint16_t t;

		ADDBIT(x[u] >> 15); // push the sign bit
		t = (uint16_t)(x[u] < 0 ? (-x[u]) : x[u]); // absolute value

		size_t nsteps = 0;
		for (int16_t idx = MAX_Q10 + t; idx > 1; ) {
			int16_t next_idx = huffman_tree_q10.p[idx];
			steps[nsteps++] = huffman_tree_q10.a[next_idx][1] == idx ? 1 : 0;
			idx = next_idx;
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

#undef ADDBIT

/* see inner.h */
size_t
Zf(comp_encode_q00)(
	void *out, size_t max_out_len,
	const int16_t *x, unsigned logn)
{
	uint8_t *buf;
	size_t n, u, v;
	uint32_t acc;
	unsigned acc_len;

	n = (size_t)1 << logn;
	buf = (uint8_t *)out;

	/*
	 * Make sure that all values are within the -2047..+2047 range.
	 * Moreover, q00 = q00^* so q00 + q00^* must be in fact an integer.
	 */
	for (u = 1; u < n/2; u ++) {
		if (x[u] < -511 || x[u] > +511 || x[u] + x[n-u] != 0) {
			return 0;
		}
	}

	acc = 0;
	acc_len = 0;
	v = 0;

	/*
	 * The first value of q00 is of the order sigma_kg^2 2d,
	 * but definitely << 8d when sigma_kg = 1.425, by the standard
	 * Laurent-Massart bound (see https://en.wikipedia.org/wiki/Chi-squared_distribution#Concentration).
	 */
	if (buf != NULL) {
		if (max_out_len < 2) return 0;
		buf[0] = (uint8_t)x[0];
		buf[1] = ((uint16_t)x[0]) >> 8;
	}
	v += 2;

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

/* see inner.h */
size_t
Zf(comp_encode_q10)(
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
	 * Make sure that all values are within the -4095..+4095 range.
	 */
	for (u = 0; u < n; u ++) {
		if (x[u] < -4095 || x[u] > +4095) {
			return 0;
		}
	}

	acc = 0;
	acc_len = 0;
	v = 0;
	for (u = 0; u < n; u ++) {
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

		const int lo_bits = 8;

		acc <<= lo_bits;
		acc |= w & ((1U << lo_bits) - 1);
		w >>= lo_bits;

		/*
		 * We pushed exactly `lo_bits + 1` bits.
		 */
		acc_len += (lo_bits + 1);

		// TODO: perhaps this still works, but perhaps we need uint64_t...

		/*
		 * Push as many zeros as necessary, then a one. Since the
		 * absolute value is at most 4095, w can only range up to
		 * 7 at this point, thus we will add at most 8 bits
		 * here. With the 10 bits above and possibly up to 7 bits
		 * from previous iterations, we may go up to 25 bits, which
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

/* see inner.h for keygen */
void
keygen_count_fails(inner_shake256_context *rng,
	int8_t *restrict f, int8_t *restrict g,
	int8_t *restrict F, int8_t *restrict G,
	fpr *restrict q00, fpr *restrict q10, fpr *restrict q11,
	unsigned logn, uint8_t *restrict tmp,
	int *normfg_even_fails, int *normfg_fails,
	int *norminvq00_fails, int *NTRU_fails)
{
	/*
	 * Algorithm is the following:
	 *
	 *  - Generate f and g with the Gaussian distribution.
	 *
	 *  - If either N(f) or N(g) is even, try again.
	 *
	 *  - Solve the NTRU equation fG - gF = 1; if the solving fails, try again.
	 *    Usual failure condition is when N(f) and N(g) are not coprime.
	 *
	 *  - Use Babai Reduction on F, G.
	 *
	 *  - Calculate the Gram matrix of the basis [[f, g], [F, G]].
	 */
	size_t n, u;
	uint8_t fg_okay;
	int32_t norm;
	prng p;
	fpr *rt1, *rt2, *rt3;

	n = MKN(logn);
	rt1 = (fpr *)tmp;
	rt2 = rt1 + n;
	rt3 = rt2 + n;

	for (;;) {
		/*
		 * The coefficients of polynomials f and g will be generated from a
		 * discrete gaussian that draws random numbers from a fast PRNG that is
		 * seeded from a SHAKE context ('rng').
		 */
		Zf(prng_init)(&p, rng);

		/*
		 * The coefficients of f and g are generated independently of each
		 * other, with a discrete Gaussian distribution of standard deviation
		 * 1.500. The expected l2-norm of (f, g) is 2n sigma^2.
		 *
		 * We require that N(f) and N(g) are both odd (the binary GCD in the
		 * NTRU solver requires it), so we require (fg_okay & 1) == 1.
		 */
		fg_okay = poly_small_mkgauss(&p, f, logn)
			& poly_small_mkgauss(&p, g, logn) & 1u;

		if (fg_okay == 0u) {
			(*normfg_even_fails)++;
			continue;
		}

		Zf(int8_to_fft)(rt2, f, logn);
		Zf(int8_to_fft)(rt3, g, logn);
		Zf(poly_invnorm2_fft)(rt1, rt2, rt3, logn);
		Zf(iFFT)(rt1, logn);

		if (logn == 9) {
			/*
			 * For n = 512, we reject a key pair if cst(1/q00) >= 0.001, as the
			 * failure probability of decompressing a signature is bounded from
			 * above by 1.9e-32 < 2^{-105}.  Experimentally this fails with
			 * probability of 9%.
			 */
			fg_okay &= fpr_lt(rt1[0], fpr_inv(fpr_of(1000)));
			if (fg_okay == 0u) {
				(*norminvq00_fails)++;
				continue;
			}
		}

		/*
		 * If the l2-norm of (f, g) is shorter than sigma_sec^2 * 2n, BKZ may
		 * return a shortest vector when given the public key much faster than
		 * other instances, so this private key is not secure to use.
		 * Thus, set fg_okay to 0 when ||(f, g)||^2 < Zf(l2bound)[logn]/4.
		 */
		norm = 0;
		for (u = 0; u < n; u++) {
			norm += (int32_t)f[u] * (int32_t)f[u];
			norm += (int32_t)g[u] * (int32_t)g[u];
		}
		norm -= (int32_t)(Zf(l2bound)[logn] >> 2);
		fg_okay &= ((uint32_t) -norm) >> 31;

		if (fg_okay == 0u) {
			(*normfg_fails)++;
			continue;
		}

		assert(fg_okay == 1u);

		/*
		 * Try to solve the NTRU equation for polynomials f and g, i.e. find
		 * polynomials F, G that satisfy
		 *
		 *     f * G - g * F = 1 (mod X^n + 1).
		 */
		if (!solve_NTRU(logn, F, G, f, g, 127, (uint32_t *)tmp)) {
			(*NTRU_fails)++;
			continue;
		}

		Zf(make_public)(f, g, F, G, q00, q10, q11, logn, tmp);
		/*
		 * A valid key pair is generated.
		 */
		break;
	}
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
const size_t logn = 9, n = MKN(logn);

void measure_keygen() {
	uint8_t b[48 << logn];
	int8_t f[n], g[n], F[n], G[n];
	fpr q00[n], q10[n], q11[n];
	int16_t q00n[n], q10n[n];
	unsigned char seed[48];
	inner_shake256_context sc;

	struct timeval t0, t1;
	const int n_repetitions = 1000;

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);

	size_t tot_h00 = 0, tot_c00 = 0;
	size_t tot_h10 = 0, tot_c10 = 0;
	size_t sq_h00 = 0, sq_c00 = 0;
	size_t sq_h10 = 0, sq_c10 = 0;

	int normfg_even_fails = 0, normfg_fails = 0;
	int norminvq00_fails = 0, NTRU_fails = 0;
	for (int i = 0; i < n_repetitions; i++) {
		// Generate key pair.
		keygen_count_fails(&sc, f, g, F, G, q00, q10, q11, logn, b,
			&normfg_even_fails, &normfg_fails, &norminvq00_fails, &NTRU_fails);

		Zf(fft_to_int16)(q00n, q00, logn);
		Zf(fft_to_int16)(q10n, q10, logn);

		size_t pubkey_sz_hq00 = Zf(huffman_encode_q00)(NULL, 0, q00n, logn);
		size_t pubkey_sz_hq10 = Zf(huffman_encode_q10)(NULL, 0, q10n, logn);
		size_t pubkey_sz_cq00 = Zf(comp_encode_q00)(NULL, 0, q00n, logn);
		size_t pubkey_sz_cq10 = Zf(comp_encode_q10)(NULL, 0, q10n, logn);

		if (!pubkey_sz_hq00 || !pubkey_sz_hq10 || !pubkey_sz_cq00 ||!pubkey_sz_cq10) {
			printf("Encoding failed at step %d\n", i);
			i--; continue;
		}

		tot_h00 += pubkey_sz_hq00;
		tot_c00 += pubkey_sz_cq00;
		tot_h10 += pubkey_sz_hq10;
		tot_c10 += pubkey_sz_cq10;
		sq_h00 += pubkey_sz_hq00*pubkey_sz_hq00;
		sq_c00 += pubkey_sz_cq00*pubkey_sz_cq00;
		sq_h10 += pubkey_sz_hq10*pubkey_sz_hq10;
		sq_c10 += pubkey_sz_cq10*pubkey_sz_cq10;
	}

	gettimeofday(&t1, NULL);
	double kg_duration = (double)time_diff(&t0, &t1) / n_repetitions; // (in us)
	printf("Average time per keygen: %.3f ms\n", kg_duration / 1000.0);

	/*
	 * Crunch analysis on what may fail during basis completion, once f, g are
	 * generated.
	 */
	int samples = n_repetitions
		+ normfg_even_fails + normfg_fails
		+ norminvq00_fails + NTRU_fails;
	printf("Pr[ N(f) or N(g) even       ] = %.2f%%\n",
		100.0 * normfg_even_fails / samples);
	printf("Pr[ || (f,g) ||^2 too small ] = %.2f%%\n",
		100.0 * normfg_fails / samples);
	printf("Pr[ cst(1/q00) too large    ] = %.2f%%\n",
		100.0 * norminvq00_fails / samples);
	printf("Pr[ NTRU_solve fails        ] = %.2f%%\n",
		100.0 * NTRU_fails / samples);
	printf("Pr[ keygen works            ] = %.2f%%\n",
		100.0 * n_repetitions / samples);

	double avg, std;
	printf("\nType | |pk| (#bytes) (h = huffman, c = falcon-compression)\n");

	avg = (double) tot_h00 / n_repetitions;
	std = sqrt( (double) sq_h00 / n_repetitions - avg*avg );
	printf("hq00 | %.1f (%.1f)\n", avg, std);
	avg = (double) tot_h10 / n_repetitions;
	std = sqrt( (double) sq_h10 / n_repetitions - avg*avg );
	printf("hq10 | %.1f (%.1f)\n", avg, std);

	avg = (double) tot_c00 / n_repetitions;
	std = sqrt( (double) sq_c00 / n_repetitions - avg*avg );
	printf("cq00 | %.1f (%.1f)\n", avg, std);
	avg = (double) tot_c10 / n_repetitions;
	std = sqrt( (double) sq_c10 / n_repetitions - avg*avg );
	printf("cq10 | %.1f (%.1f)\n", avg, std);
}

void report_invq00_fail_prob() {
	const int n_repetitions = 100 * 1000;

	int8_t f[n], g[n];
	fpr rt1[n], rt2[n], rt3[n];
	unsigned char seed[48];
	inner_shake256_context sc;
	prng p;

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);
	Zf(prng_init)(&p, &sc);

	int num_fails = 0;
	for (int i = 0; i < n_repetitions; i++) {
		poly_small_mkgauss(&p, f, logn);
		poly_small_mkgauss(&p, g, logn);

		Zf(int8_to_fft)(rt2, f, logn);
		Zf(int8_to_fft)(rt3, g, logn);
		Zf(poly_invnorm2_fft)(rt1, rt2, rt3, logn);
		Zf(iFFT)(rt1, logn);

		num_fails += !fpr_lt(rt1[0], fpr_inv(fpr_of(1000)));
	}

	printf("\nExperimental probability that cst(1/q00) >= 0.001 is %.6f\n",
		((double) num_fails) / n_repetitions);
}


int main() {
	init_huffman_trees();
	measure_keygen();
	report_invq00_fail_prob();
	return 0;
}
