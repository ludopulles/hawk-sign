/*
 * Hawk pk, sk and sig compression 
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2017-2019  Falcon Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ===========================(LICENSE END)=============================
 *
 * @author   Thomas Pornin <thomas.pornin@nccgroup.com>
 */

#include "math.h"

#include "inner.h"


#define MAX_Q00 (512) // ~6*sigma
#define ENCODING_LEN_Q00 (96) // max path in the tree
// #define MAX_Q10 (2048) // sufficient enough
// #define ENCODING_LEN_Q10 (21) // max path in the tree
#define MAX_Q10 (4096)
#define ENCODING_LEN_Q10 (56)

#define MAX_S1 (512)
#define ENCODING_LEN_S1 (64)

// a[i]: { left child of i, right child of i }
// p[i]: parent of node i
static struct {
	uint16_t a[MAX_Q00][2], p[2*MAX_Q00];
} tree_q00 = { { {0, 0}}, {} };
static struct {
	uint16_t a[MAX_Q10][2], p[2*MAX_Q10];
} tree_q10 = { { {0, 0}}, {} };
static struct {
	uint16_t a[MAX_S1][2], p[2*MAX_S1];
} tree_s1 = {{{0, 0}},{}};


static void
init_huffman_trees() {
	// TODO: make this work with multiple threads
	if (tree_q00.a[0][0]) {
		// initialization is already performed
		return;
	}

	// initialize tree_q00 and tree_q10

	float freq[2*MAX_Q10];
	uint16_t u, l, r, v;
	/* int len[2*MAX_Q10]; */

#define BUILD_TREE(T, N, sigma)                                          \
	/* calculate PDF of normal distribution */                           \
	/* memset(len, 0, sizeof len); */                                    \
	for (u = 0; u < N; u++) {                                            \
		freq[N + u] = exp((float)-u * u / (2.0 * sigma * sigma));        \
    }                                                                    \
	/* construct the tree */                                             \
	for (u = N; --u >= 1; ) {                                            \
		/* find 2 nodes with smallest frequencies */                     \
		l = r = 0;                                                       \
		for (v = 2*N; --v > u; ) {                                       \
			if (freq[v] < 0) continue;                                   \
			if (!l || freq[v] < freq[l]) r = l, l = v;                   \
			else if (!r || freq[v] < freq[r]) r = v;                     \
		}                                                                \
		freq[u] = freq[l] + freq[r];                                     \
		/* hide frequency */                                             \
		freq[l] = freq[r] = -1;                                          \
		/* len[node] = 1 + (len[l] > len[r] ? len[l] : len[r]); */       \
		T.p[l] = T.p[r] = u;                                             \
		T.a[u][0] = r;                                                   \
		T.a[u][1] = l;                                                   \
	}                                                                    \
	/* printf("Longest codeword has length %d\n", len[1]); */

	BUILD_TREE(tree_q00, MAX_Q00, 45.75);
	BUILD_TREE(tree_q10, MAX_Q10, 512.0);

	// This has a too large std.dev. to be optimal:
	// BUILD_TREE(tree_s1, MAX_S1, 58.0);
	BUILD_TREE(tree_s1, MAX_S1, 64.0);

	/* size_t len[2*MAX_Q10];
	memset(len, 0, sizeof len);
	for (u = 2*MAX_S1; --u > 1; )
		if (len[u]+1 > len[tree_s1.p[u]]) len[tree_s1.p[u]] = len[u]+1;
	printf("Longest length: %u\n", len[1]); */

/*
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

	DEBUG_TREE(tree_q00, MAX_Q00, ENCODING_LEN_Q00);
	DEBUG_TREE(tree_q10, MAX_Q10, ENCODING_LEN_Q10);
*/
	
	// mark the tree generation as done:
	tree_q00.a[0][0] = 1;
}

size_t
Zf(encode_pubkey)(void *out, size_t max_out_len,
	const int16_t *q00, const int16_t *q10, unsigned logn)
{
	uint8_t *buf;
	size_t n, u, v;
	uint8_t acc, acc_len, steps[ENCODING_LEN_Q10];


	/*
	 * Within one byte, the oldest bits are the most significant bits of
	 * the byte, while in the buffer, the oldest bytes are the first ones.
	 * Thus, byte-ordering is little-endian, while bit-ordering is
	 * big-endian.
	 */
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


	buf = (uint8_t *)out;
	n = MKN(logn);
	v = 0;
	acc = acc_len = 0;

	/*
	 * The first value of q00 is of the order sigma_kg^2 2d,
	 * but definitely << 8d when sigma_kg = 1.425, by the standard
	 * Laurent-Massart bound [1].
	 * Here, be lazy and print the whole of x[0].
	 *
	 * [1] en.wikipedia.org/wiki/Chi-squared_distribution#Concentration
	 */

	init_huffman_trees();

	for (u = 1; u < n/2; u ++)
		if (q00[u] <= -MAX_Q00 || q00[u] >= MAX_Q00) return 0;
	for (u = 0; u < n; u ++)
		if (q10[u] <= -MAX_Q10 || q10[u] >= MAX_Q10) return 0;

	/*
	 * First output q00 using q00 is self-adjoint:
	 * - output q00[0] in little-endian format (2 bytes)
	 * - output q00[1] ... q00[n/2 - 1] using the first Huffman tree
	 */
	if (buf != NULL) {
		if (max_out_len < 2) return 0;
		buf[0] = (uint8_t)q00[0];
		buf[1] = ((uint16_t)q00[0]) >> 8;
	}
	v += 2;

	for (u = 1; u < n/2; u ++) {
		uint16_t t, s;
		size_t nsteps;

		ADDBIT(q00[u] >> 15); // push the sign bit
		t = (uint16_t)(q00[u] < 0 ? (-q00[u]) : q00[u]); // absolute value
		nsteps = 0; // store the steps to go up the tree in the buffer
		for (t += MAX_Q00; t > 1; t = s) {
			s = tree_q00.p[t];
			steps[nsteps++] = (tree_q00.a[s][1] == t);
		}

		// print the bits in reverse order, i.e. from root to leaf
		while (nsteps --> 0) {
			ADDBIT(steps[nsteps]);
		}
	}

	/*
	 * Then output q10 using the second Huffman tree.
	 */
	for (u = 0; u < n; u ++) {
		uint16_t t, s;
		size_t nsteps;

		ADDBIT(q10[u] >> 15); // push the sign bit
		t = (uint16_t)(q10[u] < 0 ? (-q10[u]) : q10[u]); // absolute value
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

	// Flush remaining bits (if any) and pad with zeros.
	if (acc_len > 0) {
		if (buf != NULL) {
			if (max_out_len <= v) return 0;
			buf[v] = (uint8_t)(acc << (8 - acc_len));
		}
		v++;
	}
	return v;
}

size_t
Zf(decode_pubkey)(int16_t *q00, int16_t *q10,
	const void *in, size_t max_in_len, unsigned logn)
{
	const uint8_t *buf;
	size_t n, u, v;
	uint8_t acc, acc_len;


#define ENSUREBIT()                                                      \
	if (acc_len == 0) {                                                  \
		if (max_in_len <= v) return 0; /* not enough bits */             \
		acc = buf[v++];                                                  \
		acc_len = 8;                                                     \
	}

#define GETBIT() ((acc >> (--acc_len)) & 1)


	buf = (uint8_t *)in;
	n = MKN(logn);
	v = 0;
	acc = acc_len = 0;

	init_huffman_trees();

	/*
	 * First read q00[0].
	 */
	if (max_in_len < 2) return 0;
	uint16_t q00_0 = ((uint16_t)buf[0] << 8) | (uint16_t) buf[1];
	q00[0] = q00_0;
	v += 2;

	/*
	 * Read q00[1] ... q00[n/2 - 1].
	 */
	for (u = 1; u < n/2; u ++) {
		uint16_t s, val;

		ENSUREBIT();
		s = GETBIT();

		/*
		 * First, val is an index in the Huffman tree. After that, it is
		 * the value of the u'th coefficient of q00.
		 */
		val = 1;
		while (val < MAX_Q00) {
			ENSUREBIT();
			val = tree_q00.a[val][GETBIT()];
		}
		val -= MAX_Q00;

		/*
		 * "-0" is forbidden.
		 */
		if (s && val == 0) {
			return 0;
		}

		q00[u] = s ? -(int16_t)val : (int16_t)val;
	}
	/*
	 * Since q00 is self-adjoint, we can recover q00[n / 2] ... q00[n - 1]
	 * now.
	 */
	q00[n/2] = 0U;
	for (u = n/2 + 1; u < n; u ++)
		q00[u] = -q00[n - u];

	/*
	 * Read q10[0] ... q10[n - 1].
	 */
	for (u = 0; u < n; u ++) {
		uint16_t s, val;

		ENSUREBIT();
		s = GETBIT();

		/*
		 * First, val is an index in the Huffman tree. After that, it is
		 * the value of the u'th coefficient of q00.
		 */
		val = 1;
		while (val < MAX_Q10) {
			ENSUREBIT();
			val = tree_q10.a[val][GETBIT()];
		}
		val -= MAX_Q10;

		/*
		 * "-0" is forbidden.
		 */
		if (s && val == 0) {
			return 0;
		}

		q10[u] = s ? -(int16_t)val : (int16_t)val;
	}

	/*
	 * Unused bits in the last byte must be zero.
	 */
	if ((acc & ((1u << acc_len) - 1u)) != 0) {
		return 0;
	}

	return v;
}

/* see inner.h */
size_t
Zf(encode_sig_huffman)(
	void *out, size_t max_out_len,
	const int16_t *x, unsigned logn)
{
	uint8_t *buf;
	size_t n, u, v;
	uint8_t acc, acc_len, steps[ENCODING_LEN_S1];

	/*
	 * Within one byte, the oldest bits are the most significant bits of
	 * the byte, while in the buffer, the oldest bytes are the first ones.
	 * Thus, byte-ordering is little-endian, while bit-ordering is
	 * big-endian.
	 */
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


	buf = (uint8_t *)out;
	n = MKN(logn);
	v = 0;
	acc = acc_len = 0;

	/*
	 * The first value of q00 is of the order sigma_kg^2 2d,
	 * but definitely << 8d when sigma_kg = 1.425, by the standard
	 * Laurent-Massart bound [1].
	 * Here, be lazy and print the whole of x[0].
	 *
	 * [1] en.wikipedia.org/wiki/Chi-squared_distribution#Concentration
	 */

	init_huffman_trees();

	for (u = 0; u < n; u ++)
		if (x[u] <= -MAX_S1 || x[u] >= MAX_S1) return 0;

	for (u = 0; u < n; u ++) {
		uint16_t t, s;
		size_t nsteps;

		ADDBIT(x[u] >> 15); // push the sign bit
		t = (uint16_t)(x[u] < 0 ? (-x[u]) : x[u]); // absolute value
		nsteps = 0; // store the steps to go up the tree in the buffer
		for (t += MAX_S1; t > 1; t = s) {
			s = tree_s1.p[t];
			steps[nsteps++] = (tree_s1.a[s][1] == t);
		}

		// print the bits in reverse order, i.e. from root to leaf
		while (nsteps --> 0) {
			ADDBIT(steps[nsteps]);
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
Zf(encode_sig)(
	void *out, size_t max_out_len,
	const int16_t *x, unsigned logn, size_t lo_bits)
	// TODO: hardcode lo_bits with the best possible value
{
	uint8_t *buf;
	size_t n, u, v;
	uint32_t acc;
	unsigned acc_len;

	n = (size_t)1 << logn;
	buf = (uint8_t *)out;

	init_huffman_trees();
	// TODO: use Huffman tree

	/*
	 * Make sure that all values are within the -511..+511 range.
	 */
	for (u = 0; u < n; u ++) {
		if (x[u] <= -512 || x[u] >= 512) {
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
		 * Push the lowest `lo_bits` bits of the absolute value.
		 */
		acc <<= lo_bits;
		acc |= w & ((1U << lo_bits) - 1);
		w >>= lo_bits;
		acc_len += lo_bits + 1;

		if (acc_len + w + 1 >= 32) {
			// There will be an overflow
			return 0;
		}
		/*
		 * TODO: fix the numbers here.
		 *
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

	if (lo_bits == 5 && v >= 530) return 0;
	return v;
}
