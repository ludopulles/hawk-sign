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

#include <math.h>

#include "inner.h"


#define MAX_Q00 (512) // ~6*sigma
#define ENCODING_LEN_Q00 (52) // max path in the tree
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

#define BUILD_TREE(T, N, sigma)                                               \
	/* calculate PDF of normal distribution */                                \
	for (u = 0; u < N; u++)                                                   \
		freq[N + u] = exp((float)-u * u / (2.0 * sigma * sigma));             \
	/* construct the tree */                                                  \
	for (u = N; --u >= 1; ) {                                                 \
		l = r = 0; /* find 2 nodes with smallest frequencies */               \
		for (v = 2*N; --v > u; ) {                                            \
			if (freq[v] < 0) continue; /* v is already used */                \
			if (!l || freq[v] < freq[l]) r = l, l = v;                        \
			else if (!r || freq[v] < freq[r]) r = v;                          \
		}                                                                     \
		freq[u] = freq[l] + freq[r];                                          \
		freq[l] = freq[r] = -1; /* mark l and r as used */                    \
		T.p[l] = T.p[r] = u;                                                  \
		T.a[u][0] = r;                                                        \
		T.a[u][1] = l;                                                        \
	}

	BUILD_TREE(tree_q00, MAX_Q00, 64.00);
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
Zf(encode_pubkey_huffman)(void *out, size_t max_out_len,
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
#define ADDBIT(x) {                                                           \
	acc = (acc << 1) | (x);                                                   \
	if (++acc_len == 8) {                                                     \
		if (buf != NULL) {                                                    \
			if (max_out_len <= v) return 0;                                   \
			buf[v] = acc;                                                     \
		}                                                                     \
		acc_len = acc = 0; /* reset acc */                                    \
		v++;                                                                  \
	}                                                                         \
}

	buf = (uint8_t *)out;
	n = MKN(logn);
	v = 0;
	acc = acc_len = 0;

	/*
	 * The constant coefficient of q00 follows a chi-squared distribution, has
	 * a mean of roughly sigma_pk^2 2n and is quite concentrated by [1] so will
	 * surely fit in a uint16_t. Thus, print q00[0] in little-endian format.
	 *
	 * [1] https://en.wikipedia.org/wiki/Chi-squared_distribution#Concentration
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
Zf(decode_pubkey_huffman)(int16_t *q00, int16_t *q10,
	const void *in, size_t max_in_len, unsigned logn)
{
	const uint8_t *buf;
	size_t n, u, v;
	uint8_t acc, acc_len;

#define ENSUREBIT()                                                           \
	if (acc_len == 0) {                                                       \
		if (max_in_len <= v) return 0; /* not enough bits */                  \
		acc = buf[v++];                                                       \
		acc_len = 8;                                                          \
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
	q00[n/2] = 0;
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
#undef ENSUREBIT
#undef GETBIT
}

/* see inner.h */
size_t
Zf(encode_sig_huffman)(void *out, size_t max_out_len,
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
#define ADDBIT(x) {                                                           \
	acc = (acc << 1) | (x);                                                   \
	if (++acc_len == 8) {                                                     \
		if (buf != NULL) {                                                    \
			if (max_out_len <= v) return 0;                                   \
			buf[v] = acc;                                                     \
		}                                                                     \
		acc_len = acc = 0; /* reset acc */                                    \
		v++;                                                                  \
	}                                                                         \
}


	buf = (uint8_t *)out;
	n = MKN(logn);
	v = 0;
	acc = acc_len = 0;

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


/* =============================================================================
 * Encoding/decoding that will be used with the signature scheme, as the gain of
 * a Huffman table is not significant but makes the code more complex.
 */

#define MAX_ALLOWED_Q00(logn) (1 << ((28-logn)/2))
// n = 1,2: 2^13
// n = 3,4: 2^12
// n = 5,6: 2^11
// n = 7,8: 2^10
// n = 9: 2^9

static const size_t low_bits_q00[10] = {
	0 /* unused */, 1, 2, 2, 3, 3, 4, 4, 5, 5
};
static const size_t low_bits_q10[10] = {
	0 /* unused */, 1, 2, 3, 4, 5, 6, 7, 8, 9
};

/* see inner.h */
size_t
Zf(encode_pubkey)(void *out, size_t max_out_len,
	const int16_t *q00, const int16_t *q10, unsigned logn)
{
	uint8_t *buf;
	size_t n, u, v;
	uint16_t w, low_mask;
	uint64_t acc;
	unsigned acc_len;

	n = MKN(logn);
	buf = (uint8_t *)out;
	low_mask = (1U << low_bits_q00[logn]) - 1;

	/*
	 * Make sure coefficient 1 up to n/2 of q00 are within the -512..+511
	 * range. Moreover, we silently assume q00 is self adjoint.
	 */
	for (u = 1; u < n/2; u ++) {
		if (q00[u] < -MAX_ALLOWED_Q00(logn) || q00[u] >= MAX_ALLOWED_Q00(logn)) {
			return 0;
		}
	}

	/*
	 * Make sure all coefficients of q10 are within the -4096..+4095 range.
	 */
	for (u = 0; u < n; u ++) {
		if (q10[u] < -4096 || q10[u] >= 4096) {
			return 0;
		}
	}

	acc = 0;
	acc_len = 0;
	v = 0;


	/*
	 * Encode q00.
	 * The constant coefficient of q00 follows a chi-squared distribution, has
	 * a mean of roughly sigma_pk^2 2n and is quite concentrated by [1] so will
	 * surely fit in a uint16_t. Thus, print q00[0] in little-endian format.
	 *
	 * [1] https://en.wikipedia.org/wiki/Chi-squared_distribution#Concentration
	 */
	if (buf != NULL) {
		if (max_out_len < 2) return 0;
		buf[0] = (uint8_t)q00[0];
		buf[1] = ((uint16_t)q00[0]) >> 8;
	}
	v += 2;

	for (u = 1; u < n/2; u ++) {
		/*
		 * Push the sign bit and store |x| - [x<0] in w.
		 * Note that x = 1,0,-1,-2,... give value for w of 1,0,0,1,... resp.
		 */
		w = (uint16_t) q00[u];
		acc = (acc << 1) | (w >> 15);
		w ^= -(w >> 15);

		/*
		 * Push the low 5 bits of w, when logn = 9.
		 */
		acc = (acc << low_bits_q00[logn]) | (w & low_mask);
		w >>= low_bits_q00[logn];
		acc_len += low_bits_q00[logn] + 1;

		/*
		 * Push as many zeros as necessary, then a one.
		 * For n = 512: Since the absolute value is at most 511, w can only
		 * range up to 15 at this point, thus we will add at most 16 bits here.
		 * With the 6 bits above and possibly up to 7 bits from previous
		 * iterations, we may go up to 28 bits, which will fit in the
		 * accumulator, which is an uint32_t.
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
	 * Encode q10.
	 */
	low_mask = (1U << low_bits_q10[logn]) - 1;
	for (u = 0; u < n; u ++) {
		/*
		 * Push the sign bit and store |x| - [x<0] in w.
		 * Note that x = 1,0,-1,-2,... give value for w of 1,0,0,1,... resp.
		 */
		w = (uint16_t) q10[u];
		acc = (acc << 1) | (w >> 15);
		w ^= -(w >> 15);

		/*
		 * Push the low 8 bits of w.
		 */
		acc = (acc << low_bits_q10[logn]) | (w & low_mask);
		w >>= low_bits_q10[logn];
		acc_len += low_bits_q10[logn] + 1;

		/*
		 * Push as many zeros as necessary, then a one. Since the initial w is
		 * at most 4095, w can only range up to 15 at this point, thus we will
		 * add at most 16 bits here. With the 9 bits above and possibly up to 7
		 * bits from previous iterations, we may go up to 32 bits, which will
		 * fit in the accumulator, which is an uint32_t.
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
Zf(decode_pubkey)(int16_t *q00, int16_t *q10,
	const void *in, size_t max_in_len, unsigned logn)
{
	const uint8_t *buf;
	size_t n, u, v;
	uint64_t acc;
	uint16_t low_mask, high_inc;
	unsigned acc_len;

#define ENSUREBIT()                                                           \
	if (acc_len == 0) {                                                       \
		if (v >= max_in_len) return 0; /* not enough bits */                  \
		acc = buf[v++];                                                       \
		acc_len = 8;                                                          \
	}

#define GETBIT() ((acc >> (--acc_len)) & 1)

	buf = (uint8_t *)in;
	n = MKN(logn);
	acc = acc_len = 0;
	high_inc = (1U << low_bits_q00[logn]);
	low_mask = high_inc - 1;

	/*
	 * Decode q00.
	 * First read q00[0].
	 */
	if (max_in_len < 2) {
		return 0;
	}
	q00[0] = ((uint16_t)buf[1] << 8) | (uint16_t) buf[0];
	v = 2;

	/*
	 * Read q00[1] ... q00[n/2 - 1].
	 */
	for (u = 1; u < n/2; u ++) {
		uint16_t s, w;

		/*
		 * Get sign of the current coefficient
		 */
		ENSUREBIT();
		s = GETBIT();

		if (acc_len < low_bits_q00[logn]) {
			if (v >= max_in_len) {
				return 0;
			}
			acc = (acc << 8) | buf[v ++];
			acc_len += 8;
		}

		/*
		 * Get 5 least significant bits of w.
		 */
		w = (acc >> (acc_len -= low_bits_q00[logn])) & low_mask;

		/*
		 * Recover the most significant bits of w: count number of consecutive
		 * zeros up to the first 1 and add 2^5 to w for each one you see.
		 */
		for (;;) {
			ENSUREBIT();
			if (GETBIT() != 0) {
				break;
			}
			w += high_inc;
			if (w >= MAX_ALLOWED_Q00(logn)) {
				return 0;
			}
		}

		q00[u] = w ^ -s;
	}

	/*
	 * Since q00 is self-adjoint, we can recover q00[n / 2] ... q00[n - 1]
	 * now.
	 */
	q00[n/2] = 0;
	for (u = n/2 + 1; u < n; u ++) {
		q00[u] = -q00[n - u];
	}

	/*
	 * Decode q10.
	 */
	high_inc = (1U << low_bits_q10[logn]);
	low_mask = high_inc - 1;
	for (u = 0; u < n; u ++) {
		uint16_t s, w;

		/*
		 * Get sign of the current coefficient
		 */
		ENSUREBIT();
		s = GETBIT();

		/*
		 * Get next eight bits that make up the lowest significant bits of w.
		 */
		if (acc_len < low_bits_q10[logn]) {
			// should be true all the time
			if (v >= max_in_len) return 0;
			acc = (acc << 8) | buf[v ++];
			acc_len += 8;

			// Note: low_bits_q10[logn] may be 9.
			if (acc_len < low_bits_q10[logn]) {
				// should be true all the time
				if (v >= max_in_len) return 0;
				acc = (acc << 8) | buf[v ++];
				acc_len += 8;
			}
		}

		/*
		 * Get 8 least significant bits of w.
		 */
		w = (acc >> (acc_len -= low_bits_q10[logn])) & low_mask;

		/*
		 * Recover the most significant bits of w: count number of consecutive
		 * zeros up to the first 1 and add 2^low_bits_q10[logn] to w for each
		 * one you see.
		 */
		for (;;) {
			ENSUREBIT();
			if (GETBIT() != 0) {
				break;
			}
			w += high_inc;
			if (w >= 4096) {
				return 0;
			}
		}

		q10[u] = w ^ -s;
	}

	/*
	 * Unused bits in the last byte must be zero.
	 */
	if ((acc & ((1u << acc_len) - 1u)) != 0) {
		return 0;
	}

	return v;
#undef ENSUREBIT
#undef GETBIT
}

/* see inner.h */
size_t
Zf(encode_sig_simple)(void *out, size_t max_out_len,
	const int16_t *s0, const int16_t *s1, unsigned logn,
	size_t lo_bits_s0, size_t lo_bits_s1)
{
	uint8_t *buf;
	uint16_t w;
	size_t n, u, v;
	uint64_t acc;
	unsigned acc_len;

	n = MKN(logn);
	buf = (uint8_t *)out;

	/*
	 * Make sure that all values are within the -512..+511 range.
	 */
	for (u = 0; u < n; u ++)
		if (s0[u] < -2048 || s0[u] >= 2048) return 0;
	for (u = 0; u < n; u ++)
		if (s1[u] < -512 || s1[u] >= 512) return 0;

	acc = 0;
	acc_len = 0;
	v = 0;

	for (u = 0; u < n; u ++) {
		/*
		 * Push sign bit
		 */
		w = (uint16_t) s0[u];
		acc = (acc << 1) | (w >> 15);
		w ^= -(w >> 15);

		/*
		 * Push the lowest `lo_bits` bits of w which is equal to |x|'.
		 */
		acc <<= lo_bits_s0;
		acc |= w & ((1U << lo_bits_s0) - 1);
		w >>= lo_bits_s0;
		acc_len += lo_bits_s0 + 1;

		/*
		 * Push as many zeros as necessary, then a one.
		 */
		acc <<= (w + 1);
		acc |= 1;
		acc_len += w + 1;
		if (acc_len >= 64) return 0; // There is an overflow

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

	for (u = 0; u < n; u ++) {
		/*
		 * Push sign bit
		 */
		w = (uint16_t) s1[u];
		acc = (acc << 1) | (w >> 15);
		w ^= -(w >> 15);

		/*
		 * Push the lowest `lo_bits` bits of w which is equal to |x|'.
		 */
		acc <<= lo_bits_s1;
		acc |= w & ((1U << lo_bits_s1) - 1);
		w >>= lo_bits_s1;
		acc_len += lo_bits_s1 + 1;

		/*
		 * Push as many zeros as necessary, then a one.
		 */
		acc <<= (w + 1);
		acc |= 1;
		acc_len += w + 1;
		if (acc_len >= 64) return 0; // There is an overflow

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

size_t
Zf(decode_sig_simple)(int16_t *s0, int16_t *s1,
	const void *in, size_t max_in_len, unsigned logn,
	size_t lo_bits_s0, size_t lo_bits_s1)
{
	const uint8_t *buf;
	uint16_t s, w;
	size_t n, u, v;
	uint64_t acc;
	unsigned acc_len;

#define ENSUREBIT()                                                           \
	if (acc_len == 0) {                                                       \
		if (v >= max_in_len) return 0; /* not enough bits */                  \
		acc = buf[v++];                                                       \
		acc_len = 8;                                                          \
	}

#define GETBIT() ((acc >> (--acc_len)) & 1)

	buf = (uint8_t *)in;
	n = MKN(logn);
	acc = 0;
	acc_len = 0;
	v = 0;

	for (u = 0; u < n; u ++) {
		/*
		 * Get sign of the current coefficient
		 */
		ENSUREBIT();
		s = GETBIT();

		/*
		 * Get next lo_bits bits that make up the lowest significant bits of w.
		 */
		if (acc_len < lo_bits_s0) {
			// should be true all the time
			if (v >= max_in_len) return 0;
			acc = (acc << 8) | buf[v ++];
			acc_len += 8;
		}

		/*
		 * Get lo_bits least significant bits of w.
		 */
		w = (acc >> (acc_len -= lo_bits_s0)) & ((1U << lo_bits_s0) - 1);

		/*
		 * Recover the most significant bits of w: count number of consecutive
		 * zeros up to the first 1 and add 2^lo_bits to w for each one you see.
		 */
		for (;;) {
			ENSUREBIT();
			if (GETBIT() != 0) break;
			w += (1U << lo_bits_s0);
			if (w >= 2048) return 0;
		}

		s0[u] = w ^ -s;
	}

	for (u = 0; u < n; u ++) {
		/*
		 * Get sign of the current coefficient
		 */
		ENSUREBIT();
		s = GETBIT();

		/*
		 * Get next lo_bits bits that make up the lowest significant bits of w.
		 */
		if (acc_len < lo_bits_s1) {
			// should be true all the time
			if (v >= max_in_len) return 0;
			acc = (acc << 8) | buf[v ++];
			acc_len += 8;
		}

		/*
		 * Get lo_bits least significant bits of w.
		 */
		w = (acc >> (acc_len -= lo_bits_s1)) & ((1U << lo_bits_s1) - 1);

		/*
		 * Recover the most significant bits of w: count number of consecutive
		 * zeros up to the first 1 and add 2^lo_bits to w for each one you see.
		 */
		for (;;) {
			ENSUREBIT();
			if (GETBIT() != 0) break;
			w += (1U << lo_bits_s1);
			if (w >= 512) return 0;
		}

		s1[u] = w ^ -s;
	}

	/*
	 * Unused bits in the last byte must be zero.
	 */
	if ((acc & ((1u << acc_len) - 1u)) != 0) {
		return 0;
	}

	return v;
#undef ENSUREBIT
#undef GETBIT
}

/* see inner.h */
size_t
Zf(encode_sig)(void *out, size_t max_out_len, const int16_t *x, unsigned logn,
	size_t lo_bits)
{
	uint8_t *buf;
	size_t n, u, v;
	uint64_t acc;
	unsigned acc_len;

	n = MKN(logn);
	buf = (uint8_t *)out;

	/*
	 * Make sure that all values are within the -512..+511 range.
	 */
	for (u = 0; u < n; u ++)
		if (x[u] < -512 || x[u] >= 512) return 0;

	acc = 0;
	acc_len = 0;
	v = 0;
	for (u = 0; u < n; u ++) {
		uint16_t w;

		/*
		 * Push sign bit
		 */
		w = (uint16_t) x[u];
		acc = (acc << 1) | (w >> 15);
		w ^= -(w >> 15);

		/*
		 * Push the lowest `lo_bits` bits of w which is equal to |x| - [x < 0].
		 */
		acc <<= lo_bits;
		acc |= w & ((1U << lo_bits) - 1);
		w >>= lo_bits;
		acc_len += lo_bits + 1;

		/*
		 * TODO: assume lo_bits = 5 in the thing below.
		 * Push as many zeros as necessary, then a one. Since the initial w is
		 * at most 511, w can only range up to 15 at this point, thus we will
		 * add at most 16 bits here. With the 6 bits above and possibly up to 7
		 * bits from previous iterations, we may go up to 29 bits, which will
		 * fit in the accumulator, which is an uint32_t.
		 */
		acc <<= (w + 1);
		acc |= 1;
		acc_len += w + 1;

		// TODO: remove when this check cannot fail anymore.
		if (acc_len >= 64) return 0; // There is an overflow

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

size_t
Zf(decode_sig)(int16_t *x, const void *in, size_t max_in_len, unsigned logn,
	size_t lo_bits)
{
	const uint8_t *buf;
	size_t n, u, v;
	uint16_t acc;
	unsigned acc_len;

#define ENSUREBIT()                                                           \
	if (acc_len == 0) {                                                       \
		if (v >= max_in_len) return 0; /* not enough bits */                  \
		acc = buf[v++];                                                       \
		acc_len = 8;                                                          \
	}

#define GETBIT() ((acc >> (--acc_len)) & 1)

	buf = (uint8_t *)in;
	n = MKN(logn);
	acc = 0;
	acc_len = 0;
	v = 0;

	for (u = 0; u < n; u ++) {
		uint16_t s, w;

		/*
		 * Get sign of the current coefficient
		 */
		ENSUREBIT();
		s = GETBIT();

		/*
		 * Get next lo_bits bits that make up the lowest significant bits of w.
		 */
		if (acc_len < lo_bits) {
			// should be true all the time
			if (v >= max_in_len) return 0;
			acc = (acc << 8) | buf[v ++];
			acc_len += 8;
		}

		/*
		 * Get lo_bits least significant bits of w.
		 */
		w = (acc >> (acc_len -= lo_bits)) & ((1U << lo_bits) - 1);

		/*
		 * Recover the most significant bits of w: count number of consecutive
		 * zeros up to the first 1 and add 2^lo_bits to w for each one you see.
		 */
		for (;;) {
			ENSUREBIT();
			if (GETBIT() != 0) {
				break;
			}
			w += (1U << lo_bits);
			if (w >= 512) {
				return 0;
			}
		}

		x[u] = w ^ -s;
	}

	/*
	 * Unused bits in the last byte must be zero.
	 */
	if ((acc & ((1u << acc_len) - 1u)) != 0) {
		return 0;
	}

	return v;
#undef ENSUREBIT
#undef GETBIT
}

/*
 * Encode f, g with 5 bits and F with 8 bits.
 * Note: G can be recalculated from fG - gF = 1 so do not encode it.
 */
size_t
Zf(encode_seckey)(uint8_t *out, size_t max_out_len,
	const int8_t *f, const int8_t *g, const int8_t *F, unsigned logn)
{
	uint16_t acc;
	size_t n, u, out_len, acc_len;

	n = MKN(logn);
	acc = 0;
	acc_len = 0;
	out_len = ((n * 18) + 7) >> 3;

	for (u = 0; u < n; u++) {
		if (f[u] < -16 || f[u] >= 16 || g[u] < -16 || g[u] >= 16) {
			return 0;
		}
	}

	if (out == NULL) {
		return out_len;
	} else if (out_len > max_out_len) {
		return 0;
	}

	/*
	 * First write F literally as it is.
	 */
	for (u = 0; u < n; u++) {
		*out ++ = (uint8_t)F[u];
	}

	/*
	 * Now write f, g with 5 bits each.
	 */
#define POLY_ENC(x)                                                           \
	for (u = 0; u < n; u ++) {                                                \
		/* Push exactly len bits */                                           \
		acc = (acc << 5) | ((uint8_t) x[u] & 31u);                            \
		acc_len += 5;                                                         \
		/* Produce all full bytes. */                                         \
		if (acc_len >= 8) {                                                   \
			acc_len -= 8;                                                     \
			*out ++ = (uint8_t)(acc >> acc_len);                              \
		}                                                                     \
	}

	POLY_ENC(f);
	POLY_ENC(g);

	/*
	 * Flush remaining bits (if any).
	 */
	if (acc_len > 0) {
		*out ++ = (uint8_t)(acc << (8 - acc_len));
	}
	return out_len;
#undef POLY_ENC
}

/*
 * Decodes the f, g, and F from the encoded secret key basis [[f,g], [F,G]].
 * Moreover, G can be efficiently constructed from fG - gF = 1 (mod phi) using
 * FFT.
 */
size_t
Zf(decode_seckey)(int8_t *f, int8_t *g, int8_t *F,
	const uint8_t *in, size_t max_in_len, unsigned logn)
{
	uint8_t w;
	uint16_t acc;
	size_t n, u, in_len, acc_len;

	n = MKN(logn);
	acc = 0;
	acc_len = 0;
	in_len = ((n * 18) + 7) >> 3;

	if (in_len > max_in_len) {
		return 0;
	}

	/*
	 * First read F literally as it is.
	 */
	for (u = 0; u < n; u++) {
		F[u] = (int8_t)(*in ++);
	}

	/*
	 * Now read f, g with 5 bits each.
	 */
#define POLY_DEC(x)                                                           \
	for (u = 0; u < n; u ++) {                                                \
		/* Fill the accumulator if needed */                                  \
		if (acc_len < 5) {                                                    \
			acc = (acc << 8) | (uint16_t)(*in ++);                            \
			acc_len += 8;                                                     \
		}                                                                     \
		/* Take 5 bits from the accumulator */                                \
		w = (uint8_t)(acc >> (acc_len -= 5)) & 31u;                           \
		/* Let bits 7..5 match (sign) bit 4 */                                \
		w -= (w & 16) << 1;                                                   \
		x[u] = (int8_t)w;                                                     \
	}

	POLY_DEC(f);
	POLY_DEC(g);

	/*
	 * Extra bits in the last byte must be zero.
	 */
	if ((acc & ((1u << acc_len) - 1u)) != 0) {
		return 0;
	}
	return in_len;
#undef POLY_DEC
}
