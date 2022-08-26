/*
 * Hawk pk, sk and sig compression.
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2022 Hawk Project
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
 * @author   Ludo Pulles <ludo.pulles@cwi.nl>
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

#include "inner.h"

#define LOW_BITS_S1(logn) ((logn) == 10 ? 6U : 5U)

static const size_t low_bits_q00[11] = {
	0 /* unused */, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6
};
static const size_t low_bits_q01[11] = {
	0 /* unused */, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
};

/* see inner.h */
size_t
Zf(encode_pubkey)(void *out, size_t max_out_len,
	const int16_t *q00, const int16_t *q01, unsigned logn)
{
	uint8_t *buf;
	size_t n, u, v;
	uint16_t w, low_mask;
	int16_t bound;
	uint64_t acc;
	unsigned acc_len;

	n = MKN(logn);
	buf = (uint8_t *)out;
	low_mask = (1U << low_bits_q00[logn]) - 1;
	bound = (int16_t)(1U << Zf(bits_q00)[logn]);

	/*
	 * Make sure no coefficient is too large.
	 */
	for (u = 1; u < n/2; u ++) {
		if (q00[u] < -bound || q00[u] >= bound) {
			return 0;
		}
	}

	bound = (int16_t)(1U << Zf(bits_q01)[logn]);
	/*
	 * Make sure no coefficient is too large.
	 */
	for (u = 0; u < n; u ++) {
		if (q01[u] < -bound || q01[u] >= bound) {
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
	 * surely fit in an uint16_t. Thus, print q00[0] in little-endian format.
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
		 *
		 * HAWK-1024: Since the initial w is < 2^10, w can only range up to 15
		 * at this point, thus we will add at most 16 bits here. With possibly
		 * up to 7 bits from previous iterations and the 6+1 bits above, we
		 * may go up to 30 bits, which will fit in an uint64_t accumulator.
		 *
		 * HAWK-512: w is <2^9/2^5 at this point and we get 5+1 bits from
		 * above, so we may go up to 29 bits here.
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
	 * Encode q01.
	 */
	low_mask = (1U << low_bits_q01[logn]) - 1;
	for (u = 0; u < n; u ++) {
		/*
		 * Push the sign bit and store |x| - [x<0] in w.
		 * Note that x = 1,0,-1,-2,... give value for w of 1,0,0,1,... resp.
		 */
		w = (uint16_t) q01[u];
		acc = (acc << 1) | (w >> 15);
		w ^= -(w >> 15);

		/*
		 * Push the low 8 bits of w.
		 */
		acc = (acc << low_bits_q01[logn]) | (w & low_mask);
		w >>= low_bits_q01[logn];
		acc_len += low_bits_q01[logn] + 1;

		/*
		 * Push as many zeros as necessary, then a one.
		 *
		 * HAWK-1024: Since the initial w is < 2^14, w can only range up to 15
		 * at this point, thus we will add at most 16 bits here. With possibly
		 * up to 7 bits from previous iterations and the 10+1 bits above, we
		 * may go up to 34 bits, which will fit in an uint64_t accumulator.
		 *
		 * HAWK-512: w is <2^12/2^9 at this point and we get 9+1 bits from
		 * above, so we may go up to 25 bits here.
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

#define ENSUREBIT()                                                           \
	if (acc_len == 0) {                                                       \
		if (v >= max_in_len) return 0; /* not enough bits */                  \
		acc = buf[v++];                                                       \
		acc_len = 8;                                                          \
	}

#define GETBIT() ((acc >> (--acc_len)) & 1)

/* see inner.h */
size_t
Zf(decode_pubkey)(int16_t *q00, int16_t *q01,
	const void *in, size_t max_in_len, unsigned logn)
{
	const uint8_t *buf;
	size_t n, u, v;
	uint64_t acc;
	uint16_t low_mask, high_inc;
	int16_t bound;
	unsigned acc_len;

	buf = (uint8_t *)in;
	n = MKN(logn);
	acc = acc_len = 0;
	high_inc = (1U << low_bits_q00[logn]);
	low_mask = high_inc - 1;
	bound = (int16_t)(1U << Zf(bits_q00)[logn]);

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
		ENSUREBIT();
		while (GETBIT() == 0) {
			w += high_inc;
			if (w >= bound) {
				return 0;
			}
			ENSUREBIT();
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
	 * Decode q01.
	 */
	high_inc = (1U << low_bits_q01[logn]);
	low_mask = high_inc - 1;
	bound = (int16_t)(1U << Zf(bits_q01)[logn]);
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
		if (acc_len < low_bits_q01[logn]) {
			// should be true all the time
			if (v >= max_in_len) return 0;
			acc = (acc << 8) | buf[v ++];
			acc_len += 8;

			// Note: low_bits_q01[logn] may be 9.
			if (acc_len < low_bits_q01[logn]) {
				// should be true all the time
				if (v >= max_in_len) return 0;
				acc = (acc << 8) | buf[v ++];
				acc_len += 8;
			}
		}

		/*
		 * Get 8 least significant bits of w.
		 */
		w = (acc >> (acc_len -= low_bits_q01[logn])) & low_mask;

		/*
		 * Recover the most significant bits of w: count number of consecutive
		 * zeros up to the first 1 and add 2^low_bits_q01[logn] to w for each
		 * one you see.
		 */
		ENSUREBIT();
		while (GETBIT() == 0) {
			w += high_inc;
			if (w >= bound) {
				return 0;
			}
			ENSUREBIT();
		}

		q01[u] = w ^ -s;
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
Zf(encode_sig)(void *out, size_t max_out_len, const int16_t *s1, unsigned logn)
{
	size_t n, u, v;
	uint8_t *buf;
	int16_t bound;
	uint32_t acc;
	unsigned acc_len;

	n = MKN(logn);
	buf = (uint8_t *)out;

	/*
	 * Make sure no coefficient is too large.
	 */
	bound = (1U << Zf(bits_s1)[logn]) - 1;
	for (u = 0; u < n; u ++) {
		if (s1[u] < -bound || s1[u] >= bound) return 0;
	}


	acc = 0;
	acc_len = 0;
	v = 0;
	for (u = 0; u < n; u ++) {
		uint16_t w;

		/*
		 * Push sign bit
		 */
		w = (uint16_t) s1[u];
		acc = (acc << 1) | (w >> 15);
		w ^= -(w >> 15);

		/*
		 * Push the lowest `lo_bits` bits of w which is equal to |x| - [x < 0].
		 */
		acc <<= LOW_BITS_S1(logn);
		acc |= w & ((1U << LOW_BITS_S1(logn)) - 1);
		w >>= LOW_BITS_S1(logn);
		acc_len += LOW_BITS_S1(logn) + 1;

		/*
		 * Push as many zeros as necessary, then a one.
		 * 
		 * HAWK-1024: Since the initial w is < 2^10, w can only range up to 15
		 * at this point, thus we will add at most 16 bits here. With possibly
		 * up to 7 bits from previous iterations and the 6+1 bits from above, we
		 * may go up to 30 bits, which will fit in an uint32_t accumulator.
		 *
		 * HAWK-512: w is <2^9/2^5 at this point and we get 5+1 bits from
		 * above, so we may go up to 29 bits here.
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

size_t
Zf(decode_sig)(int16_t *s1, const void *in, size_t max_in_len, unsigned logn)
{
	const uint8_t *buf;
	size_t n, u, v, nbits;
	int16_t bound;
	uint16_t mask;
	uint32_t acc;
	unsigned acc_len;

	buf = (uint8_t *)in;
	n = MKN(logn);
	acc = 0;
	acc_len = 0;
	v = 0;
	nbits = 1 + LOW_BITS_S1(logn);
	bound = (1U << Zf(bits_s1)[logn]) - 1;
	mask = (1U << LOW_BITS_S1(logn)) - 1;

	for (u = 0; u < n; u ++) {
		uint16_t s, w;

		/*
		 * Check space in buffer
		 */
		if (acc_len < nbits) {
			if (v >= max_in_len) return 0;
			acc = (acc << 8) | buf[v ++];
			acc_len += 8;
		}

		/*
		 * Get sign of the current coefficient
		 */
		s = (acc >> (acc_len - 1)) & 1;

		/*
		 * Get lo_bits least significant bits of w.
		 */
		acc_len -= nbits;
		w = (acc >> acc_len) & mask;

		/*
		 * Recover the most significant bits of w: count number of consecutive
		 * zeros up to the first 1 and add 2^lo_bits to w for each one you see.
		 */
		ENSUREBIT();
		while (GETBIT() == 0) {
			w += mask + 1;

			/*
			 * Make sure no coefficient is too large.
			 */
			if (w >= bound) return 0;
			ENSUREBIT();
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
}

/*
 * Encode f, g with 5 bits and F with 8 bits.
 * Note: G can be recalculated from fG - gF = 1 so do not encode it.
 */
size_t
Zf(encode_seckey)(uint8_t *out, size_t max_out_len,
	const int8_t *f, const int8_t *g, const int8_t *F, unsigned logn)
{
	uint16_t acc, mask;
	size_t n, u, out_len, acc_len, lobits;

	n = MKN(logn);
	acc = 0;
	acc_len = 0;

	lobits = logn == 10 ? 6 : 5;
	mask = 1u << (lobits - 1);
	for (u = 0; u < n; u++) {
		if (f[u] < -mask || f[u] >= mask || g[u] < -mask || g[u] >= mask) {
			return 0;
		}
	}

	out_len = ((n * (8 + 2 * lobits)) + 7) >> 3;
	mask = (1u << lobits) - 1;

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
		acc = (acc << lobits) | ((uint8_t) x[u] & mask);                      \
		acc_len += lobits;                                                    \
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
	uint16_t acc, mask;
	size_t n, u, in_len, acc_len, lobits;

	n = MKN(logn);
	acc = 0;
	acc_len = 0;

	lobits = logn == 10 ? 6 : 5;
	in_len = ((n * (8 +  2 * lobits)) + 7) >> 3;
	mask = (1u << lobits) - 1;

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
		if (acc_len < lobits) {                                               \
			acc = (acc << 8) | (uint16_t)(*in ++);                            \
			acc_len += 8;                                                     \
		}                                                                     \
		/* Take 5 bits from the accumulator */                                \
		w = (uint8_t)(acc >> (acc_len -= lobits)) & mask;                     \
		/* Let bits 7..5 match (sign) bit 4 */                                \
		w -= (w << 1) & (mask + 1);                                           \
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

#undef ENSUREBIT
#undef GETBIT
