/*
 * Common functions for Hawk
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

/* see inner.h */
/* NIST-1: */
const uint32_t Zf(l2bound_512)[10] = {
	0u /* unused */, 32u, 64u, 129u, 259u, 519u, 1039u, 2079u, 4158u, 8317u
};

/* NIST-5: */
const uint32_t Zf(l2bound_1024)[11] = {
	0u /* unused */, 39u, 79u, 158u, 316u, 632u, 1265u, 2530u, 5060u, 10121u, 20243u
};

/* see inner.h */
const unsigned Zf(bits_q00)[11] = {
	0 /* unused */, 5, 6, 6, 7, 7,  8,  8,  9,  9, 10
};
const unsigned Zf(bits_q01)[11] = {
	0 /* unused */, 5, 6, 7, 8, 8,  9, 10, 11, 12, 14
};
const unsigned Zf(bits_q11)[11] = {
	0 /* unused */, 5, 6, 7, 8, 9, 10, 12, 13, 15, 17
};

/* see inner.h */
const unsigned Zf(bits_s0)[11] = {
	0 /* unused */, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
};
const unsigned Zf(bits_s1)[11] = {
	0 /* unused */, 5, 6, 6, 7, 7,  8,  8,  9,  9, 10
};

/* see inner.h */
void
Zf(int8_to_fft)(fpr *p, const int8_t *x, unsigned logn)
{
	size_t u;

	u = MKN(logn);
	while (u -- > 0) {
		p[u] = fpr_of(x[u]);
	}
	Zf(FFT)(p, logn);
}

/* see inner.h */
void
Zf(int16_to_fft)(fpr *p, const int16_t *x, unsigned logn)
{
	size_t u;

	u = MKN(logn);
	while (u -- > 0) {
		p[u] = fpr_of(x[u]);
	}
	Zf(FFT)(p, logn);
}

/* see inner.h */
int
Zf(in_positive_half)(const int16_t *s, const uint8_t *h, unsigned logn)
{
	size_t n, u, v;
	uint32_t value, flag, set, result;
	uint8_t hash;

	n = MKN(logn);
	/*
	 * flag:   indicates whether a nonzero entry is seen yet,
	 * result: indicates whether the error vector e = h - 2s is valid so far,
	 * set:    indicates if the current value is nonzero.
	 */
	flag = 0;
	result = 1;

	if (logn <= 3) {
		for (u = 0; u < n; u ++) {
			value = (uint32_t)((h[0] >> u) & 1) - 2u * (uint32_t) s[u];
			set = (value | -value) >> 31;
			result &= flag | (set ^ 1u) | (-value >> 31);
			flag |= set;
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				value = (uint32_t)(hash & 1u) - 2u * (uint32_t) s[u];
				/*
				 * set = (value != 0);
				 * if (flag == 0 && set == 1) {
				 *     result = (value > 0);
				 *     flag = 1;
				 * }
				 */
				set = (value | -value) >> 31;
				/*
				 * result = 0 is only possible when flag = 0, set = 1 and
				 * -value >= 0, which is equivalent to flag = 0 and value < 0.
				 */
				result &= flag | (set ^ 1u) | (-value >> 31);
				flag |= set;
				hash >>= 1;
			}
		}
	}
	/*
	 * Both flag must be set and result must still be equal to 1.
	 */
	return (int) flag & result;
}
