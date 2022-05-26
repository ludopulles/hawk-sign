/*
 * Common functions for Hawk
 */

#include "inner.h"

/* see inner.h */
const uint32_t Zf(l2bound)[10] = {
	0u /* unused */, 32u, 64u, 129u, 259u, 519u, 1039u, 2079u, 4158u, 8317u
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
Zf(fft_to_int16)(int16_t *x, fpr *p, unsigned logn)
{
	size_t u;

	u = MKN(logn);

	Zf(iFFT)(p, logn);
	while (u -- > 0) {
		x[u] = fpr_rint(p[u]);
	}
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
