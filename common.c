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
