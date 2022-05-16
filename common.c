/*
 * Common functions for Hawk
 */

#include "inner.h"

/* see inner.h */
const uint32_t Zf(l2bound)[10] = {
	0u /* unused */, 31u, 63u, 126u, 252u, 505u, 1011u, 2023u, 4047u, 8094u
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
