/*
 * Common functions for Hawk
 */

#include "inner.h"

/*
 * This number indicates the maximum l2-norm that is allowed for the small
 * error (from a valid signature) that is chosen around the target point during
 * signature generation.
 * The coefficients of this error are distributed according to a discrete
 * gaussian over the integers that have the same parity as the target
 * coefficient and with standard deviation 2 sigma_sig. To compute these
 * values, use:
 *
 *     l2bound(logn) = floor( (verif_margin * 2 sigma_sig)^2 * 2n ).
 *
 * Here, we have taken verif_margin = 1.1 and sigma_sig = 1.292.
 */
const uint32_t Zf(l2bound)[10] = {
	0 /* unused */, 32, 64, 129, 258, 517, 1034, 2068, 4136, 8273
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
