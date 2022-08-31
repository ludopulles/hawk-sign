#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../rng.c"
#include "../shake.c"

// =============================================================================

static const uint64_t gauss[13] = {
		6770309987939008324u, 2841792919453817158u, 824825004081786282u,
		160853309707784581u, 20707417942076380u, 1740733985516594u,
		94912702842187u, 3342501151111u, 75826385177u, 1106214542u, 10367241u,
		62372u, 240u
};

static int8_t secure_keygen(prng *rng)
{
	uint64_t r;
	uint8_t v, k, neg;

	/*
	 * Generate a 64 bit value.
	 * The highest bit determines the sign, say the other 63 bits give a value
	 * r. Then v = sum_{k=0}^{12} [[ v < gauss[k] ]].
	 */
	r = prng_get_u64(rng);

	/*
	 * Get the sign bit, and unset this bit in r.
	 */
	neg = (uint8_t)(r >> 63);
	r &= ~((uint64_t)1u << 63);

	v = 0;
	for (k = 0; k < 13; k++) {
		/*
		 * Add 1 iff r < gauss[k].
		 */
		v += (uint8_t)((uint64_t)(r - gauss[k]) >> 63);
	}

	/*
	 * Apply the sign ('neg' flag). If neg = 0, this has no effect.
	 * However, if neg = 1, this changes v into -v = (~v) + 1.
	 */
	v = (v ^ -neg) + neg;
	return *(int8_t *)&v;
}


// =============================================================================

static const uint16_t gauss_hi[10] = {
	0x580B, 0x35F9,
	0x1D34, 0x0DD7,
	0x05B7, 0x020C,
	0x00A2, 0x002B,
	0x000A, 0x0001,
};

static const uint64_t gauss_lo[26] = {
	0x0C27920A04F8F267, 0x3C689D9213449DC9,
	0x1C4FF17C204AA058, 0x7B908C81FCE3524F,
	0x5E63263BE0098FFD, 0x4EBEFD8FF4F07378,
	0x56AEDFB0876A3BD8, 0x4628BC6B23887196,
	0x061E21D588CC61CC, 0x7F769211F07B326F,
	0x2BA568D92EEC18E7, 0x0668F461693DFF8F,
	0x00CF0F8687D3B009, 0x001670DB65964485,
	0x000216A0C344EB45, 0x00002AB6E11C2552,
	0x000002EDF0B98A84, 0x0000002C253C7E81,
	0x000000023AF3B2E7, 0x0000000018C14ABF,
	0x0000000000EBCC6A, 0x000000000007876E,
	0x00000000000034CF, 0x000000000000013D,
	0x0000000000000006, 0x0000000000000000,
};

static inline int8_t
secure_sign(prng *rng, uint8_t parity)
{
	uint16_t r_hi, p_hi;
	uint64_t r_lo, p_lo;
	uint8_t c, v, k, neg;

	/*
	 * We use 80 random bits to determine the value, by looking in the
	 * cumulative probability table. However, we only use 15 bits for r_hi so
	 * we can check if r_hi < p_hi holds by computing (r_hi - p_hi) >> 15. This
	 * way is better than doing a comparison to achieve constant-time
	 * execution.
	 */
	prng_get_80_bits(rng, &r_hi, &r_lo);

	/*
	 * Get the sign bit out of the lowest part
	 */
	neg = (uint8_t)(r_lo >> 63);

	/*
	 * Unset the sign bits in the unsigned ints for convenience in comparisons
	 * later on, as we can now use the highest bit of `a - b` to check if `a <
	 * b` or not for numbers `a, b`.
	 */
	r_hi &= ~((uint16_t)1u << 15);
	r_lo &= ~((uint64_t)1u << 63);

	v = 0;
	for (k = 10; k < 26; k += 2) {
		/*
		 * Constant-time for:
		 *     p_lo = gauss_lo[k + parity];
		 */
		p_lo = (gauss_lo[k] & (parity - 1u)) | ((gauss_lo[k + 1] & -parity));

		/*
		 * Add 1 iff r_lo < p_lo.
		 */
		v += (uint8_t)((uint64_t)(r_lo - p_lo) >> 63);
	}

	/*
	 * If r_hi > 0, set v to zero, otherwise leave v as is. This is a
	 * micro-optimization as p_hi would be zero for all k >= 5.
	 */
	v = v & -((r_hi - 1) >> 15);

	for (k = 0; k < 10; k += 2) {
		/*
		 * Constant-time for:
		 *     p_lo = gauss_lo[k + parity];
		 *     p_hi = gauss_hi[k + parity];
		 */
		p_lo = (gauss_lo[k] & (parity - 1u)) | ((gauss_lo[k + 1] & -parity));
		p_hi = (gauss_hi[k] & (parity - 1u)) | ((gauss_hi[k + 1] & -parity));

		/*
		 * c = [[ r_lo < p_lo ]]
		 */
		c = (uint8_t)((uint64_t)(r_lo - p_lo) >> 63);

		/*
		 * Constant-time code to add 1 to v iff
		 *     r_hi < p_hi or (r_hi == p_hi and c is true)
		 * holds.
		 */
		c = (uint8_t)((uint16_t)(r_hi - p_hi - c) >> 15);
		v += c;
	}

	/*
	 * Multiply by two and apply the change in support:
	 * If parity = 0, then v = 0,2,4,...
	 * If parity = 1, then v = 1,3,5,...
	 */
	v = (v << 1) | parity;

	/*
	 * Apply the sign ('neg' flag). If neg = 0, this has no effect.
	 * However, if neg = 1, this changes v into -v = (~v) + 1.
	 */
	v = (v ^ -neg) + neg;
	return *(int8_t *)&v;

}

// =============================================================================

#define NUM_SAMPLES ((int64_t)(10 * 1000000))
typedef double FT;

/*
 * As we are sampling, our epsilon must be taken quite large, since sampling
 * does not give the exact correct probability, but only close up to a factor
 * of ~ 1 / sqrt(NUM_SAMPLES).
 */
const FT sigma_pk = 1.500, sigma_sig = 1.278;

FT P[100] = {};
int freq[100] = {};

FT rho(FT x, FT sigma) {
	return expl(-x * x / (2.0 * sigma * sigma));
}

int main() {
	FT eps = 2.0 / sqrt(NUM_SAMPLES);
	inner_shake256_context sc;
	unsigned char seed[48];
	prng p;
	FT norm_pk = 0.0, norm[2] = { 0.0, 0.0 };

	// Initialize PRNG
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);
	Zf(prng_init)(&p, &sc);


	// Key generation sampler:
	for (int i = -50; i < 50; i++)
		P[i + 50] = rho(i, sigma_pk), norm_pk += P[i + 50];
	// normalize:
	for (int i = -50; i < 50; i++) P[i + 50] /= norm_pk;

	memset(freq, 0, sizeof freq);
	for (int i = NUM_SAMPLES; i--; ) {
		int x = secure_keygen(&p);
		freq[x + 50]++;
	}

	FT max_abs_diff = 0.0, max_rel_diff = 0.0;
	for (int i = -50; i < 50; i++) {
		if (freq[i + 50] == 0) continue;
		FT found = (FT) freq[i + 50] / NUM_SAMPLES, expected = P[i + 50];
		// printf("P(X == %d) = %12.8f, %12.8f\n", i, found, expected);

		FT abs_diff = fabs(found - expected);
		FT rel_diff = abs_diff / expected;
		if (abs_diff > max_abs_diff)
			max_abs_diff = abs_diff;
		if (rel_diff < abs_diff && rel_diff > max_rel_diff)
			max_rel_diff = rel_diff;

		assert(abs_diff < eps || rel_diff < eps);
	}
	printf("Keygen sampler        has abs.diff. < %.5f and rel.diff. < %.5f\n",
		max_abs_diff, max_rel_diff);


	// Signing sampler:
	for (int i = -50; i < 50; i++) {
		// sample from 2Z + coset by scaling everything up by a factor of two
		P[i + 50] = rho(i, 2.0 * sigma_sig);
		norm[((unsigned) i) & 1] += P[i + 50];
	}
	// normalize:
	for (int i = -50; i < 50; i++)
		P[i + 50] /= norm[((unsigned) i) & 1];

	for (int8_t coset = 0; coset < 2; coset++) {
		memset(freq, 0, sizeof freq);

		for (int i = NUM_SAMPLES; i--; ) {
			int x = secure_sign(&p, coset);
			freq[x + 50]++;
		}

		max_abs_diff = max_rel_diff = 0.0;
		for (int i = -50; i < 50; i++) {
			if (freq[i + 50] == 0) continue;
			FT found = (FT) freq[i + 50] / NUM_SAMPLES, expected = P[i + 50];
			// printf("P(X == %d) = %12.8f, %12.8f\n", i, found, expected);

			FT abs_diff = fabs(found - expected);
			FT rel_diff = abs_diff / expected;
			if (abs_diff > max_abs_diff)
				max_abs_diff = abs_diff;
			if (rel_diff < abs_diff && rel_diff > max_rel_diff)
				max_rel_diff = rel_diff;

			assert(abs_diff < eps || rel_diff < eps);
		}
		printf("Sign sampler (2Z + %d) has abs.diff. < %.5f and rel.diff. < %.5f\n",
			coset, max_abs_diff, max_rel_diff);
	}

	return 0;
}
