#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../rng.c"
#include "../shake.c"

// =============================================================================

/*
 * Table below incarnates a discrete Gaussian distribution:
 *    D(x) = exp(-(x^2)/(2*sigma^2))
 * where sigma = 1.500.
 * Element 0 of the table is P(x = 0).
 * For k > 0, element k is P(x >= k+1 | x > 0).
 * Probabilities are scaled up by 2^63.
 *
 * To generate the values in the table below, run 'gen_table.cpp'.
 */
static const uint64_t gauss_keygen[14] = {
	2453062048915767484u,
	3871449519226705105u,
	1123680878940444328u,
	 219134710439743982u,
	  28210262150869885u,
	   2371447864901096u,
	    129302080834770u,
	      4553577562215u,
	       103300286390u,
	         1507025277u,
	           14123567u,
	              84972u,
	                328u,
	                  1u
};

/*
 * Generate a random value with a Gaussian distribution centered on 0.
 * The RNG must be ready for extraction (already flipped).
 * Distribution has standard deviation 1.5.
 */
static int
mkgauss_keygen(prng *rng)
{
	/*
	 * Each iteration generates one value with the
	 * Gaussian distribution for N = 512.
	 *
	 * We use two random 64-bit values. First value
	 * decides on whether the generated value is 0, and,
	 * if not, the sign of the value. Second random 64-bit
	 * word is used to generate the non-zero value.
	 *
	 * For constant-time code we have to read the complete
	 * table. This has negligible cost, compared with the
	 * remainder of the keygen process (solving the NTRU
	 * equation).
	 */
	uint64_t r;
	uint32_t f, v, k, neg;

	/*
	 * First value:
	 *  - flag 'neg' is randomly selected to be 0 or 1.
	 *  - flag 'f' is set to 1 if the generated value is zero,
	 *    or set to 0 otherwise.
	 */
	r = prng_get_u64(rng);
	neg = (uint32_t)(r >> 63);
	r &= ~((uint64_t)1 << 63);
	f = (uint32_t)((r - gauss_keygen[0]) >> 63);

	/*
	 * We produce a new random 63-bit integer r, and go over
	 * the array, starting at index 1. We store in v the
	 * index of the first array element which is not greater
	 * than r, unless the flag f was already 1.
	 */
	v = 0;
	r = prng_get_u64(rng);
	r &= ~((uint64_t)1 << 63);
	for (k = 1; k < 14; k ++) {
		uint32_t t;
		t = (uint32_t)((r - gauss_keygen[k]) >> 63) ^ 1;
		v |= k & -(t & (f ^ 1));
		f |= t;
	}

	/*
	 * We apply the sign ('neg' flag). If the value is zero,
	 * the sign has no effect.
	 */
	v = (v ^ -neg) + neg;
	return *(int32_t *)&v;
}


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

/*
 * Table below incarnates two discrete Gaussian distribution:
 *    D(x) = exp(-((x - mu)^2)/(2*sigma^2))
 * where sigma = 1.278 and mu is 0 or 1 / 2.
 * Element 0 of the first table is P(x = 0) and 2*P(x = 1) in the second table.
 * For k > 0, element k is P(x >= k+1 | x > 0) in the first table, and
 * P(x >= k+2 | x > 1) in the second table.
 * For constant-time principle, mu = 0 is in the even indices and
 * mu = 1 / 2 is in the odd indices.
 * Probabilities are scaled up by 2^63.
 *
 * To generate the values in the table below, run 'gen_table.cpp'.
 */
static const uint64_t gauss_sign[24] = {
	2879180808586524119u, 5334099443607918879u,
	3059393892786389767u, 2365650607671913513u,
	 598985657464043109u,  350188603152862536u,
	  66570217148287807u,   29069029197456365u,
	   4111724430703280u,    1332230758935082u,
	    139536889245539u,      33427798003753u,
	      2585835445815u,        457141300403u,
	        26080447382u,          3398953827u,
	          142905880u,            13721995u,
	             424994u,               30058u,
	                686u,                  36u,
	                  1u,                   0u
};

/*
 * Sample an integer with parity equal to double_mu from a discrete Gaussian
 * distribution with support 2\ZZ + double_mu, mean 0 and sigma 2 * 1.278.
 * That is, an integer x (== double_mu mod 2) is chosen with probability
 * proportional to:
 *
 *     exp(- x^2 / (8 1.278^2)).
 */
static inline int
mkgauss_sign(prng *rng, uint8_t double_mu)
{
	uint64_t r;
	uint32_t v, k, neg;
	int32_t w;

	/*
	 * We use two random 64-bit values. First value is used to generate the
	 * non-zero value. Second value decides on whether the generated value is 0
	 * and the sign of the value. For constant-time code we have to read the
	 * complete table.
	 */

	r = prng_get_u64(rng) & ~((uint64_t)1 << 63);
	v = 1;
	for (k = 1; k < 12; k ++) {
		/*
		 * Add 1 if r < gauss_sign[2 * k + double_mu].
		 */
		v += (r - gauss_sign[2 * k + double_mu]) >> 63;
	}

	/*
	 * First value:
	 *  - flag 'neg' is randomly selected to be 0 or 1.
	 *  - if r/2^63 <= P(X == 0), then set v to zero.
	 */
	r = prng_get_u64(rng);
	neg = (uint32_t)(r >> 63);
	r &= ~((uint64_t)1 << 63);
	v &= -((gauss_sign[double_mu] - r) >> 63);

	/*
	 * We apply the sign ('neg' flag).
	 * If mu = 0 change v to v if neg = 0 and -v if neg = 1.
	 * If mu = 1/2, change v to v + 1 if neg = 0 and -v if neg = 1.
	 */
	v = (v ^ -neg) + neg + (~neg & double_mu);

	/*
	 * Now, transform the support of the sampler from Z to 2Z - double_mu, i.e.
	 * for mu = 1/2, we have -1, 1 with equal likelihood and -2, 2 with equal
	 * likelihood, etc.
	 */
	w = *(int32_t *)&v;
	return 2 * w - (int) double_mu;
}


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
const FT eps = 2.0 / sqrt(NUM_SAMPLES);

FT P[100] = {};
int freq[100] = {};

FT rho(FT x, FT sigma) {
	return expl(-x * x / (2.0 * sigma * sigma));
}

int main() {
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
