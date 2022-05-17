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

// =============================================================================

const int num_samples = 10 * 1000 * 1000;
typedef double FT;

/*
 * As we are sampling, our epsilon must be taken quite large, since sampling
 * does not give the exact correct probability, but only close up to a factor
 * of ~ 1 / sqrt(num_samples).
 */
const FT sigma_pk = 1.500, sigma_sig = 1.278;

FT eps = 0.001, P[100] = {};
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
	for (int i = num_samples; i--; ) {
		int x = mkgauss_keygen(&p);
		freq[x + 50]++;
	}

	FT max_abs_diff = 0.0, max_rel_diff = 0.0;
	for (int i = -50; i < 50; i++) {
		if (freq[i + 50] == 0) continue;
		FT found = (FT) freq[i + 50] / num_samples, expected = P[i + 50];
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

	for (int coset = 0; coset < 2; coset++) {
		memset(freq, 0, sizeof freq);

		for (int i = num_samples; i--; ) {
			int x = mkgauss_sign(&p, coset);
			freq[x + 50]++;
		}

		max_abs_diff = max_rel_diff = 0.0;
		for (int i = -50; i < 50; i++) {
			if (freq[i + 50] == 0) continue;
			FT found = (FT) freq[i + 50] / num_samples, expected = P[i + 50];
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
