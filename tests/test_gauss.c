/*
 * Prints average and std.dev. of the used discrete gaussian in HAWK for secret
 * key generation, i.e. with sigma_pk = 1.425 and prints estimations for the probabilities of certain outputs.
 *
 * Moreover, this prints a table to be used in keygen.c
 */

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "shake.c"

/*
 * Get a random 8-byte integer from a SHAKE-based RNG. This function
 * ensures consistent interpretation of the SHAKE output so that
 * the same values will be obtained over different platforms, in case
 * a known seed is used.
 */
static inline uint64_t
get_rng_u64(inner_shake256_context *rng)
{
	/*
	 * We enforce little-endian representation.
	 */
	uint8_t tmp[8];

	inner_shake256_extract(rng, tmp, sizeof tmp);
	return (uint64_t)tmp[0]
		| ((uint64_t)tmp[1] << 8)
		| ((uint64_t)tmp[2] << 16)
		| ((uint64_t)tmp[3] << 24)
		| ((uint64_t)tmp[4] << 32)
		| ((uint64_t)tmp[5] << 40)
		| ((uint64_t)tmp[6] << 48)
		| ((uint64_t)tmp[7] << 56);
}


/*
 * Table below incarnates a discrete Gaussian distribution:
 *    D(x) = exp(-(x^2)/(2*sigma^2))
 * where sigma = 1.425.
 * Element 0 of the table is P(x = 0).
 * For k > 0, element k is P(x >= k+1 | x > 0).
 * Probabilities are scaled up by 2^63.
 */
static const uint64_t gauss_1425[14] = {
	2582170577806070936u,
	3616484622030002669u,
	 937850763665829726u,
	 155804309064628172u,
	  16270507104385775u,
	   1056136771359479u,
	     42327595817352u,
	      1043181220683u,
	        15771375580u,
	          146052920u,
	             827733u,
	               2869u,
	                  6u,
	                  0u
};

/*
 * Generate a random value with a Gaussian distribution centered on 0.
 * The RNG must be ready for extraction (already flipped).
 *
 * Distribution has standard deviation 1.425. The code is now only usable for N
 * = 512 as other values for N would require different values for sigma.
 */
static int
mkgauss(inner_shake256_context *rng, unsigned logn)
{
	unsigned u, g;
	int val;

	g = 1U << (9 - logn);
	val = 0;
	for (u = 0; u < g; u ++) {
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
		r = get_rng_u64(rng);
		neg = (uint32_t)(r >> 63);
		r &= ~((uint64_t)1 << 63);
		f = (uint32_t)((r - gauss_1425[0]) >> 63);

		/*
		 * We produce a new random 63-bit integer r, and go over
		 * the array, starting at index 1. We store in v the
		 * index of the first array element which is not greater
		 * than r, unless the flag f was already 1.
		 */
		v = 0;
		r = get_rng_u64(rng);
		r &= ~((uint64_t)1 << 63);
		for (k = 1; k < 14; k ++) {
			uint32_t t;
			t = (uint32_t)((r - gauss_1425[k]) >> 63) ^ 1;
			v |= k & -(t & (f ^ 1));
			f |= t;
		}

		/*
		 * We apply the sign ('neg' flag). If the value is zero,
		 * the sign has no effect.
		 */
		v = (v ^ -neg) + neg;

		/*
		 * Generated value is added to val.
		 */
		val += *(int32_t *)&v;
	}
	return val;
}

void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

typedef long double FT;

int main() {
	srand(time(NULL));

	FT dx[100], sum = 0.0, p63 = powl(2.0, 63);
	for (int i = 0; i < 100; i++) {
		dx[i] = expl( - i*i / (FT)(2.0 * 1.425 * 1.425));
		sum += dx[i];
	}
	FT sum_geq[101];
	sum_geq[100] = 0.0;
	for (int i = 99; i >= 0; i--)
		sum_geq[i] = dx[i] + sum_geq[i + 1];

	FT P0 = dx[0] / (2.0*sum - dx[0]);
	printf("\t%19lluu,\n", (long long)(p63 * P0));
	for (int i = 1; i < 20; i++) {
		FT Pi = sum_geq[i + 1] / sum_geq[1];
		printf("\t%19lluu,\n", (long long)roundl(p63 * Pi));
	}


	unsigned char seed[48];
	inner_shake256_context sc;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	int freq[40] = {};
	long long sumval = 0, sumsqval = 0;	
	int num_sims = 10 * 1000 * 1000;
	for (int i = 0; i < num_sims; i++) {
		if (i > 1 && (i & (i-1)) == 0) {
			double avg = (double)sumval / i;
			double std = sqrt((double)sumsqval / i - avg*avg);
			printf("Std.dev: %.8f\n", std);
		}
		int x = mkgauss(&sc, 9);
		freq[x + 20]++;
		sumval += x;
		sumsqval += x*x;
	}

	double avg = (double)sumval / num_sims;
	double std = sqrt((double)sumsqval / num_sims - avg*avg);
	assert(fabs(avg) < 0.01);
	printf("Average: %.8f\n", avg);
	printf("Std.dev: %.8f\n", std);

	for (int x = -20; x < 20; x++) {
		FT Px = ((FT)freq[x + 20]) / num_sims;
		printf("P(X = %3d) = %.8f\n", x, (double) Px);
	}

	return 0;
}
