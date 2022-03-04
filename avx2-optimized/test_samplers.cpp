#include<bits/stdc++.h>

#define restrict
#include "rng.c"
#include "shake.c"

void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}


/** To generate the values in the table below, run the following code:

#include<bits/stdc++.h>
int main() {
	long double sigma = 1.292, mu = 0, p63 = powl(2, 63), table[100], csum[100];
	unsigned long long results[2][20] = {};
	for (int i = 0; i < 2; i++, mu += 0.5) {
		for (int x = 100; x --> 0; ) {
			table[x] = expl(-0.5 * (x-mu)*(x-mu) / sigma / sigma);
			csum[x] = table[x];
			if (x < 99) csum[x] += csum[x+1];
		}
		results[i][0] = llroundl((1.0+i) * p63 * table[i] / (csum[i] + csum[1]));
		for (int x = 19; --x >= 1; )
			results[i][x] = llroundl(p63 * csum[1+i+x] / csum[1+i]);
	}
	for (int x = 0; x < 20; x++)
		printf("\t%19lluu, %19lluu\n", results[0][x], results[1][x]);
}

 */

/*
 * Table below incarnates two discrete Gaussian distribution:
 *    D(x) = exp(-((x - mu)^2)/(2*sigma^2))
 * where sigma = 1.292 and mu is 0 or 1/2.
 * Element 0 of the first table is P(x = 0) and 2*P(x = 1) in the second table.
 * For k > 0, element k is P(x >= k+1 | x > 0) in the first table, and
 * P(x >= k+2 | x > 1) in the second table.
 * For constant-time principle, mu = 0 is in the even indices and
 * mu = 1/2 is in the odd indices.
 * Probabilities are scaled up by 2^63.
 */
static const uint64_t gauss_1292[26] = {
	2847982254933138603u, 5285010687306232178u,
	3115855658194614154u, 2424313226695581870u,
	 629245045388085487u,  372648834165936922u,
	  73110737927091842u,   32559817584178793u,
	   4785625785139312u,    1592210133688742u,
	    174470148146634u,      43209976786070u,
	      3520594834759u,        647780323462u,
	        39186846585u,          5350987999u,
	          240149359u,            24322099u,
	             809457u,               60785u,
	               1500u,                  83u,
	                  2u,                   0u,
	                  0u,                   0u,
};

/*
 * Generate a random value with a Gaussian distribution centered on double_mu/2.
 * The RNG must be ready for extraction (already flipped).
 *
 * Distribution has standard deviation 1.292 sqrt(512/N).
 */
static int
mkgauss(void *samp_ctx, unsigned logn, uint8_t double_mu)
{
	unsigned u, g;
	int val;

	sampler_context *sc = (sampler_context *)samp_ctx;

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
		r = prng_get_u64(&sc->p);
		neg = (uint32_t)(r >> 63);
		r &= ~((uint64_t)1 << 63);
		f = (uint32_t)((r - gauss_1292[double_mu]) >> 63);

		/*
		 * We produce a new random 63-bit integer r, and go over
		 * the array, starting at index 1. We store in v the
		 * index of the first array element which is not greater
		 * than r, unless the flag f was already 1.
		 */
		v = 0;
		r = prng_get_u64(&sc->p);
		r &= ~((uint64_t)1 << 63);
		for (k = 1; k < 13; k ++) {
			uint32_t t;
			t = (uint32_t)((r - gauss_1292[2 * k + double_mu]) >> 63) ^ 1;
			v |= k & -(t & (f ^ 1));
			f |= t;
		}

		/*
		 * We apply the sign ('neg' flag). If the value is zero and mu = 0,
		 * the sign has no effect. Moreover, if mu = 1/2 and neg=0, add one.
		 */
		v = (v ^ -neg) + neg + (~neg & double_mu);

		/*
		 * Generated value is added to val.
		 */
		val += *(int32_t *)&v;
	}
	return val;
}

int main() {
	srand(time(NULL));
	unsigned char seed[48];
	inner_shake256_context sc;
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	sampler_context spc;
	Zf(prng_init)(&spc.p, &sc);
	void *samp_ctx = &spc;

	int sumv = 0, sqv = 0;
	int nums = 1 << 25;
	int freq[100] = {};
	for (int i = nums; i-->0; ) {
		int x = mkgauss(samp_ctx, 9, 1);
		sumv += x;
		sqv += x*x;
		freq[x + 50]++;
	}
	printf("Sum, square: %d %d\n", sumv, sqv);
	double avg = ((double)sumv) / nums;
	double std = sqrt(((double)sqv) / nums - avg*avg);
	printf("%.10f +/- %.10f\n", avg, std);

	for (int x = -20; x < 20; x++)
		printf("%d: %d\n", x, freq[x + 50]);
}
