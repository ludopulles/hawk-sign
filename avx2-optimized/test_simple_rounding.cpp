/*

f,  g  <-- D_{Z^n, pk=1.425}
x0, x1 <-- D_{Z^n + h/2, sig=1.292}

|f|^2, |g|^2 ~ n sig_pk^2.

FFT seems to scale all standard deviations up by a factor 16 for n = 512 = 2^9...
*/

#include <bits/stdc++.h>

// x86_64 specific:
#include<sys/time.h>

extern "C" {
	#ifndef restrict
		#define restrict
	#endif

	#include "inner.h"
	#include "fpr.c"
	#include "fft.c"
	#include "rng.c"
	#include "shake.c"
}

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

void random_hash(int8_t *h, unsigned logn) {
	assert(RAND_MAX == INT_MAX); // rand() should generate 31 random bits
	int x = rand();
	size_t RAND_BITS = 31, rand_bits = RAND_BITS;
	for (size_t u = MKN(logn); u -- > 0; ) {
		if (rand_bits == 0) {
			x = rand();
			rand_bits = RAND_BITS;
		}
		h[u] = (x & 1);
		x >>= 1;
		rand_bits--;
	}
}

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

// =============================================================================
// | DISCRETE GAUSSIAN DISTRIBUTION                                            |
// =============================================================================

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
static const uint64_t gauss_1292[24] = {
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
//	                  0u,                   0u,
};

/*
 * Generate a random value with a Gaussian distribution centered on double_mu/2.
 * The RNG must be ready for extraction (already flipped).
 *
 * Distribution has standard deviation 1.292 sqrt(512/N).
 * Note/TODO: Watch out with N != 512, the Renyi-divergence with the perfect
 * discrete gaussian is unknown.
 */
static int
mkgauss_1292(prng *rng, unsigned logn, uint8_t double_mu)
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
		 * We use two random 64-bit values. First value decides on whether the
		 * generated value is 0, and, if not, the sign of the value. Second
		 * random 64-bit word is used to generate the non-zero value.
		 *
		 * For constant-time code we have to read the complete table. Currently
		 * sampling takes up most time of the signing process, so this is
		 * something that should be optimized.
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
		f = (uint32_t)((r - gauss_1292[double_mu]) >> 63);

		/*
		 * We produce a new random 63-bit integer r, and go over
		 * the array, starting at index 1. We store in v the
		 * index of the first array element which is not greater
		 * than r, unless the flag f was already 1.
		 */
		v = 0;

		r = prng_get_u64(rng);

		r &= ~((uint64_t)1 << 63);
		for (k = 1; k < 12; k ++) {
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

		/*
		 * In the case that this code is run multiple times, we want to use
		 * center mu = 0 in the other runs (g > 0), as to not change the center
		 * and only scale the standard deviation. However, we might be a bit off
		 * since sum of discrete gaussians might not be a discrete gaussian (it
		 * is only true in the continuous limit).
		 */
		double_mu = 0;
	}
	return val;
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
 * Distribution has standard deviation 1.425 sqrt(512/N).
 */
static int
mkgauss_1425(prng *rng, unsigned logn)
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
		r = prng_get_u64(rng);
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
		r = prng_get_u64(rng);
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

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
constexpr size_t logn = 9, n = MKN(logn);

/*
 * Convert an integer polynomial (with small values) into the
 * representation with complex numbers.
 * Also works when r and t overlap, since the loop goes in decreasing order.
 */
static void
smallints_to_fpr(fpr *r, const int8_t *t) {
	for (size_t u = n; u --> 0; ) r[u] = fpr_of(t[u]);
}

/*
 * Generate a random polynomial with a Gaussian distribution. This function
 * also makes sure that the resultant of the polynomial with phi is odd.
 */
static void
poly_small_mkgauss(prng *rng, int8_t *f)
{
	size_t u;
	int s;
	unsigned mod2;

	mod2 = 0;
	for (u = n; u -- > 1; ) {
		s = mkgauss_1425(rng, logn);
		mod2 ^= (unsigned)(s & 1);
		f[u] = (int8_t)s;
	}

	do {
		s = mkgauss_1425(rng, logn);
	} while (mod2 == (unsigned)(s & 1));
	f[0] = (int8_t)s;
}

static void
poly_sign_mkgauss(prng *rng, int8_t *x, uint8_t *h)
{
	for (size_t u = 0; u < n; u ++) {
		uint8_t h_bit = (h[u/8] >> (u % 8)) & 1;
		x[u] = 2 * mkgauss_1292(rng, logn, h_bit) - h_bit;
	}
}

void simple_rounding_test() {
	uint8_t h[n / 4];
	int8_t f[n], g[n], x0[n], x1[n];
	fpr fp[n], gp[n], x0p[n], x1p[n], q00[n], res[n];

	// Initialize a RNG.
	unsigned char seed[48];
	inner_shake256_context sc;
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);
	prng rng;
	Zf(prng_init)(&rng, &sc);

	double avg[2*n] = {}, var[2*n] = {};
	double avgQ[n] = {}, varQ[n] = {};
	double avgres[n] = {}, varres[n] = {};

#define ADDSAMPLE(poly, A, V) { \
	for (size_t u = 0; u < n; u++) { \
		A[u] += (poly)[u].v; \
		V[u] += (poly)[u].v * (poly)[u].v; \
	} \
}

	const int num_samples = 100000;
	for (int rep = 0; rep < num_samples; rep++) {
		poly_small_mkgauss(&rng, f);
		poly_small_mkgauss(&rng, g);

		randombytes(h, sizeof h);

		poly_sign_mkgauss(&rng, x0, h);
		poly_sign_mkgauss(&rng, x1, h + n/8);
		smallints_to_fpr(fp, f);
		smallints_to_fpr(gp, g);
		smallints_to_fpr(x0p, x0);
		smallints_to_fpr(x1p, x1);
		Zf(FFT)(fp, logn);
		Zf(FFT)(gp, logn);
		Zf(FFT)(x0p, logn);
		Zf(FFT)(x1p, logn);
		Zf(poly_add_muladj_fft)(q00, fp, gp, fp, gp, logn);

		// Now sample f,g again
		/* 
		poly_small_mkgauss(&rng, f);
		poly_small_mkgauss(&rng, g);
		smallints_to_fpr(fp, f);
		smallints_to_fpr(gp, g);
		Zf(FFT)(fp, logn);
		Zf(FFT)(gp, logn); */

		// ADDSAMPLE(q00, avgQ, varQ);
		// Q00 is self-adjoint, so only take real parts
		for (size_t u = 0; u < n/2; u++)
			avgQ[u] += q00[u].v, varQ[u] += fpr_sqr(q00[u]).v;

		Zf(poly_div_autoadj_fft)(fp, q00, logn);
		Zf(poly_div_autoadj_fft)(gp, q00, logn);

		Zf(poly_add_muladj_fft)(res, x0p, x1p, fp, gp, logn); // f^* x0 + g^* x1

		Zf(iFFT)(fp, logn);
		Zf(iFFT)(gp, logn);
		Zf(iFFT)(res, logn);

		ADDSAMPLE(fp, avg, var);
		ADDSAMPLE(gp, (avg+n), (var+n));
		ADDSAMPLE(res, avgres, varres);
	}

	double tavg = 0.0, tvar = 0.0;
	for (size_t u = 0; u < n+n; u++) {
		tavg += avg[u];
		tvar += var[u];
	}

	double xavg = tavg / (2*n*num_samples);
	double xstd = sqrt(tvar / (2*n*num_samples) - xavg * xavg);
	printf("normalized f mu,sigma = %.10f +/- %.10f\n", xavg, xstd);

	tavg = tvar = 0.0;
	for (size_t u = 0; u < n; u++) {
		tavg += avgres[u];
		tvar += varres[u];
	}
	xavg = tavg / (n*num_samples);
	xstd = sqrt(tvar / (n*num_samples) - xavg * xavg);
	printf("res mu,sigma = %.10f +/- %.10f\n", xavg, xstd);

	tavg = tvar = 0.0;
	for (size_t u = 0; u < n/2; u++) {
		tavg += avgQ[u], tvar += varQ[u];
	}

	xavg = tavg * 2 / (n*num_samples);
	xstd = sqrt(tvar * 2 / (n*num_samples) - xavg * xavg);
	printf("q00 mu,sigma = %.10f +/- %.10f\n", xavg, xstd);
}

int main() {
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	simple_rounding_test();
	return 0;
}
