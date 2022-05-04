/*

f,  g  <-- D_{Z^n, pk=1.425}
x0, x1 <-- D_{Z^n + h/2, sig=1.292}

|f|^2, |g|^2 ~ n sig_pk^2.

FFT seems to scale all standard deviations up by a factor 16 for n = 512 = 2^9...
*/

#include <climits>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <atomic>
#include <thread>
#include <vector>

// x86_64 specific:
#include<sys/time.h>

extern "C" {
	#define restrict
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

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

// =============================================================================
// | DISCRETE GAUSSIAN DISTRIBUTION                                            |
// =============================================================================

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

static int
mkgauss_1292(prng *rng, uint8_t double_mu)
{
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
	return *(int32_t *)&v;
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

static int
mkgauss_1425(prng *rng)
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
	return *(int32_t *)&v;
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
constexpr size_t logn = 9, n = MKN(logn);

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
		s = mkgauss_1425(rng);
		mod2 ^= (unsigned)(s & 1);
		f[u] = (int8_t)s;
	}

	do {
		s = mkgauss_1425(rng);
	} while (mod2 == (unsigned)(s & 1));
	f[0] = (int8_t)s;
}

static void
poly_sign_mkgauss(prng *rng, int8_t *x, uint8_t *h)
{
	for (size_t u = 0; u < n; u ++) {
		// TODO: can be improved, but not bottleneck probably
		uint8_t h_bit = (h[u/8] >> (u % 8)) & 1;
		x[u] = 2 * mkgauss_1292(rng, h_bit) - h_bit;
	}
}

const int num_samples = 100'000;

int simple_rounding_test(unsigned char *seed) {
	uint8_t h[n / 4];
	int8_t f[n], g[n], x0[n], x1[n];
	fpr fp[n], gp[n], x0p[n], x1p[n], q00[n], res[n];

	// Initialize a RNG.
	inner_shake256_context sc;
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);
	prng rng;
	Zf(prng_init)(&rng, &sc);

	int num_failed = 0;
	for (int rep = 0; rep < num_samples; rep++) {
		poly_small_mkgauss(&rng, f);
		poly_small_mkgauss(&rng, g);

		Zf(int8_to_fft)(fp, f);
		Zf(int8_to_fft)(gp, g);

		Zf(poly_invnorm2_fft)(q00, fp, gp, logn);

		Zf(prng_get_bytes)(&rng, (void *)h, sizeof h);

		poly_sign_mkgauss(&rng, x0, h);
		poly_sign_mkgauss(&rng, x1, h + n / 8);

		Zf(int8_to_fft)(x0p, x0);
		Zf(int8_to_fft)(x1p, x1);

		Zf(poly_add_muladj_fft)(res, x0p, x1p, fp, gp, logn); // f^* x0 + g^* x1
		Zf(poly_mul_fft)(res, q00, logn);

		Zf(iFFT)(res, logn);

		for (size_t u = 0; u < n; u++) {
			if (res[u].v <= -0.5 || res[u].v >= 0.5) [[unlikely]] {
				num_failed++;
				break;
			}
		}
	}

	return num_failed;

	printf("# simulations = %8d\n", num_samples);
	printf("#      failed = %8d\n", num_failed);
}

std::atomic<long long> total_fails = 0;
std::vector<unsigned char *> seeds;

void work(int id)
{
	int inc = simple_rounding_test(seeds[id]);
	total_fails += (long long)inc;
}

int main() {
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	const int nthreads = 4;
	for (int i = 0; i < nthreads; i++) {
		unsigned char *new_seed = (unsigned char *)malloc(48);
		randombytes(new_seed, 48);
		seeds.push_back(new_seed);
	}
	std::thread* pool[nthreads-1];
	for (int i = 0; i < nthreads-1; i++) pool[i] = new std::thread(work, 1 + i);
	work(0);
	for (int i = 0; i < nthreads-1; i++) pool[i]->join(), delete pool[i];
	for (int i = 0; i < nthreads; i++) {
		free(seeds[i]);
	}

	printf("# simulations = %9d\n", nthreads * num_samples);
	printf("# fails       = %9lld\n", (long long) total_fails);

	double fail_prob = double(total_fails) / nthreads / num_samples;
	printf("Fail prob.    = %.9f\n", fail_prob);

	return 0;
}
