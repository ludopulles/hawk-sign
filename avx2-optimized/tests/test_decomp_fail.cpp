/*
 * Determine the probability that decompressing a signature s1 gives a
 * different (s0, s1) than the original s0.
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <thread>
#include <mutex>
#include <vector>

extern "C" {
	#define restrict
	#include "../inner.h"
}

// =============================================================================
// SIGNING
// =============================================================================
// Taken from sign.c:
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

#define SECOND_HASH(h, logn) \
	((h) + ((logn) <= 3 ? 1u : 1u << ((logn) - 3)))

static void
hash_to_fft(fpr *p, const uint8_t *h, unsigned logn)
{
	size_t n, u, v;
	uint8_t hash;

	n = MKN(logn);

	if (logn <= 3) {
		for (v = 0; v < n; v ++) {
			p[v] = fpr_of((h[0] >> v) & 1);
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				p[u] = fpr_of(hash & 1);
				hash >>= 1;
			}
		}
	}
	Zf(FFT)(p, logn);
}

/*
 * Returns the vector s = (h - noise) / 2 which gives a lattice point s that is
 * close to h / 2, assuming that noise has a small norm. This function assumes
 * that noise is in coefficient representation and integral (up to rounding
 * errors).
 */
static void
noise_to_lattice(int16_t *s, const uint8_t *h, const fpr *noise, unsigned logn)
{
	size_t n, u, v;
	uint8_t hash;
	int64_t x;

	n = MKN(logn);
	if (logn <= 3) {
		for (v = 0; v < n; v ++) {
			x = fpr_rint(noise[v]);
			s[v] = (int16_t)((int64_t)((h[0] >> v) & 1) - x) / 2;
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				x = fpr_rint(noise[u]);
				s[u] = (int16_t)((int64_t)(hash & 1) - x) / 2;
				hash >>= 1;
			}
		}
	}
}


static int
do_sign(prng *rng, const fpr *restrict expanded_seckey,
	int16_t *s0, int16_t *s1, const uint8_t *restrict h, unsigned logn,
	uint8_t *restrict tmp)
{
	size_t n, u;
	const fpr *bf, *bg, *bF, *bG, *invq00;
	fpr *x0, *x1, *res;
	int32_t z, norm;

	n = MKN(logn);
	norm = 0;

	bf = expanded_seckey;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	invq00 = bG + n;

	x0 = (fpr *)tmp;
	x1 = x0 + n;
	res = x1 + n;

	hash_to_fft(x0, h, logn);
	hash_to_fft(x1, SECOND_HASH(h, logn), logn);
	Zf(poly_matmul_fft)(bf, bF, bg, bG, x0, x1, logn);

	Zf(iFFT)(x0, logn);
	Zf(iFFT)(x1, logn);

	do {
		for (u = 0; u < n; u ++) {
			z = mkgauss_sign(rng, fpr_rint(x0[u]) & 1);
			norm += z*z;
			x0[u] = fpr_of(z);
		}
		for (u = 0; u < n; u ++) {
			z = mkgauss_sign(rng, fpr_rint(x1[u]) & 1);
			norm += z*z;
			x1[u] = fpr_of(z);
		}
	} while ((uint32_t)norm > Zf(l2bound)[logn]);

	Zf(FFT)(x0, logn);
	Zf(FFT)(x1, logn);

	Zf(poly_add_muladj_fft)(res, x0, x1, bf, bg, logn);
	Zf(poly_mul_autoadj_fft)(res, invq00, logn);
	Zf(iFFT)(res, logn);

	int num_errors = 0;
	for (u = 0; u < n; u++) {
		if (!fpr_lt(fpr_neg(fpr_one), res[u]) || !fpr_lt(res[u], fpr_one))
			num_errors++;
	}

	/* if (num_errors > 10) {
		printf("ERRORS!\n");
		for (u = 0; u < n; u++) {
			printf("%.2f ", *(double *)&res[u]);
		}
		printf("\n");
	} */

	/*
	 * Calculate B^{-1} (x0, x1) = ((-G) x0 + F x1, g x0 + (-f) x1).
	 */
	Zf(poly_neg)(x0, logn);
	Zf(poly_matmul_fft)(bG, bF, bg, bf, x0, x1, logn);
	Zf(poly_neg)(x1, logn);

	/*
	 * Extract the signature from x0, t1
	 */
	Zf(iFFT)(x0, logn);
	Zf(iFFT)(x1, logn);

	noise_to_lattice(s0, h, x0, logn);
	noise_to_lattice(s1, SECOND_HASH(h, logn), x1, logn);
	return 1 + num_errors;
}


#define LOGN (9)
#define N MKN(LOGN)
uint8_t b[48 << LOGN];
int8_t f[N], g[N], F[N], G[N];
fpr q00[N], q10[N], q11[N], iq00[N], exp_key[9 << (LOGN-1)];
constexpr size_t num_reps = 100 * 1000;

void measure_decompression_failure(inner_shake256_context *sc) {
	prng rng;
	uint8_t h[N/4], owntmp[24*N];
	int16_t s0[N], s0p[N], s1[N];

	inner_shake256_extract(sc, rng.state.d, 56);
	// initialize only once for performance
	Zf(prng_init)(&rng, sc);

	int res;
	for (size_t i = 0; i < num_reps; i++) {
		// Generate a random hash to sign.
		Zf(prng_get_bytes)(&rng, h, sizeof h);
		res = do_sign(&rng, exp_key, s0, s1, h, LOGN, owntmp);
		if (res > 1) printf("res = %d\n", res);
		if (res	> 1 && Zf(verify_simple_rounding)(h, s0p, s1, q00, q10, q11, LOGN, owntmp)) {
			printf("Expecting different signature!!!\n");
			for (size_t u = 0; u < N; u++)
				printf("%d ", s0p[u] - s0[u]);
			printf("\n");
			fflush(stdout);
			assert(!Zf(complete_verify)(h, s0, s1, q00, q10, q11, LOGN, b));
		}
	}
}

typedef long double FT;
constexpr FT sigma_sig = 1.278;
FT erfcl_s(FT x, FT sigma) { return erfcl(x / (sigma * sqrtl(2))); }
FT fail_prob(FT ciq00, long long n) {
	FT val = 1.0 - powl(erfl(0.5 / sqrtl(2*ciq00) / sigma_sig), n); /* direct approach */
	if (val < 1e-50) return n * erfcl_s(0.5, sqrt(ciq00) * sigma_sig); /* union bound */
	return val;
}

void spawn_threads(inner_shake256_context *sc, int num_threads)
{
	unsigned char seed[48];
	std::vector<std::thread*> pool;
	std::vector<inner_shake256_context*> shakes;

	for (int i = 0; i < num_threads; i++) {
		inner_shake256_context *newsc = new inner_shake256_context();
		inner_shake256_init(newsc);
		inner_shake256_inject(newsc, seed, sizeof seed);
		inner_shake256_flip(newsc);
		shakes.push_back(newsc);
	}

	for (int i = 0; i < num_threads - 1; i++) {
		pool.push_back(new std::thread(measure_decompression_failure, shakes[i]));
	}
	measure_decompression_failure(shakes[num_threads - 1]);
	for (int i = 0; i < num_threads - 1; i++) {
		pool[i]->join();
		delete pool[i];
	}
	for (int i = 0; i < num_threads; i++)
		delete shakes[i];
}

void test_decomp_prob(inner_shake256_context *sc)
{
	size_t u, hn = N/2;

	Zf(keygen)(sc, f, g, F, G, q00, q10, q11, LOGN, b);

	for (u = 0; u < hn; u ++) {
		iq00[u] = fpr_inv(q00[u]);
		iq00[u + hn] = fpr_zero;
	}
	Zf(iFFT)(iq00, LOGN);

	double cst_iq = *((double*)&iq00[0]);
	FT p = fail_prob(cst_iq, N);

	printf("Probability cst(1/Q00) = %.5f: %Le exp\n", cst_iq, p);

	Zf(expand_seckey)(exp_key, f, g, F, LOGN);

	spawn_threads(sc, 4);
}

int main() {
	unsigned char seed[48] = "seedseedseedseedseedseedseedseedseedseedseedsee";
	inner_shake256_context sc;

	// Initialize a RNG.
	// Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	test_decomp_prob(&sc);
	return 0;
}
