/*
 * Determine what the probability is that decompressing a signature, i.e.
 * reconstructing the first half of a signature fails, and check if this is
 * close to the expected probability acquired when taking the continuous limit.
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
// x86_64 specific:
#include <sys/time.h>

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

static int
do_sign(prng *rng, const fpr *restrict expanded_seckey,
	const uint8_t *restrict h, unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u;
	const fpr *bf, *bg, *bF, *bG, *invq00;
	fpr *x0, *x1, *res;
	int32_t z;

	n = MKN(logn);

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

	for (u = 0; u < n; u ++) {
		z = fpr_rint(x0[u]) & 1;
		z = mkgauss_sign(rng, z);
		x0[u] = fpr_of(z);
	}
	for (u = 0; u < n; u ++) {
		z = fpr_rint(x1[u]) & 1;
		z = mkgauss_sign(rng, z);
		x1[u] = fpr_of(z);
	}

	// if ((uint32_t)norm > Zf(l2bound)[logn]) return 0;

	Zf(FFT)(x0, logn);
	Zf(FFT)(x1, logn);

	Zf(poly_add_muladj_fft)(res, x0, x1, bf, bg, logn);
	Zf(poly_mul_autoadj_fft)(res, invq00, logn);
	Zf(iFFT)(res, logn);

	for (u = 0; u < n; u++) {
		if (!fpr_lt(fpr_neg(fpr_one), res[u]) || !fpr_lt(res[u], fpr_one)) {
			return 0;
		}
	}
	return 1;
}

// =============================================================================
// =============================================================================
#define LOGN (9)
#define N MKN(LOGN)

typedef long double FT;
const FT sigma_sig = 1.278;

// Returns the error function for a normal distributed variable with standard deviation sigma,
// i.e. the probability that a sample X has value that does NOT lie in (-x, x).
FT erfcl_s(FT x, FT sigma) {
	return erfcl(x / (sigma * sqrtl(2)));
}

FT fail_prob(FT iq00, long long n) {
	FT val = 1.0 - powl(erfl(0.5 / sqrtl(2*iq00) / sigma_sig), n);
	if (val < 1e-50) return n * erfcl_s(0.5, sqrt(iq00) * sigma_sig);
	return val;
}

uint8_t b[48 << LOGN];
int8_t f[N], g[N], F[N], G[N];
int16_t iq00[N], iq01[N];
fpr inv_q00[N], exp_key[9*N/2];

/*
 * This function measures the ACTUAL probability that signing gives a different
 * first half of the signature. Perhaps the norm is still too large, or
 * another reconstruction worked, we do NOT consider these cases (the prob. for
 * this is small anyways).
 * 
 * This assumes we have generated a key which is stored in the global variables
 * above.
 */
size_t measure_decompression_failure(inner_shake256_context *sc, size_t num_reps) {
	prng rng;
	uint8_t h[N/4];

	// initialize only once for performance
	Zf(prng_init)(&rng, sc);

	size_t num_fails = 0;
	for (size_t i = 0; i < num_reps; i++) {
		// Generate a random hash to sign.
		Zf(prng_get_bytes)(&rng, h, sizeof h);
		num_fails += !do_sign(&rng, exp_key, h, LOGN, b);
	}
	return num_fails;
}

void test_decomp_prob(inner_shake256_context *sc)
{
	size_t hn, u;
	hn = N/2;

	const size_t num_reps = 100 * 1000;

	while (1) {
		Zf(keygen)(sc, f, g, F, G, iq00, iq01, LOGN, b);

		Zf(int16_to_fft)(inv_q00, iq00, LOGN);
		Zf(FFT)(inv_q00, LOGN);
		for (u = 0; u < hn; u ++)
			inv_q00[u] = fpr_inv(inv_q00[u]);
		Zf(iFFT)(inv_q00, LOGN);

		double cst_iq = *((double*)&inv_q00[0]);
		FT p = fail_prob(cst_iq, N);

		// if (cst_iq < 0.002) continue;
		printf("Probability cst(1/Q00) = %.5f: %Le exp\n", cst_iq, p);
		continue;

		// if (num_reps * p >= 0.01) {
		Zf(expand_seckey)(exp_key, f, g, F, LOGN);
		size_t found = measure_decompression_failure(sc, num_reps);
		printf("Probability cst(1/Q00) = %.5f: %.5Lf exp, %.5f found\n", cst_iq, p, (double)found / num_reps);
		// }
	}
}

int main() {
	unsigned char seed[48];
	inner_shake256_context sc;

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	test_decomp_prob(&sc);
	return 0;
}
