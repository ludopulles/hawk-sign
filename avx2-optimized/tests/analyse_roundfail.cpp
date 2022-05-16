#include <cassert>
#include <climits>
#include <cstdio>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>

#include <atomic>
#include <thread>
#include <mutex>

// Determine prob. signature fails, when fixing given (f, g).

// x86_64 specific:
#include<sys/time.h>

extern "C" {
	#define restrict
	#define HAWK_NO_SIGN_NORM_CHECK
	#include "../inner.h"
}

using namespace std;
typedef long long ll;
typedef pair<double, double> pt;

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

bool mock_sign(prng *rng, const fpr *restrict expanded_seckey, const uint8_t *restrict h,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp, vector<double> &lens)
{
	size_t n, u, v, w;
	int8_t *x0, *x1;
	const fpr *tf, *tg, *tF, *tG, *invq00;
	fpr *tx0, *tx1, *tres;
	int32_t norm, z;

begin_sign:
	n = MKN(logn);
	norm = 0;

	tf = expanded_seckey;
	tg = tf + n;
	tF = tg + n;
	tG = tF + n;
	invq00 = tG + n;

	tx0 = (fpr *)tmp;
	tx1 = tx0 + n;
	tres = tx1 + n;
	x0 = (int8_t *)(tres + n);
	x1 = x0 + n;

	/*
	 * Set the target vector to (h0, h1) B (row-notation).
	 */
	for (u = 0, w = 0; w < n; u ++ ) {
		uint8_t h0, h1;

		h0 = h[u];
		h1 = h[n / 8 + u];
		for (v = 0; v < 8; v ++, w ++) {
			tx0[w] = fpr_of(h0 & 1);
			tx1[w] = fpr_of(h1 & 1);
			h0 >>= 1;
			h1 >>= 1;
		}
	}

	/*
	 * Sample (x0, x1) according to a discrete gaussian distribution on 2
	 * Z^{2n} + (h0, h1) B with standard deviation 2 sigma_{sig}. Gaussian
	 * smoothing is used to not reveal information on the secret basis.
	 */
	Zf(FFT)(tx0, logn);
	Zf(FFT)(tx1, logn);

	/*
	 * First, calculate the 0th component of (h0, h1) B, and sample x0.
	 */
	Zf(poly_add_mul_fft)(tres, tx0, tx1, tf, tF, logn);
	Zf(iFFT)(tres, logn);

	for (u = 0; u < n; u ++) {
		z = fpr_rint(tres[u]) & 1;
		z = mkgauss_sign(rng, z);
		x0[u] = (int8_t) z;
		norm += z*z;
	}

	/*
	 * Second, calculate the 1th component of (h0, h1) B, and sample x1.
	 */
	Zf(poly_add_mul_fft)(tres, tx0, tx1, tg, tG, logn);
	Zf(iFFT)(tres, logn);

	for (u = 0; u < n; u ++) {
		z = fpr_rint(tres[u]) & 1;
		z = mkgauss_sign(rng, z);
		x1[u] = (int8_t) z;
		norm += z*z;
	}

	/*
	 * Test whether the l2-norm of (x0, x1) is below the given bound. The
	 * code below uses only 32-bit operations to compute the squared norm,
	 * since the max. value is 2n * 128^2 <= 2^24 (when logn <= 9).
	 * For a large enough verification margin, it is unlikely that the
	 * norm of the gaussian (x0, x1) is too large.
	 */
	if ((uint32_t)norm >= bound) {
		// Norm is too large, so signature would not be valid.

		// Normally we would fail here...
		goto begin_sign;
	}

	Zf(int8_to_fft)(tx0, x0, logn);
	Zf(int8_to_fft)(tx1, x1, logn);

	/*
	 * Calculate the rounding errors that occur when we want to recover s0 from
	 * s1 during verification of this signature. These errors should be of
	 * absolute value at most 1/2.
	 * The errors are calculated by:
	 *
	 *     (f^* x0 + g^* x1) / (f^* f + g^* g).
	 *
	 * If err / 2 is not in (-.5,.5)^n, verification will reconstruct a
	 * different s0 so the signature may not be valid anymore.
	 */
	Zf(poly_add_muladj_fft)(tres, tx0, tx1, tf, tg, logn);
	Zf(poly_mul_autoadj_fft)(tres, invq00, logn);
	Zf(iFFT)(tres, logn);

	// double norm2 = 0.0;
	// for (u = 0; u < n; u++) norm2 += tres[u].v * tres[u].v;
	// lens.push_back(norm2);
	// lens.emplace_back( tres[0].v, tres[n/2].v );
	// return fpr_lt(fpr_neg(fpr_one), tres[0]) && fpr_lt(tres[0], fpr_one);

	int nr_fail = 0;
	for (u = 0; u < n; u++) {
		if (!fpr_lt(fpr_neg(fpr_one), tres[u]) || !fpr_lt(tres[u], fpr_one)) {
			nr_fail++;
		}
	}
	if (nr_fail >= 1) {
		// printf("%d ", nr_fail);
	}
	return nr_fail == 0;
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
constexpr int logn = 9, n = MKN(logn);

void output_poly(int8_t *f) {
	for (size_t u = 0; u < n; u++) printf("%d ", f[u]);
	printf("\n");
}

// int8_t f[512], g[512], F[512], G[512];

int8_t f[512] = { 0,1,1,-3,1,0,-1,-2,-1,0,0,1,1,0,0,1,-1,1,-1,-1,-1,2,0,-2,3,0,-1,0,0,0,2,-2,-1,0,0,1,0,-1,1,2,-1,0,3,-4,-2,0,1,1,1,-2,-1,0,1,2,0,2,0,0,1,0,2,0,-2,1,1,0,1,-3,1,1,1,0,-4,0,1,1,2,-4,0,2,-1,-1,-1,1,0,1,1,0,-1,1,0,1,1,1,4,1,-1,-1,1,-2,-1,0,0,-1,1,2,0,1,-1,1,-1,1,-1,-2,0,0,-1,-1,-1,0,-1,2,-1,0,-1,0,1,-1,0,1,-1,-1,1,2,0,1,1,2,3,2,-1,0,-1,-3,-2,2,-1,-1,1,-2,1,1,-2,0,-2,-1,-1,-1,1,2,0,0,0,0,-3,1,1,1,0,2,0,2,-1,-1,-3,1,2,2,1,-2,3,3,-2,0,-1,1,-1,0,-2,0,-2,-2,0,3,-1,0,0,1,-1,3,-2,-1,-2,0,0,-1,1,0,-2,3,0,2,2,0,1,3,0,0,-1,-1,-2,-1,-1,0,-1,0,0,3,-2,-1,0,2,-2,1,-1,-1,0,-1,1,0,0,-2,1,1,1,-1,2,-1,0,2,-1,-1,0,1,0,-2,-1,0,1,0,1,0,1,-1,0,0,-2,0,-2,-2,-1,0,0,2,-2,-3,-4,-2,-1,-2,2,1,1,1,-2,1,3,0,-1,-4,3,-1,1,1,1,0,-3,-2,2,1,-3,-3,0,1,-1,-1,-1,2,-1,1,1,-2,0,-1,1,1,1,2,-2,-1,2,0,-2,-1,-2,1,3,1,1,-1,0,-2,0,0,-1,1,0,0,-2,2,3,-2,0,0,-1,-1,0,2,0,-1,1,1,-2,-2,1,0,3,1,-1,-1,-1,2,1,0,-2,1,0,0,1,-4,1,1,0,1,-3,2,1,1,3,0,-3,2,2,2,1,-2,0,0,0,-1,-1,-1,-3,0,0,0,1,-1,-2,1,2,-2,1,1,1,0,-1,2,1,-1,0,1,-3,-1,1,-1,-1,2,-1,2,1,2,0,-1,1,0,1,-1,-1,2,0,0,1,0,-2,0,-1,0,1,-1,1,0,0,1,1,0,1,1,0,0,0,-1,1,-2,0,2,1,2,2,0,-1,1,2,1,1,-1,0,-4,-1,-1,1,2,0,0,1,2,-1,0,1,-1,-1,2,-1,-1,-3,1,0,-1,-2,1,2,-1,1,-2,1,2,1,-2,-1,-1,-3,0,1,0,1,1,0,3,-1,3,0,1 };
int8_t g[512] = { 1,1,1,-1,0,-1,0,-2,-1,1,1,-2,0,-1,1,-1,4,1,0,0,-1,2,1,2,-1,1,0,2,-1,-5,-2,1,0,2,3,1,0,-1,2,-2,0,-1,0,2,-1,-2,0,0,1,-1,-3,0,-2,-1,1,-1,1,1,1,0,-2,-2,1,1,1,0,-2,0,3,-2,-1,0,1,0,-1,1,-1,0,-1,1,0,-2,-1,-2,0,1,0,-1,3,-1,1,2,1,1,-1,-1,0,1,-2,0,1,0,-2,0,0,3,0,0,-1,0,-1,-1,0,-1,-1,0,2,-1,-2,-1,0,-1,-2,0,-2,1,-1,-1,2,-2,1,-2,0,2,0,1,0,1,0,1,-2,-2,0,0,1,1,3,0,4,-1,-2,0,-1,-1,0,0,0,-1,0,0,1,1,0,2,-3,-1,1,1,-2,-1,1,0,-1,0,1,3,-1,0,-2,-1,0,0,-2,0,2,0,-2,0,0,-2,-1,2,2,-1,1,2,1,-3,0,-1,2,4,1,1,0,-1,0,1,0,2,0,-1,-1,1,0,-1,-1,-1,0,-1,1,1,0,0,-1,0,0,-1,1,0,0,1,-1,4,-3,0,0,0,1,-3,0,1,-1,-1,-1,1,3,4,1,0,0,2,0,-1,0,1,-1,2,-1,0,-1,0,0,0,-1,-4,-1,0,-2,0,2,-2,0,0,0,0,-2,1,1,1,-1,0,0,0,2,1,1,0,1,1,-1,1,-1,1,-2,1,0,-3,1,-1,-1,-1,-1,2,-4,0,0,-2,-2,1,1,-3,1,-1,1,-2,0,1,1,-2,-3,-2,-1,0,-1,1,-1,-1,-1,-1,-1,0,-2,-3,2,2,-1,-1,-3,-1,0,-3,3,0,-1,0,0,1,0,1,0,0,-1,1,0,-1,0,1,0,-1,-2,2,-3,2,0,1,-2,0,-1,1,0,-1,2,1,2,0,0,-1,1,3,3,3,-1,-1,2,0,-2,1,3,1,1,-1,0,-1,-1,-2,1,1,-1,0,0,2,-1,0,0,1,0,-1,1,0,0,-1,1,1,0,0,-1,2,-1,0,0,-2,-1,0,-1,0,-1,-1,4,-1,-3,0,-1,1,-1,-1,-2,3,2,0,-3,-1,-2,0,-1,-1,0,-1,0,-1,1,1,3,3,1,0,1,1,0,-1,-1,0,0,3,-1,1,-1,3,-1,1,-1,0,0,-1,1,-1,2,-1,-1,1,1,0,2,-3,3,0,3,0,-1,0,1,-2,2,1,-1,1,-1,1,1,2,2,1,3,-2,0,2,-1,0,1,0,0,-1 };
int8_t F[512], G[512];

fpr q00[512], q10[512], q11[512], expkey[EXPANDED_SECKEY_SIZE(logn)];

#define SQR(x) ((x) * (x))
uint32_t bound = (uint32_t)(SQR(1.1 * 2 * 1.278) * (2*n));

const long long num_repetitions = 100'000;

atomic<int> result;

vector<double> lens;
mutex lens_mx;

vector<array<unsigned char, 48>> seeds;

void work(int id) {
	uint8_t tmp[50 << logn];
	uint8_t h[512/4];
	vector<double> _lens;

	inner_shake256_context rng;
	inner_shake256_init(&rng);
	inner_shake256_inject(&rng, (unsigned char *)&seeds[id][0], 48);
	inner_shake256_flip(&rng);

	prng p;
	Zf(prng_init)(&p, &rng);

	int fails = 0;
	for (int rep = num_repetitions; rep -- > 0; ) {
		// inner_shake256_extract(&rng, h, n / 4);
		Zf(prng_get_bytes)(&p, h, n / 4);

		// Make sure that sign.c may fail on generating a signature that does not decompress correctly.
		// Zf(sign)(&rng, sig, expkey, h, logn, tmp);
		// if (!Zf(verify_simple_rounding_fft)(h, sig, q00, q10, q11, logn, tmp)) fails++;
		fails += !mock_sign(&p, expkey, h, bound, logn, tmp, _lens);
	}

	result += fails;
	{ // acquire mutex
		lock_guard<mutex> lock(lens_mx);
		lens.insert(lens.end(), _lens.begin(), _lens.end());
	}
}

int main(int argc, char **argv) {
	// Initialize a RNG.
	inner_shake256_context sc;
	inner_shake256_init(&sc);
	struct timeval tv;
	gettimeofday(&tv, NULL);
	inner_shake256_inject(&sc, (unsigned char *)&tv, sizeof tv);
	inner_shake256_flip(&sc);

	int numsamples = 0;
	fpr sum_cst = fpr_zero;

	fpr cstQ00;

	do {
		uint8_t tmp[50 << logn];

		// Generate key pair.
		// assert( Zf(complete_private)(f, g, F, G, q00, q10, q11, logn, tmp) );
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, tmp);

		for (size_t u = 0; u < n/2; u++) {
			q00[u] = fpr_inv(q00[u]);
			q00[u + n/2] = fpr_zero;
		}
		Zf(iFFT)(q00, logn);

		sum_cst = fpr_add(sum_cst, q00[0]);
		numsamples++;

		cstQ00 = q00[0];
		// printf("cst(1/Q00) = %.8f\n", cstQ00.v);
		// printf("."); fflush(stdout);
	// } while (cstQ00.v < 0.010);
	} while (false);

	// output_poly(f); output_poly(g);

	Zf(FFT)(q00, logn); // restore value.
	for (size_t u = 0; u < n/2; u++) q00[u] = fpr_inv(q00[u]);
	printf("\ncst(1/Q00) = %.8f\n", cstQ00.v);
	// printf("Average cst(1/Q00) = %.8f\n", fpr_div(sum_cst, fpr_of(numsamples)).v);

	Zf(expand_seckey)(expkey, f, g, F, logn);

	// compute approximate fail rate.
	long double sigma_ei = sqrtl(cstQ00.v) * 1.278;
	long double success_prob = powl(erfl(0.5 / sqrt(2) / sigma_ei), 512);
	long double fail_prob = erfcl(0.5 / sqrt(2) / sigma_ei) * 512;
	printf("Success probability = %.20Lf, fail prob = %.20Lf VS %.20Lf\n", success_prob, 1.0 - success_prob, fail_prob);

	const int nthreads = 4;
	for (int i = 0; i < nthreads; i++) {
		array<unsigned char, 48> seed;
		inner_shake256_extract(&sc, &seed[0], 48);
		seeds.push_back(seed);
	}

	std::thread* pool[nthreads];
	for (int i = 1; i < nthreads; i++) {
		pool[i] = new std::thread(work, i);
	}
	work(0);
	for (int i = 1; i < nthreads; i++) pool[i]->join(), delete pool[i];

	int ans = result;
	long double prob = (long double) ans / (nthreads * num_repetitions);

	// Print information about the correlation of errors.
	for (auto &p : lens) {
		printf("%.4f\n", p);
	}

	printf("%d / %lld failed\n", ans, nthreads * num_repetitions);
	printf("Probability: %.20Lf\n", prob);

	return 0;
}
