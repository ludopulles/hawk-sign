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

static inline int
mkgauss_1292(prng *rng, uint8_t double_mu)
{
	uint64_t r;
	uint32_t f, v, k, neg;

	r = prng_get_u64(rng);
	neg = (uint32_t)(r >> 63);
	r &= ~((uint64_t)1 << 63);
	f = (uint32_t)((r - gauss_1292[double_mu]) >> 63);

	v = double_mu;
	r = prng_get_u64(rng);

	r &= ~((uint64_t)1 << 63);
	for (k = 1; k < 13; k ++) {
		uint32_t t;
		t = (uint32_t)((r - gauss_1292[2 * k + double_mu]) >> 63) ^ 1;
		v |= (k << 1) & -(t & (f ^ 1));
		f |= t;
	}

	v = (v ^ -neg) + neg;
	return *(int32_t *)&v;
}

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
		z = mkgauss_1292(rng, z);
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
		z = mkgauss_1292(rng, z);
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
uint32_t bound = (uint32_t)(SQR(1.1 * 2 * 1.292) * (2*n));

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
	long double sigma_ei = sqrtl(cstQ00.v) * 1.292;
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
