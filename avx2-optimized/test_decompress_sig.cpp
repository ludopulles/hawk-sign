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
	#include "inner.h"
}

// =============================================================================
// SIGNING
// =============================================================================
// Taken from sign.c:
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

static inline int mkgauss_1292(prng *rng, uint8_t double_mu) {
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
		z = mkgauss_1292(rng, z);
		x0[u] = fpr_of(z);
	}
	for (u = 0; u < n; u ++) {
		z = fpr_rint(x1[u]) & 1;
		z = mkgauss_1292(rng, z);
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
const FT sigma_sig = 1.292;

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
fpr q00[N], q10[N], q11[N], iq00[N], exp_key[9*N/2];

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
		Zf(keygen)(sc, f, g, F, G, q00, q10, q11, LOGN, b);

		for (u = 0; u < hn; u ++) {
			iq00[u] = fpr_inv(q00[u]);
			iq00[u + hn] = fpr_zero;
		}
		Zf(iFFT)(iq00, LOGN);

		double cst_iq = *((double*)&iq00[0]);
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
