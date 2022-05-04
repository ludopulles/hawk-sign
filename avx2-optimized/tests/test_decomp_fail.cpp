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
			z = mkgauss_1292(rng, fpr_rint(x0[u]) & 1);
			norm += z*z;
			x0[u] = fpr_of(z);
		}
		for (u = 0; u < n; u ++) {
			z = mkgauss_1292(rng, fpr_rint(x1[u]) & 1);
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
constexpr FT sigma_sig = 1.292;
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
