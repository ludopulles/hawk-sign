#include <bitset>
#include <cassert>
#include <cmath>
#include <climits>
#include <cstdio>
#include <vector>

// x86_64 specific:
#include<sys/time.h>

extern "C" {
	#define restrict
	#include "../inner.h"
}

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

static inline void
int8_to_ntt(uint16_t *restrict p, const int8_t *restrict f, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u++) {
		p[u] = Zf(mq_conv_small)(f[u]);
	}
	Zf(mq_NTT)(p, logn);
}

/* =================================================================== */
constexpr size_t logn = 9, n = MKN(logn);
uint8_t tmp[100 * n];

/*
 * Given the hash of a message consisting of 2n coefficients in {0,1},
 * Calculate B * (h0, h1) = (f h0 + F h1, g h0 + G h1).
 */
static void
calc_target_NTT(int8_t *x0, int8_t *x1, const int8_t *h0, const int8_t *h1,
	int8_t *f, int8_t *g, int8_t *F, int8_t *G)
{
	size_t u;
	uint16_t *bf, *bg, *bF, *bG, *t0, *t1;

	bf = (uint16_t *)tmp;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;
	t0 = bG + n;
	t1 = t0 + n;

	int8_to_ntt(t0, h0, logn);
	int8_to_ntt(t1, h1, logn);

	int8_to_ntt(bf, f, logn);
	int8_to_ntt(bg, g, logn);
	int8_to_ntt(bF, F, logn);

	int8_to_ntt(bG, G, logn);
	// This is slower, as computing the inverse takes quite some multiplications:
	/* for (u = 0; u < n; u++) {
		bG[u] = Zf(mq_mul)(bg[u], bF[u]);
		bG[u] = Zf(mq_add)(1, bG[u]);
	}
	Zf(mq_poly_div)(bG, bf, logn); */

	Zf(mq_poly_tomonty)(t0, logn);
	Zf(mq_poly_tomonty)(t1, logn);

	/*
	 * Set the target vector to (t0, t1) = B * (h0, h1), i.e.:
	 *     t0 = f h0 + F h1,
	 *     t1 = g h0 + G h1.
	 */
	for (u = 0; u < n; u ++) {
		uint32_t res0, res1;

		res0 = Zf(mq_add)(Zf(mq_montymul)(bf[u], t0[u]),
			Zf(mq_montymul)(bF[u], t1[u]));
		res1 = Zf(mq_add)(Zf(mq_montymul)(bg[u], t0[u]),
			Zf(mq_montymul)(bG[u], t1[u]));
		t0[u] = res0;
		t1[u] = res1;
	}

	/*
	 * Sample and write the result in (x0, x1). Gaussian smoothing is used to
	 * not reveal information on the secret basis.
	 */
	Zf(mq_iNTT)(t0, logn);
	Zf(mq_iNTT)(t1, logn);

	for (u = 0; u < n; u ++) {
		x0[u] = Zf(mq_conv_signed)(t0[u]) & 1;
		x1[u] = Zf(mq_conv_signed)(t1[u]) & 1;
	}
}

static void
calc_target_FFT(int8_t *x0, int8_t *x1, const int8_t *h0, const int8_t *h1,
	int8_t *f, int8_t *g, int8_t *F, int8_t *G)
{
	size_t u;
	fpr *t0, *t1, *bf, *bg, *bF, *bG;

	t0 = (fpr *)tmp;
	t1 = t0 + n;
	bf = t1 + n;
	bg = bf + n;
	bF = bg + n;
	bG = bF + n;

	/*
	 * Set the target vector to [h0, h1] * B (hm is the hashed message).
	 */
	for (u = 0; u < n; u++) {
		t0[u] = fpr_of(h0[u]);
		t1[u] = fpr_of(h1[u]);
	}
	Zf(FFT)(t0, logn);
	Zf(FFT)(t1, logn);

	Zf(int8_to_fft)(bf, f, logn);
	Zf(int8_to_fft)(bg, g, logn);
	Zf(int8_to_fft)(bF, F, logn);

	// G is unused.
	// Zf(int8_to_fft)(bG, G, logn);
	Zf(poly_prod_fft)(bG, bg, bF, logn);
	for (u = 0; u < n/2; u++) {
		bG[u] = fpr_add(bG[u], fpr_one);
	}
	Zf(poly_div_fft)(bG, bf, logn);

	Zf(poly_matmul_fft)(bf, bF, bg, bG, t0, t1, logn);

	Zf(iFFT)(t0, logn);
	Zf(iFFT)(t1, logn);

	for (u = 0; u < n; u++) {
		x0[u] = fpr_rint(t0[u]) & 1;
		x1[u] = fpr_rint(t1[u]) & 1;
	}
}

static void
calc_target_bitset(int8_t *x0, int8_t *x1, const int8_t *h0, const int8_t *h1,
	int8_t *f, int8_t *g, int8_t *F, int8_t *G)
{

/* naive O(n^2):
 *
	memset(x0, 0, n);
	memset(x1, 0, n);

	for (size_t u = 0; u < n; u++) {
		for (size_t v = 0; v < n; v++) {
			x0[(u+v)&(n-1)] ^= ((h0[u] & f[v]) ^ (h1[u] & F[v])) & 1;
			x1[(u+v)&(n-1)] ^= ((h0[u] & g[v]) ^ (h1[u] & G[v])) & 1;
		}
	}
	return;
*/


	size_t u;
	std::bitset<2*n> t0, t1, t2, t3, t4, t5, t6, t7, t8;

	/*
	 * Set the target vector to [h0, h1] * B.
	 */
	for (u = 0; u < n; u++) {
		t0.set(u, h0[u] & 1);
		t1.set(u, h1[u] & 1);
		t2.set(u, f[u] & 1);
		t3.set(u, g[u] & 1);
		t4.set(u, F[u] & 1);
		t5.set(u, G[u] & 1);
	}

	for (u = 0; u < n; u++) {
		if (t0.test(u)) {
			t6 ^= t2 << u;
			t7 ^= t3 << u;
		}
		if (t1.test(u)) {
			t6 ^= t4 << u;
			t7 ^= t5 << u;
		}
	}

	for (u = 0; u < n; u++) {
		x0[u] = t6.test(u) ^ t6.test(u + n);
		x1[u] = t7.test(u) ^ t7.test(u + n);
	}
}

constexpr int REPETITIONS = 1000;

int8_t f[n], g[n], F[n], G[n]; // private key (basis)
int8_t h0[REPETITIONS][n], h1[REPETITIONS][n]; // hashes
int8_t x0[n], x1[n]; // output signature
int8_t expx0[n], expx1[n]; // expected answer

fpr q00[n], q01[n], q11[n];

int main() {
	unsigned char seed[48];
	inner_shake256_context sc;

	struct timeval t0, t1;
	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q01, q11, logn, tmp);

	gettimeofday(&t1, NULL);
	printf("Key generation took %lld microseconds\n", time_diff(&t0, &t1));

	// =========================================================================
	// | Benchmark the signing of random messages                              |
	// =========================================================================

	inner_shake256_extract(&sc, (uint8_t *)&h0[0][0], REPETITIONS * n);
	inner_shake256_extract(&sc, (uint8_t *)&h1[0][0], REPETITIONS * n);
	for (int rep = 0; rep < REPETITIONS; rep++) {
		for (size_t u = 0; u < n; u++) {
			h0[rep][u] &= 1;
			h1[rep][u] &= 1;
		}
	}

	for (int rep = 0; rep < std::min(100, REPETITIONS); rep++) {
		memset(expx0, 0, n);
		memset(expx1, 0, n);

		for (size_t u = 0; u < n; u++) {
			for (size_t v = 0; v < n; v++) {
				expx0[(u+v)&(n-1)] ^= ((h0[rep][u] & f[v]) ^ (h1[rep][u] & F[v])) & 1;
				expx1[(u+v)&(n-1)] ^= ((h0[rep][u] & g[v]) ^ (h1[rep][u] & G[v])) & 1;
			}
		}

		calc_target_NTT(x0, x1, h0[rep], h1[rep], f, g, F, G);
		for (size_t u = 0; u < n; u++)
			assert(x0[u] == expx0[u] && x1[u] == expx1[u]);
		calc_target_FFT(x0, x1, h0[rep], h1[rep], f, g, F, G);
		for (size_t u = 0; u < n; u++)
			assert(x0[u] == expx0[u] && x1[u] == expx1[u]);
		calc_target_bitset(x0, x1, h0[rep], h1[rep], f, g, F, G);
		for (size_t u = 0; u < n; u++)
			assert(x0[u] == expx0[u] && x1[u] == expx1[u]);
	}

	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < REPETITIONS; rep++) {
		calc_target_NTT(x0, x1, h0[rep], h1[rep], f, g, F, G);
	}
	gettimeofday(&t1, NULL);
	printf("NTT %lld us\n", time_diff(&t0, &t1) / REPETITIONS);

	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < REPETITIONS; rep++) {
		calc_target_FFT(x0, x1, h0[rep], h1[rep], f, g, F, G);
	}
	gettimeofday(&t1, NULL);
	printf("FFT %lld us\n", time_diff(&t0, &t1) / REPETITIONS);

	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < REPETITIONS; rep++) {
		calc_target_bitset(x0, x1, h0[rep], h1[rep], f, g, F, G);
	}
	gettimeofday(&t1, NULL);
	printf("bitset %lld us\n", time_diff(&t0, &t1) / REPETITIONS);

	return 0;
}
