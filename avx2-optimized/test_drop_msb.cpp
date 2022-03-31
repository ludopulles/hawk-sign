// g++ test_drop_msb.cpp build/*
#include <cassert>
#include <climits>
#include <cstdio>
#include <vector>
#include <mutex>
#include <thread>

// x86_64 specific:
#include<sys/time.h>

extern "C" {
	#ifndef restrict
		#define restrict
	#endif

	#include "inner.h"
}


typedef long long ll;

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

ll time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}


/*
 * Used for decompressing signatures. However, it seems we cannot successfully
 * recover the most significant bits of s1 with some decent probability.
 */

/*
 * Perform Fast Fourier Orthogonalization for target vector t. The Gram matrix
 * is provided (G = [[g00, g01], [adj(g01), g11]]). The sampled vector
 * is written over (t0,t1). The Gram matrix is modified as well. The
 * tmp[] buffer must have room for four polynomials.
 */
void
Zf(ffBabai2_dyn)(fpr *restrict t0, fpr *restrict t1,
	fpr *restrict g00, fpr *restrict g01, fpr *restrict g11,
	unsigned logn, fpr *restrict tmp)
{
	size_t n, hn;
	fpr *z0, *z1;

	/*
	 * Deepest level: the LDL tree leaf value is just g00 (the
	 * array has length only 1 at this point);
	 */
	if (logn == 0) {
		t0[0] = fpr_of(fpr_rint(t0[0]));
		t1[0] = fpr_of(fpr_rint(t1[0]));
		return;
	}

	n = (size_t)1 << logn;
	hn = n >> 1;

	/*
	 * Decompose G into LDL. We only need d00 (identical to g00),
	 * d11, and l10; we do that in place.
	 */
	Zf(poly_LDL_fft)(g00, g01, g11, logn);

	/*
	 * Split d00 and d11 and expand them into half-size quasi-cyclic
	 * Gram matrices. We also save l10 in tmp[].
	 */
	Zf(poly_split_fft)(tmp, tmp + hn, g00, logn);
	memcpy(g00, tmp, n * sizeof *tmp);
	Zf(poly_split_fft)(tmp, tmp + hn, g11, logn);
	memcpy(g11, tmp, n * sizeof *tmp);
	memcpy(tmp, g01, n * sizeof *g01);
	memcpy(g01, g00, hn * sizeof *g00);
	memcpy(g01 + hn, g11, hn * sizeof *g00);

	/*
	 * The half-size Gram matrices for the recursive LDL tree
	 * building are now:
	 *   - left sub-tree: g00, g00+hn, g01
	 *   - right sub-tree: g11, g11+hn, g01+hn
	 * l10 is in tmp[].
	 */

	/*
	 * We split t1 and use the first recursive call on the two
	 * halves, using the right sub-tree. The result is merged
	 * back into tmp + 2*n.
	 */
	z1 = tmp + n;
	Zf(poly_split_fft)(z1, z1 + hn, t1, logn);
	Zf(ffBabai2_dyn)(z1, z1 + hn, g11, g11 + hn, g01 + hn,
		logn - 1, z1 + n);
	Zf(poly_merge_fft)(tmp + (n << 1), z1, z1 + hn, logn);

	/*
	 * Compute tb0 = t0 + (t1 - z1) * l10.
	 * At that point, l10 is in tmp, t1 is unmodified, and z1 is
	 * in tmp + (n << 1). The buffer in z1 is free.
	 *
	 * In the end, z1 is written over t1, and tb0 is in t0.
	 */
	memcpy(z1, t1, n * sizeof *t1);
	Zf(poly_sub)(z1, tmp + (n << 1), logn);
	memcpy(t1, tmp + (n << 1), n * sizeof *tmp);
	Zf(poly_mul_fft)(tmp, z1, logn);
	Zf(poly_add)(t0, tmp, logn);

	/*
	 * Second recursive invocation, on the split tb0 (currently in t0)
	 * and the left sub-tree.
	 */
	z0 = tmp;
	Zf(poly_split_fft)(z0, z0 + hn, t0, logn);
	Zf(ffBabai2_dyn)(z0, z0 + hn, g00, g00 + hn, g01,
		logn - 1, z0 + n);
	Zf(poly_merge_fft)(t0, z0, z0 + hn, logn);
}



// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
constexpr size_t logn = 9, n = MKN(logn);
const int num_samples = 10'000;

struct WorkerResult {
	ll iters, babai_fail;

	WorkerResult() : iters(0), babai_fail(0) {}

	void combine(const WorkerResult &res) {
		iters += res.iters;
		babai_fail += res.babai_fail;
	}
};

constexpr int SCALE = 1 << 8;

WorkerResult measure_signatures(fpr isigma_kg, fpr isigma_sig, uint32_t bound)
{
	uint8_t b[100 << logn];
	int8_t f[n], g[n], F[n], G[n], h[n];
	int16_t s0[n], s1[n], rec_s0[n], rec_s1[n];
	fpr q00[n], q10[n], q11[n];
	unsigned char seed[48];
	inner_shake256_context sc;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	WorkerResult result;
	result.iters = num_samples;

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, isigma_kg, logn, b);

	for (int rep = 0; rep < num_samples; rep++) {
		// make a signature of a random message...
		random_hash(h, logn);

		// Compute the signature.
		Zf(guaranteed_sign)(&sc, s1, f, g, q00, h, isigma_sig, bound, logn, b);
		assert(Zf(verify_simple_rounding)(h, s0, s1, q00, q10, q11, bound, logn, b));
		assert(Zf(verify_nearest_plane)(h, s0, s1, q00, q10, q11, bound, logn, b));

		/*
		 * Pretend, we are dropping the most significant bits of s1
		 * However, Babai's NP is deterministic and translation invariant,
		 * so adding SCALE * f(X) to s1(X) does not change the outcome: 
		 * both differences are in the fundamental parallelepiped generated by
		 * the GSO of the basis vectors.
		 */

		fpr *t0 = (fpr *)b;
		fpr *t1 = t0 + n;
		fpr *tmp_g00 = t1 + n;
		fpr *tmp_g10 = tmp_g00 + n;
		fpr *tmp_g11 = tmp_g10 + n;

		for (size_t u = 0; u < n; u++) {
			t0[u] = fpr_half(fpr_of(h[u]));
			t1[u] = fpr_div(fpr_of(-s1[u]), fpr_of(SCALE));
		}

		Zf(FFT)(t0, logn);
		Zf(FFT)(t1, logn);

		// Now t0 = h/2, but change it to t0 = h/2 - s1/SCALE q10/q00
		// memcpy(tmp_g11, t1, n * sizeof *t1);
		// Zf(poly_mul_fft)(tmp_g11, q10, logn);
		// Zf(poly_div_autoadj_fft)(tmp_g11, q00, logn);
		// Zf(poly_sub)(t0, tmp_g11, logn);

		memcpy(tmp_g00, q00, n * sizeof *q00);
		memcpy(tmp_g10, q10, n * sizeof *q10);
		memcpy(tmp_g11, q11, n * sizeof *q11);

		Zf(poly_mulconst)(tmp_g10, fpr_of(SCALE), logn);
		Zf(poly_mulconst)(tmp_g11, fpr_of(SCALE * SCALE), logn);

		Zf(ffBabai2_dyn)(t0, t1, tmp_g00, tmp_g10, tmp_g11, logn, tmp_g11 + n);

		Zf(iFFT)(t0, logn);
		Zf(iFFT)(t1, logn);

		for (size_t u = 0; u < n; u++) {
			rec_s0[u] = fpr_rint(t0[u]);
			rec_s1[u] = s1[u] + SCALE * fpr_rint(t1[u]);
		}

/*		printf("Recovery errors:\n");
		for (size_t u = 0; u < n; u++) {
			if (s1[u] != rec_s1[u]) { // <=> t1[u] not rounded to zero
				printf("At position ");
				for (size_t bit = logn; bit --> 0; ) {
					printf("%zu", (u >> bit) & 1);
				}
				printf(" (%3zu): %d should have been %d.\n", u, rec_s1[u], s1[u]);
			}
		}
		printf("\n"); */

		bool works = Zf(verify_simple_rounding)(h, rec_s0, rec_s1, q00, q10, q11, bound, logn, b);
		// printf("Valid recovery: %d\n", works);
		result.babai_fail += !works;
	}
	return result;
}

WorkerResult tot;
std::mutex mx;

constexpr fpr sigma_kg  = { v: 1.425 };
constexpr fpr sigma_sig = { v: 1.292 };
constexpr fpr verif_margin = { v: 1.1 };

const uint32_t bound = fpr_floor(fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_of(n))));

void work()
{
	WorkerResult result = measure_signatures(fpr_inv(sigma_kg), fpr_inv(sigma_sig), bound);

	{
		std::lock_guard<std::mutex> guard(mx);
		tot.combine(result);
	}
}

int main() {
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	const int nthreads = 4;
	std::thread* pool[nthreads-1];
	for (int i = 0; i < nthreads-1; i++) pool[i] = new std::thread(work);
	work();
	for (int i = 0; i < nthreads-1; i++) pool[i]->join(), delete pool[i];

	printf("\n");
	printf("# Signatures signed:      %lld\n", tot.iters);
	printf("# Babai roundings failed: %lld\n", tot.babai_fail);

	return 0;
}
