#include <cassert>
#include <cstdio>
#include <climits>
// x86_64 specific:
#include <sys/time.h>

// concurrency:
// #include <thread>
// #include <mutex>

extern "C" {
	#ifndef restrict
		#define restrict
	#endif

	#include "inner.h"
}

#include "inner.h"

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
// | TESTING CODE                                                              |
// =============================================================================
const size_t logn = 9, n = MKN(logn);

void find_forgeries(fpr isigma_kg, fpr isigma_sig, uint32_t bound)
{
	union { // use union to ensure alignment
		uint8_t b[42*512];
		uint64_t dummy_i64;
		fpr dummy_fpr;
	} tmp;
	int8_t f[n], g[n], F[n], G[n], h[n]; // h2[n];
	fpr q00[n], q10[n], q11[n];
	int16_t s0[n], s1[n];
	unsigned char seed[48];
	inner_shake256_context sc;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, isigma_kg, logn, tmp.b);

	for (int rep = 0; rep < 1000; rep++) {
		// make a signature of a random message...
		random_hash(h, logn);

		// Compute the signature.
		Zf(sign)(&sc, s1, f, g, h, isigma_sig, bound, logn, tmp.b);
		assert(Zf(verify)(h, s0, s1, q00, q10, q11, bound, logn, tmp.b));
	}
	exit(1);

	int nreps = 0, nworks = 0;
	for (size_t i = 0; i < n; i++) {
		h[i] ^= 1; // change some bytes of h
		nreps++;
		nworks += Zf(verify)(h, s0, s1, q00, q10, q11, bound, logn, tmp.b);
		h[i] ^= 1; // restore
	}
	for (size_t i = 0; i < n; i++) {
		h[i] ^= 1; // change a bit of h
		for (size_t j = 0; j < i; j++) {
			h[j] ^= 1; // change a bit of h
			nreps++;
			nworks += Zf(verify)(h, s0, s1, q00, q10, q11, bound, logn, tmp.b);
			h[j] ^= 1; // restore
		}
		h[i] ^= 1; // restore
	} // */

	/* for (int rep = 0; rep < 100'000; rep++) {
		// Generate random hash again, and check if this s1 works for it as well:
		random_hash(h2, logn);

		nreps++;
		nworks += Zf(verify)(h2, s0, s1, q00, q10, q11, bound, logn, tmp.b);
	} //*/
	printf("Succesful signatures found: %d / %d\n", nworks, nreps);
	printf("All went succesful!\n");
}

int8_t valid_sigma(fpr sigma_sig)
{
	return !fpr_lt(sigma_sig, fpr_sigma_min[logn])
		&& fpr_lt(sigma_sig, fpr_div(fpr_of(18205), fpr_of(10000)));
}

constexpr fpr sigma_kg  = { v: 1.425 };
constexpr fpr sigma_sig = { v: 1.292 };
constexpr fpr verif_margin = { v: 1.1 }; // sigma_kg.v / sigma_sig.v };

int main()
{
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	assert(valid_sigma(sigma_kg));

	// expected l_2 norm of a Gaussian of std.dev. 2sigma_sig and dimension 2n is 2 sigma_sig √2n
	// (verif_margin * 2 * sigma_sig)^2 * 2n
	uint32_t bound = fpr_floor(fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_of(n))));

	find_forgeries(fpr_inv(sigma_kg), fpr_inv(sigma_sig), bound);
	return 0;
}

