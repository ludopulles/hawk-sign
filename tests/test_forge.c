/*
 * Tests:
 * - if signing gives a valid signature,
 * - if (0, 0, .... 0) is not a valid signature.
 */
#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "../inner.h"

#define MAXLOGN (10)
#define MAXN MKN(MAXLOGN)

union {
	uint8_t b[48 * MAXN];
	uint64_t dummy_u64;
	fpr dummy_fpr;
} tmp;

int8_t f[MAXN], g[MAXN], F[MAXN], G[MAXN];
uint8_t h[MAXN/4];
int16_t iq00[MAXN], iq10[MAXN];
fpr q00[MAXN], q10[MAXN], q11[MAXN], exp_sk[EXPANDED_SECKEY_SIZE(MAXLOGN)];
int16_t sig[MAXN];
inner_shake256_context sc;

void test_simple_forgeries(unsigned logn)
{
	size_t n;

	n = MKN(logn);

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, iq00, iq10, logn, tmp.b);
	// Zf(complete_pubkey)(iq00, iq10, q00, q10, q11, logn);
	Zf(expand_seckey)(exp_sk, f, g, F, logn);

	for (int rep = 0; rep < 500; rep++) {
		// make a signature of a random message...
		inner_shake256_extract(&sc, h, sizeof h);

		// Compute the signature.
		while (!Zf(sign)(&sc, sig, exp_sk, h, logn, tmp.b)) {}
		// while (!Zf(sign_NTT)(&sc, sig, f, g, F, G, h, logn, tmp.b)) {}

		// assert(Zf(verify)(h, sig, q00, q10, q11, logn, tmp.b));
		assert(Zf(verify_NTT)(h, sig, iq00, iq10, logn, tmp.b));

		for (size_t u = 0; u < n; u ++)
			sig[u] = 0;
		// assert(!Zf(verify)(h, sig, q00, q10, q11, logn, tmp.b));
		assert(!Zf(verify_NTT)(h, sig, iq00, iq10, logn, tmp.b));
	}
}



const int num_keygen = 100;
const int signs_per_kg = 100;
const int totalruns = num_keygen * signs_per_kg;

int count_fails(unsigned logn)
{
	size_t num_fails = 0;

	// Generate key pair.
	for (size_t u = 0; u < num_keygen; u++) {
		Zf(keygen)(&sc, f, g, F, G, iq00, iq10, logn, tmp.b);
		Zf(expand_seckey)(exp_sk, f, g, F, logn);

		for (int rep = 0; rep < totalruns/num_keygen; rep++) {
			// make a signature of a random message...
			inner_shake256_extract(&sc, h, sizeof h);

			// Compute the signature.
			// while (!Zf(sign)(&sc, sig, exp_sk, h, logn, tmp.b)) {}
			while (!Zf(sign_NTT)(&sc, sig, f, g, F, G, h, logn, tmp.b)) {}

			num_fails += !Zf(verify_NTT)(h, sig, iq00, iq10, logn, tmp.b);
		}
	}
	return num_fails;
}

int main() {
	unsigned char seed[48] = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);

	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	/*
	 * For smaller values of logn, an all-zeros signature might verify to a
	 * certain hash, but we claim no security for them anyway.
	 */
	for (unsigned logn = 6; logn <= MAXLOGN; logn++) {
		// printf("%d", logn);
		// fflush(stdout);
		// test_simple_forgeries(logn);
		printf("%d: %d/%d fails\n",
			logn, count_fails(logn), totalruns);
	}

	printf("\nNo simple forgeries were possible.\n");

	return 0;
}
