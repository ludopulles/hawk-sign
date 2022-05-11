/*
 * Tests:
 * - if signing gives a valid signature,
 * - if (0, 0, .... 0) is not a valid signature.
 */
#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "../inner.h"

union {
	uint8_t b[48 * 512];
	uint64_t dummy_u64;
	fpr dummy_fpr;
} tmp;

int8_t f[512], g[512], F[512], G[512];
uint8_t h[512/4];
fpr q00[512], q10[512], q11[512], exp_sk[EXPANDED_SECKEY_SIZE(9)];
int16_t sig[512];
inner_shake256_context sc;

void test_simple_forgeries(unsigned logn)
{
	size_t n;

	n = MKN(logn);

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, tmp.b);
	Zf(expand_seckey)(exp_sk, f, g, F, logn);

	for (int rep = 0; rep < 500; rep++) {
		// make a signature of a random message...
		inner_shake256_extract(&sc, h, sizeof h);

		// Compute the signature.
		Zf(sign)(&sc, sig, exp_sk, h, logn, tmp.b);

		assert(Zf(verify_simple_rounding_fft)(h, sig, q00, q10, q11, logn, tmp.b));
		for (size_t u = 0; u < n; u ++)
			sig[u] = 0;
		assert(!Zf(verify_simple_rounding_fft)(h, sig, q00, q10, q11, logn, tmp.b));
	}
}

int main() {
	unsigned char seed[48];

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	/*
	 * For smaller values of logn, an all-zeros signature might verify to a
	 * certain hash, but we claim no security for them anyway.
	 */
	for (unsigned logn = 6; logn <= 9; logn++) {
		test_simple_forgeries(logn);
	}

	printf("No simple forgeries were possible.\n");

	return 0;
}