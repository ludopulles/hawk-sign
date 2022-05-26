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
fpr q00[512], q10[512], q11[512];
int16_t q00_i16[512], q10_i16[512], s0[512], s1[512];
int32_t q11_i32[512];
inner_shake256_context sc;

void test_verify_simple_pure(unsigned logn)
{
	size_t n;

	n = MKN(logn);

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, tmp.b);

	// make a signature of a random message...
	inner_shake256_extract(&sc, h, sizeof h);

	// Compute the signature.
	while (!Zf(uncompressed_sign)(&sc, s0, s1, f, g, F, G, h, logn, tmp.b)) {}

	/* const int max_s0 = 32, max_s1 = 16;
	for (size_t u = 0; u < n; u++) {
		s0[u] = (rand() % (2 * max_s0)) - max_s0;
		s1[u] = (rand() % (2 * max_s1)) - max_s1;
	} */

	Zf(verify_simple)(h, s0, s1, q00, q10, q11, logn, tmp.b);
	// Zf(verify_simple_rounding_fft)(h, s1, q00, q10, q11, logn, tmp.b);

	// ===== perform inverse-FFT on public key =====
	Zf(fft_to_int16)(q00_i16, q00, logn);
	Zf(fft_to_int16)(q10_i16, q10, logn);
	Zf(iFFT)(q11, logn);

	Zf(complete_pubkey2)(q00_i16, q10_i16, q11_i32, logn, tmp.b);
	for (size_t u =0 ; u < n; u++) {
		// printf("%d vs %d\n", (int) fpr_rint(q11[u]), q11_i32[u]);
		assert(fpr_rint(q11[u]) == q11_i32[u]);
	}

	uint64_t norm = Zf(verify_simple_NTT)(
		h, s0, s1, q00_i16, q10_i16, logn, tmp.b);
	printf("Norm = %lld\n", norm);
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
	for (unsigned logn = 1; logn <= 9; logn++) {
		test_verify_simple_pure(logn);
	}

	return 0;
}
