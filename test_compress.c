/*
 * Tests whether first encoding and then decoding a public key gives the same
 * result back, and similar for the secret key.
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
// x86_64 specific:
#include <sys/time.h>

#include "inner.h"

#define LOGN (9)
#define N MKN(LOGN)

#define BUFLEN (10240)
uint8_t b[48 << LOGN], outbuf[BUFLEN];
int8_t f[N], g[N], F[N], G[N];
fpr q00[N], q10[N], q11[N];
int16_t iq00[N], iq10[N];

int8_t _f[N], _g[N], _F[N], _G[N];
int16_t _iq00[N], _iq10[N];
fpr _q00[N], _q10[N], _q11[N], exp_seckey[EXPANDED_SECKEY_SIZE(LOGN)];

int poly_eq(int8_t *p, int8_t *q, unsigned logn) {
	/* for (size_t u = 0; u < MKN(logn); u++) printf("%d ", (int) p[u]);
	printf("\n");
	for (size_t u = 0; u < MKN(logn); u++) printf("%d ", (int) q[u]);
	printf("\n"); */
	for (size_t u = 0; u < MKN(logn); u++) {
		if (p[u] != q[u]) return 0;
	}
	return 1;
}

void test_encode_decode(unsigned logn) {
	unsigned char seed[48];
	inner_shake256_context sc;
	const int n_repetitions = 100;

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	for (int i = 0; i < n_repetitions; i++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b);

		size_t sk_len = Zf(encode_seckey)(outbuf, BUFLEN, f, g, F, logn);
		assert(sk_len != 0);
		assert(sk_len == Zf(decode_seckey)(_f, _g, _F, outbuf, sk_len, logn));

		assert(poly_eq(f, _f, logn));
		assert(poly_eq(g, _g, logn));
		assert(poly_eq(F, _F, logn));

		Zf(expand_seckey)(exp_seckey, _f, _g, _F, logn);
		Zf(iFFT)(exp_seckey + (3u << logn), logn);
		for (size_t u = 0; u < MKN(logn); u++) {
			_G[u] = fpr_rint(exp_seckey[(3u << logn) + u]);
		}
		assert(poly_eq(G, _G, logn));

		Zf(fft_to_int16)(iq00, q00, logn);
		Zf(fft_to_int16)(iq10, q10, logn);
		Zf(FFT)(q00, logn);
		Zf(FFT)(q10, logn);

		size_t pk_len = Zf(encode_pubkey)(outbuf, BUFLEN, iq00, iq10, logn);
		assert(pk_len != 0);
		assert(pk_len == Zf(decode_pubkey)(_iq00, _iq10, outbuf, pk_len, logn));
		Zf(complete_pubkey)(_iq00, _iq10, _q00, _q10, _q11, logn);

		for (size_t u = 0; u < MKN(logn); u++) {
			assert(fpr_rint(fpr_sub(q00[u], _q00[u])) == 0);
			assert(fpr_rint(fpr_sub(q10[u], _q10[u])) == 0);
			assert(fpr_rint(fpr_sub(q11[u], _q11[u])) == 0);
		}
	}
}

int main() {
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	for (unsigned logn = 1; logn <= LOGN; logn++) {
		printf("Testing logn = %u\n", logn);
		test_encode_decode(logn);
	}
	return 0;
}
