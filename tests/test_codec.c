/*
 * Tests whether first encoding and then decoding a public key gives the same
 * result back, and similar for the secret key.
 */

#include <assert.h>
#include <stdio.h>

#include "../inner.h"

int poly_eq(int8_t *p, int8_t *q, unsigned logn) {
	for (size_t u = 0; u < MKN(logn); u++)
		if (p[u] != q[u]) return 0;
	return 1;
}

int poly16_eq(int16_t *p, int16_t *q, unsigned logn) {
	for (size_t u = 0; u < MKN(logn); u++)
		if (p[u] != q[u]) return 0;
	return 1;
}


#define LOGN (10)
#define N MKN(LOGN)
#define BUFLEN (10240)

/*
 * Maximum practical signature size (in bytes) when using the compact format
 * and the uncompressed scheme.
 */
#define HAWK_UNCOMPRESSED_SIG_COMPACT_MAXSIZE(logn) \
	((logn) == 10 ? 2863u : 1309u)

/*
 * Maximum practical signature size (in bytes) when using the COMPACT
 * format.
 */
#define HAWK_SIG_COMPACT_MAXSIZE(logn) \
	((logn) == 10 ? 1263u : 593u)

uint8_t b[48 << LOGN], outbuf[BUFLEN];
int8_t f[N], g[N], F[N], G[N];
fpr q00[N], q01[N], q11[N];
int16_t iq00[N], iq01[N];

uint8_t h[N/4];
int16_t s0[N], s1[N], _s0[N], _s1[N];

int8_t _f[N], _g[N], _F[N], _G[N];
int16_t _iq00[N], _iq01[N];
fpr exp_seckey[EXPANDED_SECKEY_SIZE(LOGN)];

inner_shake256_context sc;
const int n_repetitions = 100;

void test_encode_decode(unsigned logn) {
	size_t n;
	size_t ls, lus;

	n = MKN(logn);
	ls = HAWK_SIG_COMPACT_MAXSIZE(logn);
	lus = HAWK_UNCOMPRESSED_SIG_COMPACT_MAXSIZE(logn);

	for (int i = 0; i < n_repetitions; i++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, iq00, iq01, logn, b);

		// Test secret key encoding and decoding
		size_t sk_len = Zf(encode_seckey)(outbuf, BUFLEN, f, g, F, logn);
		assert(sk_len != 0);
		assert(sk_len == Zf(decode_seckey)(_f, _g, _F, outbuf, sk_len, logn));

		assert(poly_eq(f, _f, logn));
		assert(poly_eq(g, _g, logn));
		assert(poly_eq(F, _F, logn));

		// Test expanding secret key
		Zf(expand_seckey)(exp_seckey, _f, _g, _F, logn);
		Zf(iFFT)(exp_seckey + (3u << logn), logn);
		for (size_t u = 0; u < n; u++) {
			_G[u] = fpr_rint(exp_seckey[(3u << logn) + u]);
		}
		assert(poly_eq(G, _G, logn));

		// Test public key encoding and decoding
		Zf(fft_to_int16)(iq00, q00, logn);
		Zf(fft_to_int16)(iq01, q01, logn);
		Zf(FFT)(q00, logn);
		Zf(FFT)(q01, logn);

		size_t pk_len = Zf(encode_pubkey)(outbuf, BUFLEN, iq00, iq01, logn);
		assert(pk_len != 0);
		assert(pk_len == Zf(decode_pubkey)(_iq00, _iq01, outbuf, pk_len, logn));

		assert(poly16_eq(iq00, _iq00, logn));
		assert(poly16_eq(iq01, _iq01, logn));

		// Test complete signature encoding and decoding
		inner_shake256_extract(&sc, h, n <= 8 ? 2 : n / 4);

		while (!Zf(uncompressed_sign)(&sc, s0, s1, f, g, F, G, h, logn, b)) {}
		size_t sig_len = Zf(encode_uncomp_sig)(outbuf, lus, s0, s1,
			logn);
		assert(sig_len != 0);
		assert(sig_len == Zf(decode_uncomp_sig)(_s0, _s1, outbuf, sig_len,
			logn));

		assert(poly16_eq(s0, _s0, logn));
		assert(poly16_eq(s1, _s1, logn));

		// Test (compressed) signature encoding and decoding
		while (!Zf(sign_dyn)(&sc, s1, f, g, F, G, h, logn, b)) {}
		sig_len = Zf(encode_sig)(outbuf, ls, s1, logn);
		assert(sig_len != 0);
		assert(sig_len == Zf(decode_sig)(_s1, outbuf, sig_len, logn));

		assert(poly16_eq(s1, _s1, logn));
	}
}

int main() {
	// Initialize a RNG.
	unsigned char seed[48];
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	for (unsigned logn = 1; logn <= LOGN; logn++) {
		printf("%u", logn);
		fflush(stdout);
		test_encode_decode(logn);
	}
	printf("\n");
	return 0;
}
