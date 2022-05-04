// c99 -W -Wall -O3 -march=native test_sizes.c -o test_sizes build/* -lm
#include <assert.h>
#include <stdio.h>
#include <stdint.h>

#include "../inner.h"

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
#define MAX_LOGN (9)
#define MAX_N MKN(MAX_LOGN)

long long pk_sum[MAX_LOGN+1] = {}, pk_sumsq[MAX_LOGN+1] = {}, pk_min[MAX_LOGN+1] = {}, pk_max[MAX_LOGN+1] = {};
long long sk_sum[MAX_LOGN+1] = {}, sk_sumsq[MAX_LOGN+1] = {}, sk_min[MAX_LOGN+1] = {}, sk_max[MAX_LOGN+1] = {};

long long fg_sumsq[MAX_LOGN+1] = {}, FG_sumsq[MAX_LOGN+1] = {};
long long q00_sum[MAX_LOGN+1] = {}, q00_sumsq[MAX_LOGN+1] = {};
long long q10_sum[MAX_LOGN+1] = {}, q10_sumsq[MAX_LOGN+1] = {};
long long q11_sum[MAX_LOGN+1] = {}, q11_sumsq[MAX_LOGN+1] = {};

uint8_t b[48 << MAX_LOGN];
int8_t f[MAX_N], g[MAX_N], F[MAX_N], G[MAX_N];
fpr q00[MAX_N], q10[MAX_N], q11[MAX_N];
int16_t q00n[MAX_N], q10n[MAX_N], q11n[MAX_N];
inner_shake256_context sc;

void measure_keygen(size_t n_repetitions, size_t logn) {
	for (size_t i = 0; i < n_repetitions; i++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b);

		for (size_t u = 0; u < MKN(logn); u++) {
			fg_sumsq[logn] += f[u]*f[u] + g[u]*g[u];
			FG_sumsq[logn] += F[u]*F[u] + G[u]*G[u];
		}

		Zf(fft_to_int16)(q00n, q00, logn);
		Zf(fft_to_int16)(q10n, q10, logn);
		Zf(fft_to_int16)(q11n, q11, logn);

		int pk_sz = Zf(encode_pubkey)(NULL, 0, q00n, q10n, logn);
		int sk_sz = Zf(encode_seckey)(NULL, 0, f, g, F, logn);

		pk_sum[logn] += pk_sz;
		pk_sumsq[logn] += pk_sz * pk_sz;
		if (pk_sz < pk_min[logn]) pk_min[logn] = pk_sz;
		if (pk_sz > pk_max[logn]) pk_max[logn] = pk_sz;

		sk_sum[logn] += sk_sz;
		sk_sumsq[logn] += sk_sz * sk_sz;
		if (sk_sz < sk_min[logn]) sk_min[logn] = sk_sz;
		if (sk_sz > sk_max[logn]) sk_max[logn] = sk_sz;

		for  (size_t u = 0; u < MKN(logn); u++) {
			long long x = fpr_rint(q00[u]);
			long long y = fpr_rint(q10[u]);
			long long z = fpr_rint(q11[u]);

			if (u != 0 && u < MKN(logn - 1)) {
				// prevent any biases from self-adjointness of q00 and q11
				q00_sum[logn] += x;
				q00_sumsq[logn] += x*x;
				q11_sum[logn] += z;
				q11_sumsq[logn] += z*z;
			}
			q10_sum[logn] += y;
			q10_sumsq[logn] += y*y;
		}
	}
}

int main() {
	unsigned char seed[48];
	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	size_t nreps[MAX_LOGN + 1];
	for (size_t logn = 1; logn <= MAX_LOGN; logn++) {
		nreps[logn] = 1000;
	}

	for (size_t logn = 1; logn <= MAX_LOGN; logn++) {
		pk_min[logn] = sk_min[logn] = 1000 * 1000;
		measure_keygen(nreps[logn], logn);
		printf(".");
		fflush(stdout);
	}
	printf("\n");

	printf("logn | sigma f,g | sigma F,G | sigma q00 | sigma q10 | sigma q11\n");
	for (size_t logn = 1; logn <= MAX_LOGN; logn++) {
		// printf("logn = %d: sigma_{f,g} ~ %.8f, sigma_{F,G} ~ %.8f, sig_{q00} ~ %.8f, sig_{q10} ~ %.8f, sig_{q11} ~ %.8f\n",
		printf("%4d | %.7f | %.7f | %9.4f | %9.4f | %9.4f\n",
			(int) logn,
			sqrt(fg_sumsq[logn] / (2.0 * nreps[logn] * MKN(logn))),
			sqrt(FG_sumsq[logn] / (2.0 * nreps[logn] * MKN(logn))),
			sqrt(q00_sumsq[logn] / (double) nreps[logn] / (MKN(logn-1) - 1)),
			sqrt(q10_sumsq[logn] / (double) nreps[logn] / MKN(logn)),
			sqrt(q11_sumsq[logn] / (double) nreps[logn] / (MKN(logn-1) - 1))
		);
	}
	for (size_t logn = 1; logn <= MAX_LOGN; logn++) {
		printf("averages Q: %.5f %.5f %.5f\n",
			q00_sum[logn] / (double) nreps[logn] / (MKN(logn-1) - 1),
			q10_sum[logn] / (double) nreps[logn] / MKN(logn),
			q11_sum[logn] / (double) nreps[logn] / (MKN(logn-1) - 1)
		);
		fflush(stdout);
	}

	printf("logn | Average +/- stddev | min  | max \n");
	double avg, std;

	const size_t header_byte = 1;

	printf("-----+- Secret key -------+------+-----\n");
	// Secret key
	for (size_t logn = 1; logn <= MAX_LOGN; logn++) {
		avg = (double) sk_sum[logn] / nreps[logn];
		std = sqrt( (double) sk_sumsq[logn] / nreps[logn] - avg*avg );
		printf("%4d | %7.2f +/- %6.2f | %4lld | %4lld\n",
			(int) logn, header_byte + avg, std, header_byte + sk_min[logn], header_byte + sk_max[logn]);
	}

	printf("-----+- Public key -------+------+-----\n");
	// Public key
	for (size_t logn = 1; logn <= MAX_LOGN; logn++) {
		avg = (double) pk_sum[logn] / nreps[logn];
		std = sqrt( (double) pk_sumsq[logn] / nreps[logn] - avg*avg );
		printf("%4d | %7.2f +/- %6.2f | %4lld | %4lld\n",
			(int) logn, header_byte + avg, std, header_byte + pk_min[logn], header_byte + pk_max[logn]);
	}
	return 0;
}
