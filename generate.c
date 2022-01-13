#include <assert.h>
#include <stdio.h>
// x86_64 specific:
#include<sys/time.h>

#include "inner.h"

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
int valid_sigma(int logn, fpr sigma_sig) {
	return !fpr_lt(sigma_sig, fpr_sigma_min[logn])
		&& fpr_lt(sigma_sig, fpr_div(fpr_of(18205), fpr_of(10000)));
}

void output_poly(int8_t *p, int logn) {
	for (size_t u = 0; u < MKN(logn); u++) {
		if (u) printf(" ");
		printf("%d", p[u]);
	}
	printf("\n");
}

int main(int argc, char **argv) {
	uint8_t b[28 * 512];
	int8_t f[512], g[512], F[512], G[512];
	fpr q00[512], q10[512], q11[512];
	unsigned char shakeseed[48];
	inner_shake256_context sc;

	if (argc != 2) {
		printf("Usage: %s logn\n", argv[0]);
		exit(1);
	}
	size_t logn = atoi(argv[1]);
	if (logn < 1 || logn > 9) {
		printf("Error: %s only works with 1 <= logn <= 9\n", argv[0]);
		exit(1);
	}

	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	srand(1000000 * tv.tv_sec + tv.tv_usec);

	const fpr sigma_kg = fpr_div(fpr_of(1425), fpr_of(1000));
	assert(valid_sigma(logn, sigma_kg));

	// Initialize a RNG.
	randombytes(shakeseed, sizeof shakeseed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, shakeseed, sizeof shakeseed);
	inner_shake256_flip(&sc);

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b, fpr_inv(sigma_kg));

	// Output secret key (f, g, F, G)
	output_poly(f, logn);
	output_poly(g, logn);
	output_poly(F, logn);
	output_poly(G, logn);

	return 0;
}
