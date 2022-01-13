#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// x86_64 specific:
#include<sys/time.h>

#include "codec.c"
#include "common.c"
#include "fft.c"
#include "fpr.c"
#include "rng.c"
#include "shake.c"
#include "lilipu_keygen.c"

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
const size_t logn = 9, n = MKN(logn);

int valid_sigma(fpr sigma_sig) {
	return !fpr_lt(sigma_sig, fpr_sigma_min[logn])
		&& fpr_lt(sigma_sig, fpr_div(fpr_of(18205), fpr_of(10000)));
}

void output_poly(int8_t *p) {
	for (size_t u = 0; u < n; u++) {
		if (u) printf(" ");
		printf("%d", p[u]);
	}
	printf("\n");
}

int main() {
	uint8_t b[28 * 512];
	int8_t f[512], g[512], F[512], G[512];
	fpr q00[512], q10[512], q11[512];
	unsigned char shakeseed[48];
	inner_shake256_context sc;

	// Output secret key (f, g, F, G)
	unsigned seed = time(NULL);
	srand(seed);

	const fpr sigma_kg = fpr_div(fpr_of(1425), fpr_of(1000));
	assert(valid_sigma(sigma_kg));

	// Initialize a RNG.
	randombytes(shakeseed, sizeof shakeseed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, shakeseed, sizeof shakeseed);
	inner_shake256_flip(&sc);

	// Generate key pair.
	lilipu_keygen(&sc, f, g, F, G, q00, q10, q11, logn, b, fpr_inv(sigma_kg));

	// Output secret key (f, g, F, G)
	output_poly(f);
	output_poly(g);
	output_poly(F);
	output_poly(G);

	return 0;
}
