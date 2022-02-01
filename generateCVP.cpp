// g++ build/ffo.o build/fft.o build/fpr.o build/keygen.o build/rng.o build/sampler.o build/shake.o build/sign.o build/vrfy.o generateCVP.cpp
#include <assert.h>
#include <stdio.h>
// x86_64 specific:
#include <sys/time.h>

#include <fplll.h>
#include <vector>

extern "C" {
	#define restrict
	#include "inner.h"
}

using namespace std;
using namespace fplll;


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
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, fpr_inv(sigma_kg), logn, b);

	// Output secret key (f, g, F, G)
	output_poly(f, logn);
	output_poly(g, logn);
	output_poly(F, logn);
	output_poly(G, logn);


	size_t n = MKN(logn);
	ZZ_mat<mpz_t> m(n, 2*n);
	vector<Z_NR<mpz_t>> target(2*n), sol_coord(2*n);

	for (size_t u = 0; u < n; u ++) {
		// (f,g) x^u
		for (size_t v = 0; v < u; v++) {
			m(u, v) = -f[n-u+v];
			m(u, n + v) = -g[n-u+v];
		}
		for (size_t v = u; v < n; v++) {
			m(u, v) = f[v-u];
			m(u, n + v) = g[v-u];
		}
	}
	for (size_t v = 0; v < n; v++) {
		target[v] = F[v];
		target[n + v] = G[v];
	}

const char *const RED_STATUS_STR[RED_STATUS_MAX] = {
    "success",
    "",
    "infinite number in GSO",
    "infinite loop in babai",
    "infinite loop in LLL",
    "error in SVP solver",
    "error in BKZ",
    "time limit exceeded in BKZ",
    "loops limit exceeded in BKZ",
    "error in HLLL",
    "increase of the norm",
    "error in weak size reduction",
    "Please see https://github.com/fplll/fplll/wiki/fplll-errors-FAQ for more information."
};


	int status;
	status = lll_reduction(m, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER, FT_DEFAULT, 0, LLL_DEFAULT);
	cerr << "Status LLL: " << RED_STATUS_STR[status] << endl;
	if (status != RED_SUCCESS)
		return status;

	status = fplll::closest_vector(m, target, sol_coord, CVPM_FAST, CVP_VERBOSE);
	cerr << "Status CVP: " << RED_STATUS_STR[status] << endl;
	if (status != RED_SUCCESS)
		return status;

	for (size_t v = 0; v < n; v++) {
		cout << "(" << sol_coord[v] << "," << sol_coord[n+v] << ")";
	}
	cout << endl;

	int8_t newF[512], newG[512];
	for (size_t v = 0; v < n; v++) {
		newF[v] = F[v] - (int) sol_coord[v].get_si();
		newG[v] = G[v] - (int) sol_coord[n + v].get_si();
	}

	output_poly(newF, logn);
	output_poly(newG, logn);

	return 0;
}
