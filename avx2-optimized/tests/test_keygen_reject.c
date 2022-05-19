/*
 * Estimates the probability that keygen fails because cst(1 / q00) is too large.
 */
#include <assert.h>
#include <stdio.h>

#include "../keygen.c"

int num_samples = 0, num_fails = 0;

int num_samples_fg = 0, num_fails_fg = 0;

void check_fg_condition(inner_shake256_context *rng,
	int8_t *restrict f, int8_t *restrict g, // secret key
	fpr *restrict r0, fpr *restrict r1, fpr *restrict r2, // public key
	unsigned logn)
{
	prng p;
	Zf(prng_init)(&p, rng);

	poly_small_mkgauss(&p, f, logn);
	poly_small_mkgauss(&p, g, logn);

	Zf(int8_to_fft)(r1, f, logn);
	Zf(int8_to_fft)(r2, g, logn);
	Zf(poly_invnorm2_fft)(r0, r1, r2, logn);
	Zf(iFFT)(r0, logn); // 1 / Q_00

	/*
	 * For n = 512, we reject a key pair if cst(1 / Q_{00}) > 0.001,
	 * as the failure probability of decompressing a signature is bounded
	 * from above by 9.98263 10^{-32} < 2^{-64}.
	 *
	 * Note that experimentally, cst(1 / Q_{00}) ~ 0.00097 +/- 0.00016, so
	 * rejection happens regularly because of this criterion.
	 */
	num_samples_fg++;
	if (logn == 9 && fpr_lt(fpr_inv(fpr_of(1000)), r0[0]))
		num_fails_fg++;
}

/* see inner.h for keygen */
void
keygen_count_fails(inner_shake256_context *rng,
	int8_t *restrict f, int8_t *restrict g, // secret key
	int8_t *restrict F, int8_t *restrict G, // secret key
	fpr *restrict q00, fpr *restrict q10, fpr *restrict q11, // public key
	unsigned logn, uint8_t *restrict tmp)
{
	/*
	 * Algorithm is the following:
	 *
	 *  - Generate f and g with the Gaussian distribution.
	 *
	 *  - If either N(f) or N(g) is even, try again.
	 *
	 *  - Solve the NTRU equation fG - gF = 1; if the solving fails, try again.
	 *    Usual failure condition is when N(f) and N(g) are not coprime.
	 *
	 *  - Use Babai Reduction on F, G.
	 *
	 *  - Calculate the Gram matrix of the basis [[f, g], [F, G]].
	 */
	prng p;

	for (;;) {
		/*
		 * The coefficients of polynomials f and g will be generated from a
		 * discrete gaussian that draws random numbers from a fast PRNG that is
		 * seeded from a SHAKE context ('rng').
		 */
		Zf(prng_init)(&p, rng);

		/*
		 * The coefficients of f and g are generated independently of each other,
		 * with a discrete Gaussian distribution of standard deviation 1.500. The
		 * expected l2-norm of (f, g) is 2n 1.500^2.
		 *
		 * We require that N(f) and N(g) are both odd (the NTRU equation solver
		 * requires it).
		 */
		poly_small_mkgauss(&p, f, logn);
		poly_small_mkgauss(&p, g, logn);

		Zf(int8_to_fft)(q10, f, logn);
		Zf(int8_to_fft)(q11, g, logn);
		Zf(poly_invnorm2_fft)(q00, q10, q11, logn);
		Zf(iFFT)(q00, logn); // 1 / Q_00

		/*
		 * For n = 512, we reject a key pair if cst(1 / Q_{00}) > 0.001,
		 * as the failure probability of decompressing a signature is bounded
		 * from above by 9.98263 10^{-32} < 2^{-64}.
		 *
		 * Note that experimentally, cst(1 / Q_{00}) ~ 0.00097 +/- 0.00016, so
		 * rejection happens regularly because of this criterion.
		 */
		num_samples++;
		if (logn == 9 && fpr_lt(fpr_inv(fpr_of(1000)), q00[0])) {
			num_fails++;
			continue;
		}

		/*
		 * Try to solve the NTRU equation for polynomials f and g, i.e. find
		 * polynomials F, G that satisfy
		 *
		 *     f * G - g * F = 1 (mod X^n + 1).
		 */
		if (Zf(complete_private)(f, g, F, G, q00, q10, q11, logn, tmp)) {
			/*
			 * A valid key pair is generated.
			 */
			return;
		}
	}
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
const size_t logn = 9, n = MKN(logn);

void measure_keygen_fails() {
	uint8_t b[48 << logn];
	int8_t f[n], g[n], F[n], G[n];
	fpr q00[n], q10[n], q11[n];
	unsigned char seed[48];
	inner_shake256_context sc;

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	for (int i = 0; i < 100 * 1000; i++) {
		// Generate key pair.
		check_fg_condition(&sc, f, g, q00, q10, q11, logn);
	}
	printf("Keygen sampled %d times and failed %d times, so probability on cst(1/Q00) > 0.001 is %.6f\n",
		num_samples_fg, num_fails_fg, (double)num_fails_fg / num_samples_fg);

	for (int i = 0; i < 1000; i++) {
		// Generate key pair.
		keygen_count_fails(&sc, f, g, F, G, q00, q10, q11, logn, b);
	}
	printf("Keygen sampled %d times and failed %d times, so probability on cst(1/Q00) > 0.001 is %.6f\n",
		num_samples, num_fails, (double)num_fails / num_samples);

}

int main() {
	measure_keygen_fails();
	return 0;
}
