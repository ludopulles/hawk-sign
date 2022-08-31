/*
 * Estimates the probability that keygen fails and has to start over due to
 * various reasons.
 *
 * Moreover, it prints the size of the compressed public key, including the
 * size when using Huffman encoding compared to the compressed-gaussian
 * technique.
 */

#include <cassert>
#include <climits>
#include <cstdio>
#include <vector>
#include <mutex>
#include <thread>

// x86_64 specific:
#include <sys/time.h>

extern "C" {
	#define restrict
	#include "../keygen.c"
}

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

void poly_output(fpr *p, size_t logn) {
	for (size_t u = 0; u < MKN(logn); u++) {
		if (u) printf(" ");
		printf("%ld", fpr_rint(p[u]));
	}
	printf("\n");
}

/* see inner.h for keygen */
void
keygen_count_fails(inner_shake256_context *sc,
	int8_t *restrict f, int8_t *restrict g,
	int8_t *restrict F, int8_t *restrict G,
	int16_t *restrict iq00, int16_t *restrict iq01,
	unsigned logn, uint8_t *restrict tmp,
	long long *inv_fails, long long *normfg_fails, long long *norminvq00_fails,
	long long *NTRU_fails, long long *coeff_fails)
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
	size_t n, hn, u;
	uint8_t fg_okay;
	int16_t *iq11, *atmp;
	int32_t norm, bound, x;
	uint32_t p, p0i, *gm, *igm, *rf1, *rf2;
	prng rng;
	fpr *rt1, *rt2, *rt3;

	n = MKN(logn);
	hn = n >> 1;

	iq11 = (int16_t *)tmp;
	atmp = iq11 + n;

	rt1 = (fpr *)tmp;
	rt2 = rt1 + n;
	rt3 = rt2 + n;

	gm = (uint32_t *)tmp;
	igm = gm + n;
	rf1 = igm + n;
	rf2 = rf1 + n;

	p = PRIMES[0].p;
	p0i = modp_ninv31(p);

	for (;;) {
		/*
		 * The coefficients of polynomials f and g will be generated from a
		 * discrete gaussian that draws random numbers from a fast PRNG that is
		 * seeded from a SHAKE context ('sc').
		 */
		Zf(prng_init)(&rng, sc);

		/*
		 * The coefficients of f and g are generated independently of each
		 * other, with a discrete Gaussian distribution of standard deviation
		 * 1.500. The expected l2-norm of (f, g) is 2n sigma^2.
		 *
		 * We require that N(f) and N(g) are both odd (the binary GCD in the
		 * NTRU solver requires it), so we require (fg_okay & 1) == 1.
		 */
		fg_okay = poly_small_mkgauss(&rng, f, logn)
			& poly_small_mkgauss(&rng, g, logn) & 1u;
		if (fg_okay == 0U) { continue; }

		fg_okay &= Zf(mf_is_invertible)(f, logn, tmp);
		if (fg_okay == 0U) { (*inv_fails)++; continue; }

		Zf(int8_to_fft)(rt2, f, logn);
		Zf(int8_to_fft)(rt3, g, logn);
		Zf(poly_invnorm2_fft)(rt1, rt2, rt3, logn);
		Zf(iFFT)(rt1, logn);

		if (logn == 9) {
			/*
			 * For n = 512, we reject a key pair if cst(1/q00) >= 1/1000, as
			 * the failure probability of decompressing a signature is bounded
			 * from above by 1.9e-32 < 2^{-105}. Experimentally this fails
			 * with a probability of 9%.
			 */
			fg_okay &= fpr_lt(rt1[0], fpr_inv(fpr_of(1000)));
		} else if (logn == 10) {
			/*
			 * For n = 1024, we reject a key pair if cst(1/q00) >= 1/3000, as
			 * the failure probability of decompressing a signature is bounded
			 * from above by 1.2e-95 < 2^{-315}. Experimentally this fails
			 * with a probability of 0.9%.
			 */
			fg_okay &= fpr_lt(rt1[0], fpr_inv(fpr_of(3000)));
		}

		if (fg_okay == 0U) { (*norminvq00_fails)++; continue; }

		/*
		 * If NTT_p(q00) has a zero, we cannot invert it in
		 * Zf(uncompressed_verify_NTT), so remove this key.
		 */
		modp_mkgm2(gm, igm, logn, PRIMES[0].g, p, p0i);
		for (u = 0; u < n; u++) {
			rf1[u] = modp_set((int32_t)f[u], p);
			rf2[u] = modp_set((int32_t)g[u], p);
		}
		modp_NTT2(rf1, gm, logn, p, p0i);
		modp_NTT2(rf2, gm, logn, p, p0i);
		for (u = 0; u < hn; u++) {
			uint32_t prod;

			prod = modp_add(modp_montymul(rf1[u], rf1[n - 1 - u], p, p0i),
							modp_montymul(rf2[u], rf2[n - 1 - u], p, p0i), p);
			fg_okay &= (1U - prod) >> 31;
		}
		if (fg_okay == 0U) { (*inv_fails)++; continue; }

		/*
		 * If the l2-norm of (f, g) is shorter than sigma_sec^2 * 2n, BKZ may
		 * return a shortest vector when given the public key much faster than
		 * other instances, so this secret key is not secure to use.
		 * Thus, set fg_okay to 0 when ||(f, g)||^2 < Zf(l2bound)[logn]/4.
		 *
		 * For NIST-5, ssec != sver so use a different bound array.
		 */
		norm = 0;
		for (u = 0; u < n; u++) {
			norm += (int32_t)f[u] * (int32_t)f[u];
			norm += (int32_t)g[u] * (int32_t)g[u];
		}

		if (logn == 10) {
			norm -= l2bound_ssec_1024[logn];
		} else {
			norm -= (int32_t)(Zf(l2bound_512)[logn] >> 2);
		}

		fg_okay &= ((uint32_t) -norm) >> 31;
		if (fg_okay == 0U) { (*normfg_fails)++; continue; }

		if (fg_okay == 0) {
			assert(0);
			/*
			 * Generation of (f, g) failed because:
			 * 1) N(f) or N(g) was even,
			 * 2) NTT(f) had zero coefficient,
			 * 3) NTT(q00) had zero coefficient,
			 * 4) cst(q00) = ||(f,g)||^2 < l2bound(logn)/4, or
			 * 5) cst(1/q00) was too large.
			 *
			 * Thus, resample f and g.
			 */
			continue;
		}

		/*
		 * Try to solve the NTRU equation for polynomials f and g, i.e. find
		 * polynomials F, G that satisfy
		 *
		 *     f * G - g * F = 1 (mod X^n + 1).
		 */
		if (!solve_NTRU(logn, F, G, f, g, 127, (uint32_t *)tmp)) {
			(*NTRU_fails)++;
			continue;
		}

		Zf(make_public)(f, g, F, G, iq00, iq01, iq11, logn, (uint8_t *)atmp);

		/*
		 * Check the bounds on q00 and q11.
		 */
		bound = (int32_t)(1U << Zf(bits_q00)[logn]);
		for (u = 1; u < hn; u++) {
			x = iq00[u];
			fg_okay &= (+x - bound) >> 31;
			fg_okay &= (-x - bound) >> 31;
			fg_okay &= x == -iq00[n - u];
		}

		bound = (int32_t)(1U << Zf(bits_q11)[logn]);
		for (u = 1; u < hn; u++) {
			x = iq11[u];
			fg_okay &= (+x - bound) >> 31;
			fg_okay &= (-x - bound) >> 31;
			fg_okay &= x == -iq11[n - u];
		}

		bound = (int32_t)(1U << Zf(bits_q01)[logn]);
		for (u = 0; u < n; u++) {
			x = iq01[u];
			fg_okay &= (+x - bound) >> 31;
			fg_okay &= (-x - bound) >> 31;
		}

		if (fg_okay == 0) {
			(*coeff_fails)++;
			/*
			 * There was a coefficient that was too large, i.e. there
			 * was an x in the above with
			 *   x <= -bound or x >= bound,
			 * or q00 and q11 failed to be selfadjoint (should not happen).
			 */
			continue;
		}

		/*
		 * A valid key pair is generated.
		 */
		break;
	}
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
#define LOGN (10)
#define N (MKN(LOGN))

struct WorkerResult {
	long long iters, time_diff,
		inv_fails, normfg_fails, norminvq00_fails, NTRU_fails, coeff_fails;

	WorkerResult() : iters(0), time_diff(0),
		inv_fails(0), normfg_fails(0), norminvq00_fails(0), NTRU_fails(0), coeff_fails(0) {}

	void combine(const WorkerResult &res) {
		iters += res.iters;
		time_diff += res.time_diff;
		inv_fails += res.inv_fails;
		normfg_fails += res.normfg_fails;
		norminvq00_fails += res.norminvq00_fails;
		NTRU_fails += res.NTRU_fails;
		coeff_fails += res.coeff_fails;
	}
};

WorkerResult measure_keygen(unsigned logn) {
	uint8_t b[48 << LOGN];
	int8_t f[N], g[N], F[N], G[N];
	int16_t iq00[N], iq01[N];
	unsigned char seed[48];
	inner_shake256_context sc;

	struct timeval t0, t1;

	WorkerResult res;
	res.iters = logn == 9 ? 400 : 100;

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);
	for (int i = 0; i < res.iters; i++) {
		// Generate key pair.
		keygen_count_fails(&sc, f, g, F, G, iq00, iq01, logn, b,
			&res.inv_fails, &res.normfg_fails, &res.norminvq00_fails,
			&res.NTRU_fails, &res.coeff_fails);
	}
	gettimeofday(&t1, NULL);

	res.time_diff = time_diff(&t0, &t1);
	return res;
}

WorkerResult tot;
std::mutex mx;

void work(unsigned logn)
{
	WorkerResult result = measure_keygen(logn);

	/* acquire mutex lock */ {
		std::lock_guard<std::mutex> guard(mx);
		tot.combine(result);
	}
}

void measure_keygen_multithreaded(unsigned logn)
{
	tot = WorkerResult();
	const int nthreads = 4;
	std::thread* pool[nthreads-1];
	for (int i = 0; i < nthreads-1; i++) pool[i] = new std::thread(work, logn);
	work(logn);
	for (int i = 0; i < nthreads-1; i++) pool[i]->join(), delete pool[i];

	// Collect results
	double kg_duration = ((double) tot.time_diff) / tot.iters;
	printf("Average time per keygen HAWK-%zu: %.3f ms (%lld samples)\n", MKN(logn), kg_duration / 1000.0, tot.iters);

	/*
	 * Crunch analysis on what may fail during basis completion, once f, g are
	 * generated.
	 */
	int samples = tot.iters + tot.inv_fails + tot.normfg_fails +
		tot.norminvq00_fails + tot.NTRU_fails + tot.coeff_fails;
	printf("Pr[ NTT(f) or NTT(q00) has 0 ] = %.2f%%\n", 100.0 * tot.inv_fails / samples);
	printf("Pr[ || (f,g) ||^2 too small  ] = %.2f%%\n", 100.0 * tot.normfg_fails / samples);
	printf("Pr[ cst(1/q00) too large     ] = %.2f%%\n", 100.0 * tot.norminvq00_fails / samples);
	printf("Pr[ NTRU_solve fails         ] = %.2f%%\n", 100.0 * tot.NTRU_fails / samples);
	printf("Pr[ coeff too large          ] = %.2f%%\n", 100.0 * tot.coeff_fails / samples);
	printf("Pr[ keygen works             ] = %.2f%%\n", 100.0 * tot.iters / samples);
}

void report_invq00_fail_prob(unsigned logn) {
	const int n_repetitions = 100 * 1000;

	int8_t f[N], g[N];
	fpr rt1[N], rt2[N], rt3[N], avg_invq00 = fpr_zero;
	unsigned char seed[48];
	inner_shake256_context sc;
	prng p;

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);
	Zf(prng_init)(&p, &sc);

	int inv_nu = logn == 10 ? 3000 : 1000;

	int num_fails = 0;
	for (int i = 0; i < n_repetitions; i++) {
		poly_small_mkgauss(&p, f, logn);
		poly_small_mkgauss(&p, g, logn);

		Zf(int8_to_fft)(rt2, f, logn);
		Zf(int8_to_fft)(rt3, g, logn);
		Zf(poly_invnorm2_fft)(rt1, rt2, rt3, logn);
		Zf(iFFT)(rt1, logn);

		num_fails += !fpr_lt(rt1[0], fpr_inv(fpr_of(inv_nu)));
		avg_invq00 = fpr_add(avg_invq00, rt1[0]);
	}
	avg_invq00 = fpr_div(avg_invq00, fpr_of(n_repetitions));

	printf("\nAverage value cst(1/q00) = %.6f\n", *(double*)&avg_invq00);
	printf("Experimental probability that cst(1/q00) >= %.6f is %.6f\n",
		((double) 1.0) / inv_nu,
		((double) num_fails) / n_repetitions);
}


int main() {
	measure_keygen_multithreaded( 9);
	measure_keygen_multithreaded(10);

	report_invq00_fail_prob( 9);
	report_invq00_fail_prob(10);
	return 0;
}
