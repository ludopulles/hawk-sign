#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// #include "api.h"
// #include "inner.h"
#include "codec.c"
#include "common.c"
#include "fft.c"
#include "fpr.c"
#include "keygen.c"
// #include "nist.c"
#include "rng.c"
#include "shake.c"
#include "sign.c"
#include "vrfy.c"

// See nist.c for info on how to use TEMPALLOC
#define TEMPALLOC
// See nist.c
#define NONCELEN   40
// See api.h
#define CRYPTO_SECRETKEYBYTES   2305
#define CRYPTO_PUBLICKEYBYTES   1793
#define CRYPTO_BYTES            1330
// #define CRYPTO_ALGNAME          "Falcon-1024"

// Perhaps not the way to go:
// void randombytes_init(unsigned char *entropy_input, unsigned char *personalization_string, int security_strength);
int randombytes(unsigned char *x, unsigned long long xlen) {
	// TODO: make proper randomness generator
	while (xlen -- > 0) {
		*x = ((unsigned char) rand());
		++x;
	}
	return 1;
}

// =============================================================================
// | Modified methods from keygen.c, sign.c and vrfy.c                         |
// =============================================================================
// KEY GENERATION:
// - solve_NTRU_deepest
// - solve_NTRU
// - lilipu_keygen
// SIGNING:
// - lilipu_inner_mulmod2
// - lilipu_inner_do_sign
// - lilipu_sign

// =============================================================================
/*
 * Solving the NTRU equation for q = 1, deepest level: compute the resultants
 * of f and g with X^N+1, and use binary GCD. The F and G values are returned
 * in tmp[].
 *
 * Returned value: 1 on success, 0 on error.
 */
static int
lilipu_solve_NTRU_deepest(unsigned logn_top,
	const int8_t *f, const int8_t *g, uint32_t *tmp)
{
	size_t len;
	uint32_t *Fp, *Gp, *fp, *gp, *t1;
	const small_prime *primes;

	len = MAX_BL_SMALL[logn_top];
	primes = PRIMES;

	Fp = tmp;
	Gp = Fp + len;
	fp = Gp + len;
	gp = fp + len;
	t1 = gp + len;

	make_fg(fp, f, g, logn_top, logn_top, 0);

	/*
	 * We use the CRT to rebuild the resultants as big integers.
	 * There are two such big integers. The resultants are always
	 * nonnegative.
	 */
	zint_rebuild_CRT(fp, len, len, 2, primes, 0, t1);

	/*
	 * Apply the binary GCD. The zint_bezout() function works only
	 * if both inputs are odd.
	 *
	 * We can test on the result and return 0 because that would
	 * imply failure of the NTRU solving equation, and the (f,g)
	 * values will be abandoned in that case.
	 */
	if (!zint_bezout(Gp, Fp, fp, gp, len, t1)) {
		return 0;
	}

	return 1;
}

// =============================================================================
/*
 * Solve the NTRU equation, but now for q = 1. Returned value is 1 on success,
 * 0 on error.  G can be NULL, in which case that value is computed but not
 * returned.  If any of the coefficients of F and G exceeds lim (in absolute
 * value), then 0 is returned.
 */
static int
lilipu_solve_NTRU(unsigned logn, int8_t *F, int8_t *G,
	const int8_t *f, const int8_t *g, int lim, uint32_t *tmp)
{
	size_t n, u;
	uint32_t *ft, *gt, *Ft, *Gt, *gm;
	uint32_t p, p0i, r;
	const small_prime *primes;

	n = MKN(logn);

	if (!lilipu_solve_NTRU_deepest(logn, f, g, tmp)) {
		return 0;
	}

	/*
	 * For logn <= 2, we need to use solve_NTRU_intermediate()
	 * directly, because coefficients are a bit too large and
	 * do not fit the hypotheses in solve_NTRU_binary_depth0().
	 */
	if (logn <= 2) {
		unsigned depth;

		depth = logn;
		while (depth -- > 0) {
			if (!solve_NTRU_intermediate(logn, f, g, depth, tmp)) {
				return 0;
			}
		}
	} else {
		unsigned depth;

		depth = logn;
		while (depth -- > 2) {
			if (!solve_NTRU_intermediate(logn, f, g, depth, tmp)) {
				return 0;
			}
		}
		if (!solve_NTRU_binary_depth1(logn, f, g, tmp)) {
			return 0;
		}
		if (!solve_NTRU_binary_depth0(logn, f, g, tmp)) {
			return 0;
		}
	}

	/*
	 * If no buffer has been provided for G, use a temporary one.
	 */
	if (G == NULL) {
		G = (int8_t *)(tmp + 2 * n);
	}

	/*
	 * Final F and G are in fk->tmp, one word per coefficient
	 * (signed value over 31 bits).
	 */
	if (!poly_big_to_small(F, tmp, lim, logn)
		|| !poly_big_to_small(G, tmp + n, lim, logn))
	{
		return 0;
	}

	/*
	 * Verify that the NTRU equation is fulfilled for q = 1. Since all elements
	 * have short lengths, verifying modulo a small prime p works, and allows
	 * using the NTT.
	 *
	 * We put Gt[] first in tmp[], and process it first, so that it does
	 * not overlap with G[] in case we allocated it ourselves.
	 */
	Gt = tmp;
	ft = Gt + n;
	gt = ft + n;
	Ft = gt + n;
	gm = Ft + n;

	primes = PRIMES;
	p = primes[0].p;
	p0i = modp_ninv31(p);
	modp_mkgm2(gm, tmp, logn, primes[0].g, p, p0i);
	for (u = 0; u < n; u ++) {
		Gt[u] = modp_set(G[u], p);
	}
	for (u = 0; u < n; u ++) {
		ft[u] = modp_set(f[u], p);
		gt[u] = modp_set(g[u], p);
		Ft[u] = modp_set(F[u], p);
	}
	modp_NTT2(ft, gm, logn, p, p0i);
	modp_NTT2(gt, gm, logn, p, p0i);
	modp_NTT2(Ft, gm, logn, p, p0i);
	modp_NTT2(Gt, gm, logn, p, p0i);
	// Changed: use q=1
	r = modp_montymul(1, 1, p, p0i);
	for (u = 0; u < n; u ++) {
		uint32_t z;

		z = modp_sub(modp_montymul(ft[u], Gt[u], p, p0i),
			modp_montymul(gt[u], Ft[u], p, p0i), p);
		if (z != r) {
			return 0;
		}
	}

	return 1;
}

// =============================================================================
void
lilipu_keygen(inner_shake256_context *rng,
	int8_t *f, int8_t *g, int8_t *F, int8_t *G,
	unsigned logn, uint8_t *tmp)
{
	/*
	 * Algorithm is the following:
	 *
	 *  - Generate f and g with the Gaussian distribution.
	 *
	 *  - If either Res(f,phi) or Res(g,phi) is even, try again.
	 *
	 *  - If ||(f,g)|| is too large, try again.
	 *
	 *  - If ||B~_{f,g}|| is too large, try again.
	 *
	 *  - If f is not invertible mod phi mod q, try again.
	 *
	 *  - Solve the NTRU equation fG - gF = 1; if the solving fails,
	 *    try again. Usual failure condition is when Res(f,phi)
	 *    and Res(g,phi) are not prime to each other.
	 */
	size_t n, u;
	RNG_CONTEXT *rc;

	n = MKN(logn);
	rc = rng;

	/*
	 * We need to generate f and g randomly, until we find values
	 * such that the norm of (g,-f), and of the orthogonalized
	 * vector, are satisfying. The orthogonalized vector is:
	 *   (q*adj(f)/(f*adj(f)+g*adj(g)), q*adj(g)/(f*adj(f)+g*adj(g)))
	 * (it is actually the (N+1)-th row of the Gram-Schmidt basis).
	 *
	 * In the binary case, coefficients of f and g are generated
	 * independently of each other, with a discrete Gaussian
	 * distribution of standard deviation 1.17*sqrt(q/(2*N)). Then,
	 * the two vectors have expected norm 1.17*sqrt(q), which is
	 * also our acceptance bound: we require both vectors to be no
	 * larger than that (this will be satisfied about 1/4th of the
	 * time, thus we expect sampling new (f,g) about 4 times for that
	 * step).
	 *
	 * We require that Res(f,phi) and Res(g,phi) are both odd (the
	 * NTRU equation solver requires it).
	 */
	for (;;) {
		fpr *rt1, *rt2, *rt3;
		fpr bnorm;
		uint32_t normf, normg, norm;
		int lim;

		/*
		 * The poly_small_mkgauss() function makes sure
		 * that the sum of coefficients is 1 modulo 2
		 * (i.e. the resultant of the polynomial with phi
		 * will be odd).
		 */
		poly_small_mkgauss(rc, f, logn);
		poly_small_mkgauss(rc, g, logn);

		/*
		 * Verify that all coefficients are within the bounds
		 * defined in max_fg_bits. This is the case with
		 * overwhelming probability; this guarantees that the
		 * key will be encodable with FALCON_COMP_TRIM.
		 */
		lim = 1 << (Zf(max_fg_bits)[logn] - 1);
		for (u = 0; u < n; u ++) {
			/*
			 * We can use non-CT tests since on any failure
			 * we will discard f and g.
			 */
			if (f[u] >= lim || f[u] <= -lim
				|| g[u] >= lim || g[u] <= -lim)
			{
				lim = -1;
				break;
			}
		}
		if (lim < 0) {
			continue;
		}

		/*
		 * Bound is 1.17*sqrt(q). We compute the squared
		 * norms. With q = 12289, the squared bound is:
		 *   (1.17^2)* 12289 = 16822.4121
		 * Since f and g are integral, the squared norm
		 * of (g,-f) is an integer.
		 */
		normf = poly_small_sqnorm(f, logn);
		normg = poly_small_sqnorm(g, logn);
		norm = (normf + normg) | -((normf | normg) >> 31);
		if (norm >= 16823) {
			continue;
		}

		/*
		 * We compute the orthogonalized vector norm.
		 */
		rt1 = (fpr *)tmp;
		rt2 = rt1 + n;
		rt3 = rt2 + n;
		poly_small_to_fp(rt1, f, logn);
		poly_small_to_fp(rt2, g, logn);
		Zf(FFT)(rt1, logn);
		Zf(FFT)(rt2, logn);
		Zf(poly_invnorm2_fft)(rt3, rt1, rt2, logn);
		Zf(poly_adj_fft)(rt1, logn);
		Zf(poly_adj_fft)(rt2, logn);
		Zf(poly_mulconst)(rt1, fpr_q, logn);
		Zf(poly_mulconst)(rt2, fpr_q, logn);
		Zf(poly_mul_autoadj_fft)(rt1, rt3, logn);
		Zf(poly_mul_autoadj_fft)(rt2, rt3, logn);
		Zf(iFFT)(rt1, logn);
		Zf(iFFT)(rt2, logn);
		bnorm = fpr_zero;
		for (u = 0; u < n; u ++) {
			bnorm = fpr_add(bnorm, fpr_sqr(rt1[u]));
			bnorm = fpr_add(bnorm, fpr_sqr(rt2[u]));
		}
		if (!fpr_lt(bnorm, fpr_bnorm_max)) {
			continue;
		}

		// Changed: do not calculate h.

		/*
		 * Solve the NTRU equation to get F and G.
		 */
		lim = (1 << (Zf(max_FG_bits)[logn] - 1)) - 1;
		if (!lilipu_solve_NTRU(logn, F, G, f, g, lim, (uint32_t *)tmp)) {
			continue;
		}

		/*
		 * Key pair is generated.
		 */
		break;
	}
}

// =============================================================================
// Calculate (f*g) mod (2, phi) and store the result (mod 2) in f.
// tmp must be of size at least 2n.
static void
lilipu_inner_mulmod2(int8_t *restrict f, const int8_t *restrict g, unsigned logn, uint16_t *tmp)
{
	size_t n, u;
	uint16_t *ft, *gt;

	n = MKN(logn);
	ft = tmp;
	gt = ft + n;

	for (u = 0; u < n; u ++) {
		ft[u] = (uint16_t)f[u] & 1;
	}
	for (u = 0; u < n; u ++) {
		gt[u] = (uint16_t)g[u] & 1;
	}

	/*
	 * Use NTT with q = 12289 from vrfy.c
	 */
	mq_NTT(ft, logn);
	mq_NTT(gt, logn);
	mq_poly_tomonty(gt, logn);
	mq_poly_montymul_ntt(ft, gt, logn);
	mq_iNTT(ft, logn);

	/*
	 * Reduce the output to {0,1} where (Q-x) for 0 < x < Q/2 is cast to x%2.
	 */
	for (u = 0; u < n; u ++) {
		f[u] = (ft[u] ^ -(((Q >> 1) - (uint32_t)ft[u]) >> 31)) & 1;
	}
}

// =============================================================================
/*
 * Compute a signature: the signature contains two vectors, s0 and s1.
 * The s0 vector is not returned. The squared norm of (s0,s1) is
 * computed, and if it is short enough, then s1 is returned into the
 * s1[] buffer, and 1 is returned; otherwise, s1[] is untouched and 0 is
 * returned; the caller should then try again.
 *
 * tmp[] must have room for at least 28 * 2^logn bytes
 */
static int
lilipu_inner_do_sign(samplerZ samp, void *samp_ctx, int16_t *s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const uint16_t *hm, unsigned logn, fpr isigma_sig, int8_t *restrict tmp)
{
	size_t n, u;
	int8_t *x0, *x1;
	fpr *t0, *t1, *t2;

	n = MKN(logn);
	x0 = tmp;
	x1 = tmp + n;
	t0 = align_fpr(tmp, x1 + n);
	t1 = t0 + n;
	t2 = t1 + n;

	for (u = 0; u < n; u ++) {
		x0[u] = hm[u] & 1;
	}
	for (u = 0; u < n; u ++) {
		x1[u] = hm[u] & 1;
	}

	/*
	 * Set the target vector to [hm, 0] * B (hm is the hashed message).
	 */
	lilipu_inner_mulmod2(x0, f, logn, (uint16_t *)t0);
	lilipu_inner_mulmod2(x1, g, logn, (uint16_t *)t0);

	/*
	 * Apply sampling; result is written over (x0,x1).
	 * Perform Gaussian smoothing to not reveal information on the secret basis.
	 */
	for (u = 0; u < n; u ++) {
		int16_t T;
		T = 2*samp(samp_ctx, fpr_scaled(x0[u] & 1, -1), isigma_sig) - (x0[u] & 1);
		assert((T - x0[u]) % 2 == 0); // TODO remove debugging
		x0[u] = T;
	}
	for (u = 0; u < n; u ++) {
		int16_t T;
		T = 2*samp(samp_ctx, fpr_scaled(x1[u] & 1, -1), isigma_sig) - (x1[u] & 1);
		assert((T - x1[u]) % 2 == 0); // TODO remove debugging
		x1[u] = T;
	}

	/*
	 * Get the signature corresponding to that tiny vector, i.e.
	 * s = x * B^{-1}. Thus s0 = x0 G - x1 F and s1 = -x0 g + x1 f.
	 */
	smallints_to_fpr(t0, f, logn);
	smallints_to_fpr(t1, x1, logn);
	falcon_inner_poly_mul_fft(t0, t1, logn);
	smallints_to_fpr(t1, g, logn);
	smallints_to_fpr(t2, x0, logn);
	falcon_inner_poly_mul_fft(t1, t2, logn);
	falcon_inner_poly_sub(t0, t1, logn); // s1 = x1 f - x0 g.

	/*
	 * Extract the signature from t0
	 */

	/*
	 * TODO: check if this signature actually works...
	 * With "normal" degrees (e.g. 512 or 1024), it is very
	 * improbable that the computed vector is not short enough;
	 * however, it may happen in practice for the very reduced
	 * versions (e.g. degree 16 or below). In that case, the caller
	 * will loop, and we must not write anything into s1[] because
	 * s1[] may overlap with the hashed message hm[] and we need
	 * hm[] for the next iteration.
	 */
	falcon_inner_iFFT(t0, logn);
	for (u = 0; u < n; u ++) {
		s1[u] = (int32_t)fpr_rint(t0[u]);
	}
	return 1;
}

// =============================================================================
/*
 * Compute a signature: the signature contains two vectors, s0 and s1.
 * The s0 vector is not returned. The squared norm of (s0,s1) is
 * computed, and if it is short enough, then s1 is returned into the
 * s1[] buffer, and 1 is returned; otherwise, s1[] is untouched and 0 is
 * returned; the caller should then try again.
 *
 * tmp[] must have room for at least nine polynomials.
 */
void
lilipu_sign(int16_t *sig, inner_shake256_context *rng,
	const int8_t *restrict f, const int8_t *restrict g,
	const uint16_t *hm, unsigned logn, fpr isigma_sig, int8_t *restrict tmp)
{
	for (;;) {
		/*
		 * Signature produces short vectors s0 and s1. The
		 * signature is acceptable only if the aggregate vector
		 * s0,s1 is short; we must use the same bound as the
		 * verifier.
		 *
		 * If the signature is acceptable, then we return only s1
		 * (the verifier recomputes s0 from s1, the hashed message,
		 * and the public key).
		 */
		sampler_context spc;
		samplerZ samp;
		void *samp_ctx;

		/*
		 * Normal sampling. We use a fast PRNG seeded from our
		 * SHAKE context ('rng').
		 */
		spc.sigma_min = fpr_sigma_min[logn];
		falcon_inner_prng_init(&spc.p, rng);
		samp = Zf(sampler);
		samp_ctx = &spc;

		/*
		 * Do the actual signature.
		 */
		if (lilipu_inner_do_sign(samp, samp_ctx, sig, f, g, hm, logn,
				isigma_sig, tmp))
			break;
	}
}

// =============================================================================
int
crypto_lilipu_sign(unsigned char *sm, unsigned long long *smlen,
	const unsigned char *m, unsigned long long mlen,
	const int8_t *restrict f, const int8_t *restrict g, fpr isigma_sig)
{
	TEMPALLOC union {
		int8_t b[28 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	TEMPALLOC union {
		int16_t sig[1024];
		uint16_t hm[1024];
	} r;
	TEMPALLOC unsigned char seed[48], nonce[NONCELEN];
	TEMPALLOC unsigned char esig[CRYPTO_BYTES - 2 - sizeof nonce];
	TEMPALLOC inner_shake256_context sc;
	size_t sig_len;

	/*
	 * Create a random nonce (40 bytes).
	 */
	randombytes(nonce, sizeof nonce);

	/*
	 * Hash message nonce + message into a vector.
	 */
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, nonce, sizeof nonce);
	inner_shake256_inject(&sc, m, mlen);
	inner_shake256_flip(&sc);
	Zf(hash_to_point_vartime)(&sc, r.hm, 10);

	/*
	 * Initialize a RNG.
	 */
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);


	/*
	 * Compute the signature.
	 */
	lilipu_sign(r.sig, &sc, f, g, r.hm, 10, isigma_sig, tmp.b);


	/*
	 * Encode the signature and bundle it with the message. Format is:
	 *   signature length     2 bytes, big-endian
	 *   nonce                40 bytes
	 *   message              mlen bytes
	 *   signature            slen bytes
	 */
	esig[0] = 0x20 + 10;
	sig_len = Zf(comp_encode)(esig + 1, (sizeof esig) - 1, r.sig, 10);
	if (sig_len == 0) {
		return -1;
	}
	sig_len ++;
	memmove(sm + 2 + sizeof nonce, m, mlen);
	sm[0] = (unsigned char)(sig_len >> 8);
	sm[1] = (unsigned char)sig_len;
	memcpy(sm + 2, nonce, sizeof nonce);
	memcpy(sm + 2 + (sizeof nonce) + mlen, esig, sig_len);
	*smlen = 2 + (sizeof nonce) + mlen + sig_len;
	return 0;
}

// =============================================================================
// | Code for testing the above functions                                      |
// =============================================================================

void run() {
	const size_t logn = 10, n = MKN(logn);
	const fpr isigma_sig = fpr_inv(fpr_scaled(1331, -10)); // 1.3 ~ 1331 / 1024

	TEMPALLOC union {
		uint8_t b[FALCON_KEYGEN_TEMP_10];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	TEMPALLOC int8_t f[n], g[n], F[n], G[n];
	TEMPALLOC uint16_t h[n];
	TEMPALLOC unsigned char seed[48];
	TEMPALLOC inner_shake256_context rng;

	/*
	 * Generate key pair.
	 */
	randombytes(seed, sizeof seed);
	inner_shake256_init(&rng);
	inner_shake256_inject(&rng, seed, sizeof seed);
	inner_shake256_flip(&rng);
	lilipu_keygen(&rng, f, g, F, G, logn, tmp.b);

	printf("f, g, F, G:\n");
	for (size_t i = 0; i < n; i++) printf("%d,", f[i]);
	printf("\n");
	for (size_t i = 0; i < n; i++) printf("%d,", g[i]);
	printf("\n");
	for (size_t i = 0; i < n; i++) printf("%d,", F[i]);
	printf("\n");
	for (size_t i = 0; i < n; i++) printf("%d,", G[i]);
	printf("\n");

	randombytes((unsigned char *)h, sizeof h);
	// make a signature of a random message...

	TEMPALLOC union {
		int8_t b[28 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp2;
	TEMPALLOC union {
		int16_t sig[1024];
	} r;
	TEMPALLOC inner_shake256_context sc;

	/*
	 * Initialize a RNG.
	 */
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	/*
	 * Compute the signature.
	 */
	lilipu_sign(r.sig, &sc, f, g, h, logn, isigma_sig, tmp2.b);
	for (size_t i = 0; i < n; i++)
		printf("%d,", r.sig[i]);
	printf("\n");
}

void testmod2() {
	const size_t logn = 10, n = MKN(logn);
	int8_t f[n], g[n], check[n];
	uint16_t tmp[2*n];

	for (size_t i=0; i < n; i++)
		f[i] = rand();
	for (size_t i=0; i < n; i++)
		g[i] = rand();
	
	for (size_t i=0; i < n; i++)
		check[i] = 0;
	for (size_t i=0; i < n; i++)
		for (size_t j=0; j < n; j++)
			check[(i+j)%n] ^= (f[i]*g[j])&1;

	lilipu_inner_mulmod2(f, g, logn, tmp);

	for (size_t i=0; i < n; i++)
		assert( f[i] == check[i] );
}

int main() {
	srand(time(NULL));

	// fpr x = fpr_inv(fpr_scaled(1331, -10));
	// long long N = fpr_rint(fpr_mul(x, fpr_of(1000000)));
	// printf("%lld.%6lld\n", N / 1000000, N % 1000000);

	run();
	// testmod2();

	return 0;
}
