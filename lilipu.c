#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// x86_64 specific:
#include<sys/time.h>

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

#include "lilipu64.c"
#include "vrfy.c"


// See nist.c for info on how to use TEMPALLOC
#define TEMPALLOC
// See nist.c
#define NONCELEN   40
// See api.h
#define CRYPTO_SECRETKEYBYTES   2305
#define CRYPTO_PUBLICKEYBYTES   1793
#define CRYPTO_BYTES            1330
#define CRYPTO_ALGNAME          "Lilipu-1024"

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

void to_sage(const char *varname, const int8_t *f, unsigned logn) {
	printf("%s = %d", varname, f[0]);
	for (size_t u = 1; u < MKN(logn); u ++) {
		if (f[u] > 0) printf("+%d*z^%zu", f[u], u);
		if (f[u] < 0) printf("-%d*z^%zu", -f[u], u);
	}
	printf("\n");
}

void to_sage16(const char *varname, const int16_t *f, unsigned logn) {
	printf("%s = %d", varname, f[0]);
	for (size_t u = 1; u < MKN(logn); u ++) {
		if (f[u] > 0) printf("+%d*z^%zu", f[u], u);
		if (f[u] < 0) printf("-%d*z^%zu", -f[u], u);
	}
	printf("\n");
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
// | FUNCTIONS FOR KEY GENERATION                                              |
// =============================================================================
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

	if (!lilipu_solve_NTRU_deepest64(logn, f, g, tmp)) {
		fprintf(stderr, "lilipu_solve_NTRU_deepest64 failed\n");
		return 0;
	}

/*
	printf("Initial normDown(F), normDown(G):\n");
	size_t len = MAX_BL_SMALL64[logn];
	for (u = 0; u < len; u++) printf("%u,", tmp[u]);
	printf("\n");
	for (u = 0; u < len; u++) printf("%u,", tmp[len + u]);
	printf("\n");
*/

	unsigned depth = logn;
	while (depth -- > 0) {
		if (!solve_NTRU_intermediate64(logn, f, g, depth, tmp)) {
			fprintf(stderr, "solve_NTRU_intermediate64(%d) failed\n", depth);
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
		printf("Values for F, G are too big!\n");
		printf("Values(F):");
		for (u = 0; u < n; u++)
			printf(" %d", (int32_t)tmp[u]);
		printf("\n");
		printf("Values(G):");
		for (u = 0; u < n; u++)
			printf(" %d", (int32_t)tmp[u+n]);
		printf("\n");
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
/*
 * Generate a random polynomial with a Gaussian distribution. This function
 * also makes sure that the resultant of the polynomial with phi is odd.
 */
static void
lilipu_poly_small_mkgauss(samplerZ samp, void *samp_ctx, int8_t *f, unsigned logn, fpr isigma_kg, int lim)
{
	size_t n, u;
	int s;
	unsigned mod2;

	n = MKN(logn);
	mod2 = 0;

	for (u = n; u -- > 1; ) {
		do {
			s = samp(samp_ctx, fpr_zero, isigma_kg);
			/*
			 * We need the coefficient to fit within -127..+127;
			 * realistically, this is always the case except for
			 * the very low degrees (N = 2 or 4), for which there
			 * is no real security anyway.
			 */
		} while (s <= -lim || s >= lim);
		mod2 ^= (unsigned)(s & 1);
		f[u] = (int8_t)s;
	}

	do {
		s = samp(samp_ctx, fpr_zero, isigma_kg);
		/*
		 * We need the sum of all coefficients to be 1; otherwise,
		 * the resultant of the polynomial with X^N+1 will be even,
		 * and the binary GCD will fail.
		 */
	} while (s <= -lim || s >= lim || mod2 == (unsigned)(s & 1));
	f[0] = (int8_t)s;
}

// =============================================================================
void
lilipu_keygen(inner_shake256_context *rng,
	int8_t *f, int8_t *g, fpr *q00, fpr *q10, fpr *q11,
	unsigned logn, uint8_t *tmp, fpr isigma_kg, uint8_t tosage)
{
	/*
	 * Algorithm is the following:
	 *
	 *  - Generate f and g with the Gaussian distribution.
	 *
	 *  - If either Res(f,phi) or Res(g,phi) is even, try again.
	 *
	 *  - Solve the NTRU equation fG - gF = 1; if the solving fails,
	 *    try again. Usual failure condition is when Res(f,phi)
	 *    and Res(g,phi) are not prime to each other.
	 */
	size_t n;

	n = MKN(logn);

	/*
	 * In the binary case, coefficients of f and g are generated
	 * independently of each other, with a discrete Gaussian
	 * distribution of standard deviation 1/isigma_kg. Then,
	 * the two vectors have expected norm 2n/isigma_kg.
	 *
	 * We require that Res(f,phi) and Res(g,phi) are both odd (the
	 * NTRU equation solver requires it).
	 */
	for (;;) {
		int8_t *F, *G, *rec_tmp;
		fpr *qxx, *rt1, *rt2, *rt3, *rt4;
		int lim;

		// Normal sampling. We use a fast PRNG seeded from our SHAKE context ('rng').
		sampler_context spc;
		samplerZ samp;
		void *samp_ctx;
		spc.sigma_min = fpr_sigma_min[logn];
		falcon_inner_prng_init(&spc.p, rng);
		samp = Zf(sampler);
		samp_ctx = &spc;

		/*
		 * Verify that all coefficients are within the bounds
		 * defined in max_fg_bits. This is the case with
		 * overwhelming probability; this guarantees that the
		 * key will be encodable with FALCON_COMP_TRIM.
		 */
		lim = 1 << (Zf(max_fg_bits)[logn] - 1);
		lilipu_poly_small_mkgauss(samp, samp_ctx, f, logn, isigma_kg, lim);
		lilipu_poly_small_mkgauss(samp, samp_ctx, g, logn, isigma_kg, lim);
		// when using Falcon sigma (~2.88):
		// poly_small_mkgauss(rng, f, logn);
		// poly_small_mkgauss(rng, g, logn);

		// since we store F, G in tmp, we need 30*1024 instead of 28*1024
		// temporary bytes.
		F = (int8_t *)tmp;
		G = F + n;
		rec_tmp = G + n;

		// Solve the NTRU equation to get F and G.
		lim = (1 << (Zf(max_FG_bits)[logn] - 1)) - 1;
		if (!lilipu_solve_NTRU(logn, F, G, f, g, lim, (uint32_t *)rec_tmp)) {
			continue;
		}

		if (tosage) {
			to_sage("f", f, logn);
			to_sage("g", g, logn);
			to_sage("F", F, logn);
			to_sage("G", G, logn);
		}

		// TODO: use less memory by first calculating q10 and storing
		//       f, g, F, G halfway in the q's.
		// Calculate q00, q10, q11 (in FFT representation) using
		// Q = B * adj(B^{T})
		qxx = (fpr *)rec_tmp;
		rt1 = qxx + n;
		rt2 = rt1 + n;
		rt3 = rt2 + n;
		rt4 = rt3 + n;
		poly_small_to_fp(rt1, f, logn);
		poly_small_to_fp(rt2, g, logn);
		poly_small_to_fp(rt3, F, logn);
		poly_small_to_fp(rt4, G, logn);
		Zf(FFT)(rt1, logn);
		Zf(FFT)(rt2, logn);
		Zf(FFT)(rt3, logn);
		Zf(FFT)(rt4, logn);

		memcpy(q00, rt1, n * sizeof *rt1);
		Zf(poly_mulselfadj_fft)(q00, logn);
		memcpy(qxx, rt2, n * sizeof *rt2);
		Zf(poly_mulselfadj_fft)(qxx, logn);
		Zf(poly_add)(q00, qxx, logn); // q00 = f*bar(f) + g*bar(g)

		memcpy(q10, rt3, n * sizeof *rt3);
		Zf(poly_muladj_fft)(q10, rt1, logn);
		memcpy(qxx, rt4, n * sizeof *rt4);
		Zf(poly_muladj_fft)(qxx, rt2, logn);
		Zf(poly_add)(q10, qxx, logn); // q10 = F*bar(f) + G*bar(g)

		memcpy(q11, rt3, n * sizeof *rt3);
		Zf(poly_mulselfadj_fft)(q11, logn);
		memcpy(qxx, rt4, n * sizeof *rt4);
		Zf(poly_mulselfadj_fft)(qxx, logn);
		Zf(poly_add)(q11, qxx, logn); // q11 = F*bar(F) + G*bar(G)

		/*
		 * Key pair is generated.
		 */
		break;
	}
}

// =============================================================================
// | FUNCTIONS TO CREATE SIGNATURES                                            |
// =============================================================================
// Calculate (f*g) mod (2, phi) and store the result (mod 2) in fg.
// tmp must be of size at least 2n.
static void
lilipu_inner_mulmod2(int8_t *restrict ab, const int8_t *restrict a,
		const uint16_t *restrict b, unsigned logn, uint16_t *tmp)
{
	size_t n, u;

	n = MKN(logn);
	if (logn < 5) {
		size_t v;
		memset(ab, (int8_t)0, n);
		for (u = 0; u < n; u ++)
			for (v = 0; v < n; v ++)
				ab[modp_add(u, v, n)] ^= a[u] & (int8_t)b[v];
		for (u = 0; u < n; u ++)
			ab[u] &= 1;
	} else {
		uint16_t *at, *bt;

		n = MKN(logn);
		at = tmp;
		bt = at + n;

		for (u = 0; u < n; u ++) {
			at[u] = (uint16_t)a[u] & 1;
		}
		for (u = 0; u < n; u ++) {
			bt[u] = (uint16_t)b[u] & 1;
		}

		// Use NTT with q = 12289 from vrfy.c
		mq_NTT(at, logn);
		mq_NTT(bt, logn);
		mq_poly_tomonty(bt, logn);
		mq_poly_montymul_ntt(at, bt, logn);
		mq_iNTT(at, logn);

		// Reduce the output to {0,1} where (Q-x) for 0 < x < Q/2 is cast to x%2.
		for (u = 0; u < n; u ++) {
			ab[u] = (at[u] ^ -(((Q >> 1) - (uint32_t)at[u]) >> 31)) & 1;
		}
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
	const uint16_t *hm, unsigned logn, fpr isigma_sig, uint8_t *restrict tmp)
{
	size_t n, u;
	int8_t *x0, *x1;
	fpr *t0, *t1, *t2;

	n = MKN(logn);
	x0 = (int8_t *)tmp;
	x1 = x0 + n;
	t0 = align_fpr(tmp, x1 + n);
	t1 = t0 + n;
	t2 = t1 + n;

	/*
	 * Set the target vector to [hm, 0] * B (hm is the hashed message).
	 */
	lilipu_inner_mulmod2(x0, f, hm, logn, (uint16_t *)t0);
	lilipu_inner_mulmod2(x1, g, hm, logn, (uint16_t *)t0);

	/*
	 * Apply sampling; result is written over (x0,x1).
	 * Perform Gaussian smoothing to not reveal information on the secret basis.
	 */

	// to_sage("x0", x0, logn);
	// to_sage("x1", x1, logn);

	for (u = 0; u < n; u ++) {
		x0[u] = 2*samp(samp_ctx, fpr_half(fpr_of(x0[u])), isigma_sig) - (x0[u]);
	}
	for (u = 0; u < n; u ++) {
		x1[u] = 2*samp(samp_ctx, fpr_half(fpr_of(x1[u])), isigma_sig) - (x1[u]);
	}

	// to_sage("x0", x0, logn);
	// to_sage("x1", x1, logn);

	/*
	 * Get the signature corresponding to that tiny vector, i.e.
	 * s = x * B^{-1}. Thus s0 = x0 G - x1 F and s1 = -x0 g + x1 f.
	 */
	smallints_to_fpr(t0, f, logn);
	smallints_to_fpr(t1, x1, logn);
	falcon_inner_FFT(t0, logn);
	falcon_inner_FFT(t1, logn);
	falcon_inner_poly_mul_fft(t0, t1, logn);
	smallints_to_fpr(t1, g, logn);
	smallints_to_fpr(t2, x0, logn);
	falcon_inner_FFT(t1, logn);
	falcon_inner_FFT(t2, logn);
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
		s1[u] = (int16_t)fpr_rint(t0[u]);
	}

	// to_sage16("s1", s1, logn);

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
	const uint16_t *hm, unsigned logn, fpr isigma_sig, uint8_t *restrict tmp)
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
// | FUNCTIONS FOR SIGNATURE VERIFICATION                                      |
// =============================================================================
/*
 * Add to polynomial its own adjoint. This function works only in FFT
 * representation.
 */
void
lilipu_inner_poly_addselfadj_fft(fpr *a, unsigned logn)
{
	/*
	 * Since its own conjugate is added to each coefficient,
	 * the result contains only real values.
	 */
	size_t hn, u;

	hn = MKN(logn) >> 1;
	for (u = 0; u < hn; u ++) {
		a[u] = fpr_double(a[u]);
	}
	for (u = 0; u < hn; u ++) {
		a[u + hn] = fpr_zero;
	}
}

// Note: q00, q11 are self adjoint.
// tmp must have size at least 32 * 1024 bytes
int
lilipu_verify(const uint16_t *hm, int16_t *s0, const int16_t *s1,
	const fpr *q00, const fpr *q10, const fpr *q11, unsigned logn, const fpr verif_bound, fpr *tmp)
{
	size_t u, n;
	fpr *t0, *t1, *t2, *t3, trace;

	n = MKN(logn);
	t0 = tmp;
	t1 = t0 + n;
	t2 = t1 + n;
	t3 = t2 + n;

	// if s1 is a valid signature, then s1 == 0 (mod 2)
	for (u = 0; u < n; u ++) {
		if (s1[u] & 1) return 0;
	}

	// Reduce s1 elements modulo q ([0..q-1] range).
	for (u = 0; u < n; u ++) {
		t0[u] = fpr_of(s1[u]);
	}
	for (u = 0; u < n; u ++) {
		t2[u] = fpr_of(hm[u] & 1);
	}

	// Compute s0 = h%2 + 2 round(-q10 s1 / (2 q00))
	falcon_inner_FFT(t0, logn);
	falcon_inner_FFT(t2, logn);
	// copy s1 for later.
	memcpy(t1, t0, n * sizeof *t0);

	// falcon_inner_poly_mulconst(t0, fpr_onehalf, logn);
	falcon_inner_poly_neg(t0, logn);
	falcon_inner_poly_mul_fft(t0, q10, logn); // -q10 s1
	// Note: q00 is self adjoint
	falcon_inner_poly_div_autoadj_fft(t0, q00, logn); // -s1 q10/q00
	falcon_inner_poly_sub(t0, t2, logn); // -s1 q10/q00 - h%2
	falcon_inner_iFFT(t0, logn);

	for (u = 0; u < n; u ++) {
		s0[u] = (hm[u] & 1) + fpr_rint(fpr_half(t0[u]))*2;
	}

	// to_sage16("s0", s0, logn);

	// Currently in memory: s0, s1 (in FFT representation)
	for (u = 0; u < n; u ++) {
		t0[u] = fpr_of(s0[u]);
	}
	falcon_inner_FFT(t0, logn);

	// Currently in memory: s0, s1, s1, s0 (in FFT representation)
	memcpy(t2, t1, n * sizeof *t0);
	memcpy(t3, t0, n * sizeof *t0);

	// Compute s0 q00 s0* + s0 q01 s1* + s1 q10 s0* + s1 q11 s1*
	falcon_inner_poly_mulselfadj_fft(t2, logn);
	falcon_inner_poly_mulselfadj_fft(t3, logn);
	falcon_inner_poly_mul_autoadj_fft(t2, q11, logn); // t2 = s1 q11 s1*
	falcon_inner_poly_mul_autoadj_fft(t3, q00, logn); // t3 = s0 q00 s0*
	falcon_inner_poly_muladj_fft(t1, t0, logn); // t1 = s1 s0*
	falcon_inner_poly_mul_fft(t1, q10, logn); // t1 = s1 q10 s0*

	lilipu_inner_poly_addselfadj_fft(t1, logn); // t1 = s1 q10 s0* + s0 q01 s1*
	falcon_inner_poly_add(t1, t2, logn);
	falcon_inner_poly_add(t1, t3, logn);

	trace = fpr_zero;
	for (u = 0; u < n/2; u ++) {
		trace = fpr_add(trace, t1[u]);
	}

/*
	fpr w = fpr_zero;
	for (u = n/2; u < n; u ++) {
		w = fpr_add(w, t1[u]);
	}
	assert(fpr_lt(w, fpr_inv(fpr_of(1000 * 1000)))); // imag. embeddings < 10^{-6}
*/

	// note: only n/2 embeddings are stored,
	// the others are simply the conjugate embeddings.
	// TODO: this can be optimized in the verif_bound, cancelling with 2 in (2d).
	trace = fpr_double(trace);

	// falcon_inner_iFFT(t1, logn);
	// printf("t1[0] = %.8f\n", t1[0]);
	// printf("value for v: %.8f <=> %.8f\n", v.v, verif_bound.v);
	/*
	 * Signature is valid if and only if `v` is short enough and s1 == 0 (mod 2).
	 */
	return fpr_lt(trace, verif_bound);
}

// =============================================================================
// | TODO: make this into an actually NIST-compatible function like `nist.c`   |
// =============================================================================
int
crypto_lilipu_sign(unsigned char *sm, unsigned long long *smlen,
	const unsigned char *m, unsigned long long mlen,
	const int8_t *restrict f, const int8_t *restrict g, fpr isigma_sig)
{
	TEMPALLOC union {
		uint8_t b[28 * 1024];
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
	Zf(hash_to_point_vartime)(&sc, (uint16_t *)r.hm, 10);

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
// | TESTING CODE                                                              |
// =============================================================================
const size_t logn = 9, n = MKN(logn);

void benchmark_lilipu(fpr isigma_kg, fpr isigma_sig) {
	TEMPALLOC union {
		// uint8_t b[42 * 1024];
		uint8_t b[100 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	TEMPALLOC int8_t f[1024], g[1024];
	TEMPALLOC fpr q00[1024], q10[1024], q11[1024];
	TEMPALLOC uint16_t h[1024];
	TEMPALLOC int16_t sig[1024];
	TEMPALLOC unsigned char seed[48];
	TEMPALLOC inner_shake256_context sc;

	struct timeval t0, t1;
	const int n_repetitions = 1000;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);

	// Generate key pair.
	lilipu_keygen(&sc, f, g, q00, q10, q11, logn, tmp.b, isigma_kg, 0);

	gettimeofday(&t1, NULL);
	printf("Key generation took %lld microseconds\n", time_diff(&t0, &t1));

	// =========================================================================
	// | Benchmark the signing of random messages                              |
	// =========================================================================
	gettimeofday(&t0, NULL);

	for (int rep = 0; rep < n_repetitions; rep++) {
		// make a signature of a random message...
		randombytes((unsigned char *)h, sizeof h);

		// Compute the signature.
		lilipu_sign(sig, &sc, f, g, h, logn, isigma_sig, tmp.b);
	}

	gettimeofday(&t1, NULL);
	double sign_ps = 1000000LL * n_repetitions / (double)time_diff(&t0, &t1);
	printf("Lilipu sign/s = %.1f\n", sign_ps);
}

void benchmark_falcon() {
	TEMPALLOC union {
		uint8_t b[FALCON_KEYGEN_TEMP_10];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	TEMPALLOC union {
		uint8_t b[72 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp2;
	TEMPALLOC union {
		int16_t sig[1024];
		uint16_t hm[1024];
	} r;
	TEMPALLOC int8_t f[1024], g[1024], F[1024], G[1024];
	TEMPALLOC uint16_t h[1024];
	TEMPALLOC unsigned char seed[48];
	TEMPALLOC inner_shake256_context sc;

	struct timeval t0, t1;
	const int n_repetitions = 1000;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);

	// Generate key pair.
	falcon_inner_keygen(&sc, f, g, F, G, h, logn, tmp.b);

	gettimeofday(&t1, NULL);
	printf("Key generation took %lld microseconds\n", time_diff(&t0, &t1));

	// =========================================================================
	// | Start the signing                                                     |
	// =========================================================================
	gettimeofday(&t0, NULL);

	for (int rep = 0; rep < n_repetitions; rep++) {
		// make a signature of a random message...
		randombytes((unsigned char *)h, sizeof h);

		// Compute the signature.
		falcon_inner_sign_dyn(r.sig, &sc, f, g, F, G, r.hm, logn, tmp2.b);
	}

	gettimeofday(&t1, NULL);
	double sign_ps = 1000000LL * n_repetitions / (double)time_diff(&t0, &t1);
	printf("Falcon sign/s = %.1f\n", sign_ps);
}

void test_lilipu_valid_signature(fpr isigma_kg, fpr isigma_sig, fpr verif_bound) {
	TEMPALLOC union {
		uint8_t b[42 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	TEMPALLOC int8_t f[1024], g[1024];
	TEMPALLOC fpr q00[1024], q10[1024], q11[1024];
	TEMPALLOC int16_t sig[1024], s0[1024];
	TEMPALLOC uint16_t h[1024];
	TEMPALLOC unsigned char seed[48];
	TEMPALLOC inner_shake256_context sc;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	// Generate key pair.
	lilipu_keygen(&sc, f, g, q00, q10, q11, logn, tmp.b, isigma_kg, 0);

	for (int rep = 0; rep < 1000; rep++) {
		// make a signature of a random message...
		randombytes((unsigned char *)h, sizeof h);
		// to_sage16("h", (int16_t *)h, logn);

		// Compute the signature.
		lilipu_sign(sig, &sc, f, g, h, logn, isigma_sig, tmp.b);
		// to_sage16("s1", sig, logn);

		assert(lilipu_verify(h, s0, sig, q00, q10, q11, logn, verif_bound, (fpr *)tmp.b));
		// randombytes((unsigned char *)sig, sizeof sig);
		for (size_t u = 0; u < n; u ++)
			sig[u] = 0;
		assert(!lilipu_verify(h, s0, sig, q00, q10, q11, logn, verif_bound, (fpr *)tmp.b));
	}

	printf("Valid signatures were signed.\n");
	printf("No simple forgeries were possible.\n");
}


// try to find a forgery, put answer in s1
int try_forge(int16_t *s0, int16_t *s1, const fpr *q00, const fpr *q10, const fpr *q11,
	const uint16_t *hm, fpr isigma_sig, const fpr verif_bound, fpr *tmp)
{
	size_t u;
	int16_t s0p[1024], s1p[1024];
	fpr *t0, *t1, *t2, *t3, trace;

	t0 = tmp;
	t1 = t0 + n;
	t2 = t1 + n;
	t3 = t2 + n;

	for (u = 0; u < n; u ++)
		s1[u] = 0;

	if (lilipu_verify(hm, s0, s1, q00, q10, q11, logn, verif_bound, tmp))
		return 1; // we are done

	unsigned char seed[48];
	inner_shake256_context sc;
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	sampler_context spc;
	samplerZ samp;
	void *samp_ctx;
	spc.sigma_min = fpr_sigma_min[logn];
	falcon_inner_prng_init(&spc.p, &sc);
	samp = Zf(sampler);
	samp_ctx = &spc;

	fpr sigma = fpr_div(fpr_of(28660196105754623LL), fpr_of(10000000000000000LL));
	while (1) {
		for (u = 0; u < n; u ++)
			s0p[u] = 2*samp(samp_ctx, fpr_half(fpr_of(hm[u] & 1)), sigma) - (hm[u] & 1);
		for (u = 0; u < n; u ++)
			s1p[u] = 2*mkgauss(&sc, logn);

		for (int rep = 10; rep-->0; ) {
			for (u = 0; u < n; u ++)
				s0p[u] += 2*mkgauss(&sc, logn);
			for (u = 0; u < n; u ++)
				s1p[u] += 2*mkgauss(&sc, logn);
		}

		for (u = 0; u < n; u++) t0[u] = fpr_of(s0p[u]);
		falcon_inner_FFT(t0, logn);
		memcpy(t3, t0, n * sizeof *t0);

		for (u = 0; u < n; u++) t1[u] = fpr_of(s1p[u]);
		falcon_inner_FFT(t1, logn);
		memcpy(t2, t1, n * sizeof *t0);
		// Currently in memory: s0, s1, s1, s0 (in FFT representation)

		// Compute s0 q00 s0* + s0 q01 s1* + s1 q10 s0* + s1 q11 s1*
		falcon_inner_poly_mulselfadj_fft(t2, logn);
		falcon_inner_poly_mulselfadj_fft(t3, logn);
		falcon_inner_poly_mul_autoadj_fft(t2, q11, logn); // t2 = s1 q11 s1*
		falcon_inner_poly_mul_autoadj_fft(t3, q00, logn); // t3 = s0 q00 s0*
		falcon_inner_poly_muladj_fft(t1, t0, logn); // t1 = s1 s0*
		falcon_inner_poly_mul_fft(t1, q10, logn); // t1 = s1 q10 s0*

		lilipu_inner_poly_addselfadj_fft(t1, logn); // t1 = s1 q10 s0* + s0 q01 s1*
		falcon_inner_poly_add(t1, t2, logn);
		falcon_inner_poly_add(t1, t3, logn);

		trace = fpr_zero;
		for (u = 0; u < n/2; u ++)
			trace = fpr_add(trace, t1[u]);

		// note: only n/2 embeddings are stored,
		// the others are simply the conjugate embeddings.
		// TODO: this can be optimized in the verif_bound, cancelling with 2 in (2d).
		trace = fpr_double(trace);

		// Signature is valid if and only if `v` is short enough and s2%2 == 0 (<=> s1mod2 = 0).
		if (fpr_lt(trace, verif_bound)) {
			memcpy(s1, s1p, sizeof s1p);
			return 1; // we are done
		}
#ifdef AVX2
		printf("Nope, trace = %.2f\tvs %.2f.\n", trace.v, verif_bound.v);
		printf("s1p = ");
		for (u = 0; u < n; u++) printf("%d,", s1p[u]);
		printf("\n");
#endif
	}
}

void test_forge_signature(fpr isigma_kg, fpr isigma_sig, fpr verif_bound) {
	TEMPALLOC union {
		uint8_t b[42 * 1024];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;
	TEMPALLOC int8_t f[1024], g[1024];
	TEMPALLOC fpr q00[1024], q10[1024], q11[1024];
	TEMPALLOC int16_t sig[1024], s0[1024];
	TEMPALLOC uint16_t h[1024];
	TEMPALLOC unsigned char seed[48];
	TEMPALLOC inner_shake256_context sc;
	int res;


	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	// Generate key pair.
	lilipu_keygen(&sc, f, g, q00, q10, q11, logn, tmp.b, isigma_kg, 0);

	// generate random message
	randombytes((unsigned char *)h, sizeof h);
	// to_sage16("h", (int16_t *)h, logn);

	lilipu_sign(sig, &sc, f, g, h, logn, isigma_sig, tmp.b);
	printf("s1 = ");
	for (size_t u = 0; u < n; u++) printf("%d,", sig[u]);
	printf("\n");
	// to_sage16("s1", sig, logn);

	// try to forge a signature for 'h', and put result in s0
	res = try_forge(s0, sig, q00, q10, q11, h, isigma_sig, verif_bound, (fpr *)tmp.b);
	assert(res);

	assert(!lilipu_verify(h, s0, sig, q00, q10, q11, logn, verif_bound, (fpr *)tmp.b));
	printf("No forgery found.\n");
}

void testmod2() {
	int8_t f[1024], fg[1024], check[1024];
	uint16_t g[1024], tmp[2*1024];

	for (size_t rep = 0; rep < 10; rep ++) {
		for (size_t u = 0; u < n; u ++)
			f[u] = rand();
		for (size_t u = 0; u < n; u ++)
			g[u] = rand();
		for (size_t u = 0; u < n; u ++)
			check[u] = 0;

		for (size_t u = 0; u < n; u ++)
			for (size_t v = 0; v < n; v ++)
				check[modp_add(u, v, n)] ^= f[u] & g[v];
		for (size_t u = 0; u < n; u ++)
			check[u] &= 1;

		lilipu_inner_mulmod2(fg, f, g, logn, tmp);

		for (size_t u = 0; u < n; u ++)
			assert(fg[u] == check[u]);
	}

/*
	struct timeval t0, t1;
	gettimeofday(&t0, NULL);
	for (size_t rep = 0; rep < 100; rep ++) {
		for (size_t u = 0; u < n; u ++) f[u] = rand();
		for (size_t u = 0; u < n; u ++) g[u] = rand();
		lilipu_inner_mulmod2_old(fg, f, g, logn, tmp);
	}
	gettimeofday(&t1, NULL);
	printf("d log d took: %d μs\n", time_diff(&t0, &t1));

	gettimeofday(&t0, NULL);
	for (size_t rep = 0; rep < 100; rep ++) {
		for (size_t u = 0; u < n; u ++) f[u] = rand();
		for (size_t u = 0; u < n; u ++) g[u] = rand();
		lilipu_inner_mulmod2(fg, f, g, logn);
	}
	gettimeofday(&t1, NULL);
	printf("d^2 took: %d μs\n", time_diff(&t0, &t1));
*/
}

int8_t valid_sigma(fpr sigma_sig) {
	return !fpr_lt(sigma_sig, fpr_sigma_min[logn])
		&& fpr_lt(sigma_sig, fpr_div(fpr_of(18205), fpr_of(10000)));
}

int main() {
	unsigned seed = time(NULL);

	const fpr sigma_kg  = fpr_div(fpr_of(1425), fpr_of(1000));
	const fpr sigma_sig = fpr_div(fpr_of(1292), fpr_of(1000));
	// verif_margin = 1 + √(64 * ln(2) / 1024)   (see scheme.sage)
	const fpr verif_margin = fpr_add(fpr_one, fpr_sqrt(fpr_mul(fpr_log2, fpr_div(fpr_of(64), fpr_of(n)))));
	// verif_bound = (verif_margin * 2 * sigma_sig)^2 * (2*d) * d
	// Here, the vector (x0, x1) \in Z^{2d} is sampled from a Discrete Gaussian with sigma equal to 2*sigma_sig
	// and lattice coset (h%2) + 2Z^{2d}, so it has a SQUARED norm of around ~(2sigma_sig)^2 * 2d.
	// Using trace(s Q s^H) = trace(x x^H) = ||x||^2 [K:\QQ] = ||x||^2 d, we arrive at the verif_bound.
	const fpr verif_bound = fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_sqr(fpr_of(n))));

#ifdef __AVX2__
	const fpr sigma_FALCON = fpr_sqrt(fpr_div(fpr_of(117*117*Q), fpr_of(100*100*2*n))); // 1.17 √(q/2n)
	printf("Sigmas: %.2f (lilipu) vs %.2f (falcon)\n", sigma_kg.v, sigma_FALCON.v);
	printf("Verif margin: %.2f, bound: %.2f\n", verif_margin.v, verif_bound.v);
#else
#endif

	// sigmas used in FALCON:
	// for (int i = 1; i <= 10; i++) printf("%d: %.2f\n", i, 1.17 * sqrt(Q / (2 << i)));

	printf("Seed: %u\n", seed);
	srand(seed);
	// prepare_random();
	testmod2();
	assert(valid_sigma(sigma_kg) && valid_sigma(sigma_sig));

	// test_forge_signature(fpr_inv(sigma_kg), fpr_inv(sigma_sig), verif_bound);
	benchmark_lilipu(fpr_inv(sigma_kg), fpr_inv(sigma_sig));
	// benchmark_falcon();

	// test_lilipu_valid_signature(fpr_inv(sigma_kg), fpr_inv(sigma_sig), verif_bound);
	return 0;
}
