/*
 * Hawk signature verification.
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2017-2019  Falcon Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ===========================(LICENSE END)=============================
 *
 * @author   Thomas Pornin <thomas.pornin@nccgroup.com>
 */

#include <assert.h>
#include <math.h>
#include "inner.h"

/*
 * The hash of a message has two parts: (h0, h1), where each is a polynomial of
 * length n with coefficients in {0,1}. However, these bits are collected into
 * bytes for saving RAM usage. This function returns the index in the hash
 * array at which the bits start for h1.
 */
#define SECOND_HASH(h, logn) \
	((h) + ((logn) <= 3 ? 1u : 1u << ((logn) - 3)))

/*
 * If s != NULL, set the polynomial p equal to h - 2 * s.
 * Otherwise, set p equal to h.
 * Returns the polynomial in FFT format.
 */
static void
hash_to_fft(fpr *p, const uint8_t *h, const int16_t *s, unsigned logn)
{
	size_t n, u, v;
	uint8_t hash;

	n = MKN(logn);
	if (logn <= 3) {
		for (u = 0; u < n; u ++) {
			p[u] = fpr_of(((h[0] >> u) & 1) - (s == NULL ? 0 : 2 * s[u]));
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				p[u] = fpr_of((hash & 1) - (s == NULL ? 0 : 2 * s[u]));
				hash >>= 1;
			}
		}
	}
	Zf(FFT)(p, logn);
}

/*
 * Returns whether squared geometric norm of (t0, t1) with respect to Q is at
 * most the specified allowed l2-bound. The t0, t1 are assumed to be in FFT
 * representation.
 */
static int
has_short_trace(const fpr *restrict t0, const fpr *restrict t1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	unsigned logn)
{
	size_t u, hn;
	fpr trace;

	hn = MKN(logn) >> 1;
	trace = fpr_zero;

	/*
	 * Calculate trace( (t0, t1)^* Q (t0, t1) ) / n, and determine if this is
	 * short. The trace of a polynomial in FFT representation is the sum of its
	 * `n` complex embeddings, and since half of the embeddings are stored, we
	 * only need to sum the real value for each pair of complex embeddings.
	 * 
	 * Add the contribution of t0 q00 t0*.
	 */
	for (u = 0; u < hn; u ++) {
		fpr norm_t0;

		norm_t0 = fpr_add(fpr_sqr(t0[u]), fpr_sqr(t0[u + hn]));
		trace = fpr_add(trace, fpr_mul(norm_t0, q00[u]));
	}

	/*
	 * Add the contribution of t1 q11 t1*.
	 */
	for (u = 0; u < hn; u ++) {
		fpr norm_t1;

		norm_t1 = fpr_add(fpr_sqr(t1[u]), fpr_sqr(t1[u + hn]));
		trace = fpr_add(trace, fpr_mul(norm_t1, q11[u]));
	}

	/*
	 * Add the contribution of t1 q10 t0* + t0 q01 t1*.
	 */
	for (u = 0; u < hn; u ++) {
		fpr re, im;

		// t0* t1
		re = fpr_add(fpr_mul(t0[u], t1[u]), fpr_mul(t0[u + hn], t1[u + hn]));
		im = fpr_sub(fpr_mul(t0[u], t1[u + hn]), fpr_mul(t0[u + hn], t1[u]));

		trace = fpr_add(trace, fpr_double(
			fpr_sub(fpr_mul(re, q10[u]), fpr_mul(im, q10[u + hn]))
		));
	}

	/*
	 * Renormalize to obtain the squared geometric norm of (t0, t1) w.r.t Q.
	 */
	trace = fpr_div(trace, fpr_of(hn));

	/*
	 * Check whether the norm is in range [0, 2^31). Signature is valid iff
	 * squared norm of (t0, t1) w.r.t. Q is at most bound.
	 */
	return !fpr_lt(trace, fpr_zero)
		&& fpr_lt(trace, fpr_ptwo31m1)
		&& (uint32_t)fpr_rint(trace) <= L2BOUND(logn);
}

/* =================================================================== */

/* see inner.h */
void
Zf(complete_pubkey)(const int16_t *restrict iq00, const int16_t *restrict iq10,
	fpr *restrict q00, fpr *restrict q10, fpr *restrict q11, unsigned logn)
{
	size_t u, hn;

	hn = MKN(logn - 1);

	/*
	 * Doing this in reverse, allows iq00, iq10 to overlap with the begin of
	 * q00.
	 */
	Zf(int16_to_fft)(q10, iq10, logn);
	Zf(int16_to_fft)(q00, iq00, logn);

	/*
	 * Reconstruct q11 using q11 = (1 + q10 adj(q10)) / q00.
	 */
	Zf(poly_prod_selfadj_fft)(q11, q10, logn);
	for (u = 0; u < hn; u ++) {
		q11[u] = fpr_add(q11[u], fpr_one);
	}
	Zf(poly_div_autoadj_fft)(q11, q00, logn);
}

/* see inner.h */
int
Zf(uncompressed_verify)(const uint8_t *restrict h,
	const int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n;
	fpr *t0, *t1;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;

	/*
	 * Put (t0, t1) = (h0 - 2 s0, h1 - 2 s1) in FFT representation.
	 */
	hash_to_fft(t0, h, s0, logn);
	hash_to_fft(t1, SECOND_HASH(h, logn), s1, logn);

	return Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn)
		&& has_short_trace(t0, t1, q00, q10, q11, logn);
}

/* see inner.h */
int
Zf(recover_and_verify)(const uint8_t *restrict h,
	int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n, u, v, w;
	uint8_t h0;
	fpr *t0, *t1;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;

	hash_to_fft(t1, SECOND_HASH(h, logn), s1, logn);
	memcpy(t0, t1, n * sizeof *t1);
	Zf(poly_mul_fft)(t0, q10, logn);
	Zf(poly_div_autoadj_fft)(t0, q00, logn);
	Zf(iFFT)(t0, logn);

	/*
	 * Recover s0 with s0 = round(h0 / 2 + (h1 / 2 - s1) q10 / q00).
	 * Put (t0, t1) = (h0 - 2 * s0, h1 - 2 * s1) in FFT representation.
	 */
	if (logn <= 3) {
		h0 = h[0];
		for (u = 0; u < n; u ++) {
			s0[u] = fpr_rint(fpr_half(fpr_add(fpr_of(h0 & 1), t0[u])));
			h0 >>= 1;
		}
	} else {
		for (u = 0, w = 0; w < n; u ++) {
			h0 = h[u];
			for (v = 0; v < 8; v ++, w ++) {
				s0[w] = fpr_rint(fpr_half(fpr_add(fpr_of(h0 & 1), t0[w])));
				h0 >>= 1;
			}
		}
	}

	hash_to_fft(t0, h, s0, logn);
	return Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn)
		&& has_short_trace(t0, t1, q00, q10, q11, logn);
}

/* see inner.h */
int
Zf(verify_nearest_plane)(const uint8_t *restrict h,
	int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	unsigned logn, uint8_t *restrict tmp)
{
	/*
	 * This works slightly better than simple rounding, but is also slower.
	 * Reconstruct s0, by running Babai's NP algorithm with target
	 *     h0 / 2 + (h1 / 2 - s1) * q10 / q00.
	 */

	size_t n;
	fpr *t0, *t1;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;

	hash_to_fft(t0, h, NULL, logn);
	hash_to_fft(t1, SECOND_HASH(h, logn), s1, logn);

	Zf(poly_mul_fft)(t1, q10, logn);
	Zf(poly_div_autoadj_fft)(t1, q00, logn);
	Zf(poly_add)(t0, t1, logn);
	Zf(poly_mulconst)(t0, fpr_onehalf, logn);

	memcpy(t1, q00, n * sizeof *q00);
	/*
	 * Run Babai with target t0 and Gram-matrix q00.
	 */
	Zf(ffNearestPlane_dyn)(t0, t1, logn, t1 + n);
	Zf(fft_to_int16)(s0, t0, logn);

	/*
	 * Now run the casual verification.
	 */
	return Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn)
		&& Zf(uncompressed_verify)(h, s0, s1, q00, q10, q11, logn, tmp);
}

/* see inner.h */
int
Zf(verify)(const uint8_t *restrict h, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t u, v, w, n;
	uint8_t h0;
	int16_t s0w;
	fpr *t0, *t1;

	n = MKN(logn);
	t0 = (fpr *)tmp;
	t1 = t0 + n;

	hash_to_fft(t1, SECOND_HASH(h, logn), s1, logn);
	Zf(poly_prod_fft)(t0, t1, q10, logn);
	Zf(poly_div_autoadj_fft)(t0, q00, logn);
	Zf(iFFT)(t0, logn);

	/*
	 * Recover s0 with s0 = round(h0 / 2 + (h1 / 2 - s1) q10 / q00).
	 * Put (t0, t1) = (h0 - 2 * s0, h1 - 2 * s1) in FFT representation.
	 */
	if (logn <= 3) {
		h0 = h[0];
		for (u = 0; u < n; u ++) {
			s0w = fpr_rint(fpr_half(fpr_add(fpr_of(h0 & 1), t0[u])));
			t0[u] = fpr_of((h0 & 1) - 2 * s0w);
			h0 >>= 1;
		}
	} else {
		for (u = 0, w = 0; w < n; u ++) {
			h0 = h[u];
			for (v = 0; v < 8; v ++, w ++) {
				s0w = fpr_rint(fpr_half(fpr_add(fpr_of(h0 & 1), t0[w])));
				t0[w] = fpr_of((h0 & 1) - 2 * s0w);
				h0 >>= 1;
			}
		}
	}
	Zf(FFT)(t0, logn);

	return Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn)
		&& has_short_trace(t0, t1, q00, q10, q11, logn);
}

/* ============================================================================
 * Below is all the code necessary for Zf(uncompressed_verify_NTT).
 * Note that a lot of it is duplicate code with the code in keygen.c
 * ==========================================================================*/

/*
 * Precomputed small primes. Each element contains the following:
 *
 *  p   The prime itself.
 *
 *  g   A primitive root of phi = X^N+1 (in field Z_p).
 *
 *  s   The inverse of the product of all previous primes in the array,
 *      computed modulo p and in Montgomery representation.
 *
 * All primes are such that p = 1 mod 2048, and are lower than 2^31. They
 * are listed in decreasing order.
 */

typedef struct {
	uint32_t p, g, s;
} small_prime;

static const small_prime PRIMES[] = {
	{ 2147473409,  383167813,      10239 },
	{ 2147389441,  211808905,  471403745 },
	{ 2147387393,   37672282, 1329335065 },
/*	{ 2147377153, 1977035326,  968223422 },
	{ 2147358721, 1067163706,  132460015 },
	{ 2147352577, 1606082042,  598693809 },
	{ 2147346433, 2033915641, 1056257184 },
	{ 2147338241, 1653770625,  421286710 },
	{ 2147309569,  631200819, 1111201074 },
	{ 2147297281, 2038364663, 1042003613 },
	{ 0, 0, 0 } */
};

/*
 * Reduce a small signed integer modulo a small prime. The source
 * value x MUST be such that -p < x < p.
 */
static inline uint32_t
modp_set(int32_t x, uint32_t p)
{
	uint32_t w;

	w = (uint32_t)x;
	w += p & -(w >> 31);
	return w;
}

/*
 * Normalize a modular integer around 0.
 */
static inline int32_t
modp_norm(uint32_t x, uint32_t p)
{
	return (int32_t)(x - (p & (((x - ((p + 1) >> 1)) >> 31) - 1)));
}

/*
 * Compute -1/p mod 2^31. This works for all odd integers p that fit
 * on 31 bits.
 */
static uint32_t
modp_ninv31(uint32_t p)
{
	uint32_t y;

	y = 2 - p;
	y *= 2 - p * y;
	y *= 2 - p * y;
	y *= 2 - p * y;
	y *= 2 - p * y;
	return (uint32_t)0x7FFFFFFF & -y;
}

/*
 * Compute R = 2^31 mod p.
 */
static inline uint32_t
modp_R(uint32_t p)
{
	/*
	 * Since 2^30 < p < 2^31, we know that 2^31 mod p is simply
	 * 2^31 - p.
	 */
	return ((uint32_t)1 << 31) - p;
}

/*
 * Addition modulo p.
 */
static inline uint32_t
modp_add(uint32_t a, uint32_t b, uint32_t p)
{
	uint32_t d;

	d = a + b - p;
	d += p & -(d >> 31);
	return d;
}

/*
 * Subtraction modulo p.
 */
static inline uint32_t
modp_sub(uint32_t a, uint32_t b, uint32_t p)
{
	uint32_t d;

	d = a - b;
	d += p & -(d >> 31);
	return d;
}

/*
 * Montgomery multiplication modulo p. The 'p0i' value is -1/p mod 2^31.
 * It is required that p is an odd integer.
 */
static inline uint32_t
modp_montymul(uint32_t a, uint32_t b, uint32_t p, uint32_t p0i)
{
	uint64_t z, w;
	uint32_t d;

	z = (uint64_t)a * (uint64_t)b;
	w = ((z * p0i) & (uint64_t)0x7FFFFFFF) * p;
	d = (uint32_t)((z + w) >> 31) - p;
	d += p & -(d >> 31);
	return d;
}

/*
 * Compute R2 = 2^62 mod p.
 */
static uint32_t
modp_R2(uint32_t p, uint32_t p0i)
{
	uint32_t z;

	/*
	 * Compute z = 2^31 mod p (this is the value 1 in Montgomery
	 * representation), then double it with an addition.
	 */
	z = modp_R(p);
	z = modp_add(z, z, p);

	/*
	 * Square it five times to obtain 2^32 in Montgomery representation
	 * (i.e. 2^63 mod p).
	 */
	z = modp_montymul(z, z, p, p0i);
	z = modp_montymul(z, z, p, p0i);
	z = modp_montymul(z, z, p, p0i);
	z = modp_montymul(z, z, p, p0i);
	z = modp_montymul(z, z, p, p0i);

	/*
	 * Halve the value mod p to get 2^62.
	 */
	z = (z + (p & -(z & 1))) >> 1;
	return z;
}

/*
 * Division modulo p. If the divisor (b) is 0, then 0 is returned.
 * This function computes proper results only when p is prime.
 * Parameters:
 *   a     dividend
 *   b     divisor
 *   p     odd prime modulus
 *   p0i   -1/p mod 2^31
 *   R     2^31 mod R
 */
static uint32_t
modp_div(uint32_t a, uint32_t b, uint32_t p, uint32_t p0i, uint32_t R)
{
	uint32_t z, e;
	int i;

	e = p - 2;
	z = R;
	for (i = 30; i >= 0; i --) {
		uint32_t z2;

		z = modp_montymul(z, z, p, p0i);
		z2 = modp_montymul(z, b, p, p0i);
		z ^= (z ^ z2) & -(uint32_t)((e >> i) & 1);
	}

	/*
	 * The loop above just assumed that b was in Montgomery
	 * representation, i.e. really contained b*R; under that
	 * assumption, it returns 1/b in Montgomery representation,
	 * which is R/b. But we gave it b in normal representation,
	 * so the loop really returned R/(b/R) = R^2/b.
	 *
	 * We want a/b, so we need one Montgomery multiplication with a,
	 * which also remove one of the R factors, and another such
	 * multiplication to remove the second R factor.
	 */
	z = modp_montymul(z, 1, p, p0i);
	return modp_montymul(a, z, p, p0i);
}

/*
 * Bit-reversal index table.
 */
static const uint16_t REV10[] = {
	   0,  512,  256,  768,  128,  640,  384,  896,   64,  576,  320,  832,
	 192,  704,  448,  960,   32,  544,  288,  800,  160,  672,  416,  928,
	  96,  608,  352,  864,  224,  736,  480,  992,   16,  528,  272,  784,
	 144,  656,  400,  912,   80,  592,  336,  848,  208,  720,  464,  976,
	  48,  560,  304,  816,  176,  688,  432,  944,  112,  624,  368,  880,
	 240,  752,  496, 1008,    8,  520,  264,  776,  136,  648,  392,  904,
	  72,  584,  328,  840,  200,  712,  456,  968,   40,  552,  296,  808,
	 168,  680,  424,  936,  104,  616,  360,  872,  232,  744,  488, 1000,
	  24,  536,  280,  792,  152,  664,  408,  920,   88,  600,  344,  856,
	 216,  728,  472,  984,   56,  568,  312,  824,  184,  696,  440,  952,
	 120,  632,  376,  888,  248,  760,  504, 1016,    4,  516,  260,  772,
	 132,  644,  388,  900,   68,  580,  324,  836,  196,  708,  452,  964,
	  36,  548,  292,  804,  164,  676,  420,  932,  100,  612,  356,  868,
	 228,  740,  484,  996,   20,  532,  276,  788,  148,  660,  404,  916,
	  84,  596,  340,  852,  212,  724,  468,  980,   52,  564,  308,  820,
	 180,  692,  436,  948,  116,  628,  372,  884,  244,  756,  500, 1012,
	  12,  524,  268,  780,  140,  652,  396,  908,   76,  588,  332,  844,
	 204,  716,  460,  972,   44,  556,  300,  812,  172,  684,  428,  940,
	 108,  620,  364,  876,  236,  748,  492, 1004,   28,  540,  284,  796,
	 156,  668,  412,  924,   92,  604,  348,  860,  220,  732,  476,  988,
	  60,  572,  316,  828,  188,  700,  444,  956,  124,  636,  380,  892,
	 252,  764,  508, 1020,    2,  514,  258,  770,  130,  642,  386,  898,
	  66,  578,  322,  834,  194,  706,  450,  962,   34,  546,  290,  802,
	 162,  674,  418,  930,   98,  610,  354,  866,  226,  738,  482,  994,
	  18,  530,  274,  786,  146,  658,  402,  914,   82,  594,  338,  850,
	 210,  722,  466,  978,   50,  562,  306,  818,  178,  690,  434,  946,
	 114,  626,  370,  882,  242,  754,  498, 1010,   10,  522,  266,  778,
	 138,  650,  394,  906,   74,  586,  330,  842,  202,  714,  458,  970,
	  42,  554,  298,  810,  170,  682,  426,  938,  106,  618,  362,  874,
	 234,  746,  490, 1002,   26,  538,  282,  794,  154,  666,  410,  922,
	  90,  602,  346,  858,  218,  730,  474,  986,   58,  570,  314,  826,
	 186,  698,  442,  954,  122,  634,  378,  890,  250,  762,  506, 1018,
	   6,  518,  262,  774,  134,  646,  390,  902,   70,  582,  326,  838,
	 198,  710,  454,  966,   38,  550,  294,  806,  166,  678,  422,  934,
	 102,  614,  358,  870,  230,  742,  486,  998,   22,  534,  278,  790,
	 150,  662,  406,  918,   86,  598,  342,  854,  214,  726,  470,  982,
	  54,  566,  310,  822,  182,  694,  438,  950,  118,  630,  374,  886,
	 246,  758,  502, 1014,   14,  526,  270,  782,  142,  654,  398,  910,
	  78,  590,  334,  846,  206,  718,  462,  974,   46,  558,  302,  814,
	 174,  686,  430,  942,  110,  622,  366,  878,  238,  750,  494, 1006,
	  30,  542,  286,  798,  158,  670,  414,  926,   94,  606,  350,  862,
	 222,  734,  478,  990,   62,  574,  318,  830,  190,  702,  446,  958,
	 126,  638,  382,  894,  254,  766,  510, 1022,    1,  513,  257,  769,
	 129,  641,  385,  897,   65,  577,  321,  833,  193,  705,  449,  961,
	  33,  545,  289,  801,  161,  673,  417,  929,   97,  609,  353,  865,
	 225,  737,  481,  993,   17,  529,  273,  785,  145,  657,  401,  913,
	  81,  593,  337,  849,  209,  721,  465,  977,   49,  561,  305,  817,
	 177,  689,  433,  945,  113,  625,  369,  881,  241,  753,  497, 1009,
	   9,  521,  265,  777,  137,  649,  393,  905,   73,  585,  329,  841,
	 201,  713,  457,  969,   41,  553,  297,  809,  169,  681,  425,  937,
	 105,  617,  361,  873,  233,  745,  489, 1001,   25,  537,  281,  793,
	 153,  665,  409,  921,   89,  601,  345,  857,  217,  729,  473,  985,
	  57,  569,  313,  825,  185,  697,  441,  953,  121,  633,  377,  889,
	 249,  761,  505, 1017,    5,  517,  261,  773,  133,  645,  389,  901,
	  69,  581,  325,  837,  197,  709,  453,  965,   37,  549,  293,  805,
	 165,  677,  421,  933,  101,  613,  357,  869,  229,  741,  485,  997,
	  21,  533,  277,  789,  149,  661,  405,  917,   85,  597,  341,  853,
	 213,  725,  469,  981,   53,  565,  309,  821,  181,  693,  437,  949,
	 117,  629,  373,  885,  245,  757,  501, 1013,   13,  525,  269,  781,
	 141,  653,  397,  909,   77,  589,  333,  845,  205,  717,  461,  973,
	  45,  557,  301,  813,  173,  685,  429,  941,  109,  621,  365,  877,
	 237,  749,  493, 1005,   29,  541,  285,  797,  157,  669,  413,  925,
	  93,  605,  349,  861,  221,  733,  477,  989,   61,  573,  317,  829,
	 189,  701,  445,  957,  125,  637,  381,  893,  253,  765,  509, 1021,
	   3,  515,  259,  771,  131,  643,  387,  899,   67,  579,  323,  835,
	 195,  707,  451,  963,   35,  547,  291,  803,  163,  675,  419,  931,
	  99,  611,  355,  867,  227,  739,  483,  995,   19,  531,  275,  787,
	 147,  659,  403,  915,   83,  595,  339,  851,  211,  723,  467,  979,
	  51,  563,  307,  819,  179,  691,  435,  947,  115,  627,  371,  883,
	 243,  755,  499, 1011,   11,  523,  267,  779,  139,  651,  395,  907,
	  75,  587,  331,  843,  203,  715,  459,  971,   43,  555,  299,  811,
	 171,  683,  427,  939,  107,  619,  363,  875,  235,  747,  491, 1003,
	  27,  539,  283,  795,  155,  667,  411,  923,   91,  603,  347,  859,
	 219,  731,  475,  987,   59,  571,  315,  827,  187,  699,  443,  955,
	 123,  635,  379,  891,  251,  763,  507, 1019,    7,  519,  263,  775,
	 135,  647,  391,  903,   71,  583,  327,  839,  199,  711,  455,  967,
	  39,  551,  295,  807,  167,  679,  423,  935,  103,  615,  359,  871,
	 231,  743,  487,  999,   23,  535,  279,  791,  151,  663,  407,  919,
	  87,  599,  343,  855,  215,  727,  471,  983,   55,  567,  311,  823,
	 183,  695,  439,  951,  119,  631,  375,  887,  247,  759,  503, 1015,
	  15,  527,  271,  783,  143,  655,  399,  911,   79,  591,  335,  847,
	 207,  719,  463,  975,   47,  559,  303,  815,  175,  687,  431,  943,
	 111,  623,  367,  879,  239,  751,  495, 1007,   31,  543,  287,  799,
	 159,  671,  415,  927,   95,  607,  351,  863,  223,  735,  479,  991,
	  63,  575,  319,  831,  191,  703,  447,  959,  127,  639,  383,  895,
	 255,  767,  511, 1023
};

/*
 * Compute the roots for NTT and inverse NTT (binary case). Input
 * parameter g is a primitive 2048-th root of 1 modulo p (i.e. g^1024 =
 * -1 mod p). This fills gm[] and igm[] with powers of g and 1/g:
 *   gm[rev(i)] = g^i mod p
 *   igm[rev(i)] = (1/g)^i mod p
 * where rev() is the "bit reversal" function over 10 bits. It fills
 * the arrays only up to N = 2^logn values.
 *
 * The values stored in gm[] and igm[] are in Montgomery representation.
 *
 * p must be a prime such that p = 1 mod 2048.
 */
static void
modp_mkgm2(uint32_t *restrict gm, uint32_t *restrict igm, unsigned logn,
	uint32_t g, uint32_t p, uint32_t p0i)
{
	size_t u, n;
	unsigned k;
	uint32_t ig, x1, x2, R2;

	n = MKN(logn);

	/*
	 * We want g such that g^(2N) = 1 mod p, but the provided
	 * generator has order 2048. We must square it a few times.
	 */
	R2 = modp_R2(p, p0i);
	g = modp_montymul(g, R2, p, p0i);
	for (k = logn; k < 10; k ++) {
		g = modp_montymul(g, g, p, p0i);
	}

	ig = modp_div(R2, g, p, p0i, modp_R(p));
	k = 10 - logn;
	x1 = x2 = modp_R(p);
	for (u = 0; u < n; u ++) {
		size_t v;

		v = REV10[u << k];
		gm[v] = x1;
		igm[v] = x2;
		x1 = modp_montymul(x1, g, p, p0i);
		x2 = modp_montymul(x2, ig, p, p0i);
	}
}

/*
 * Compute the NTT over a polynomial a (binary case).
 */
static void
modp_NTT2(uint32_t *a, const uint32_t *gm, unsigned logn,
	uint32_t p, uint32_t p0i)
{
	size_t t, m, n;

	if (logn == 0) {
		return;
	}
	n = MKN(logn);
	t = n;
	for (m = 1; m < n; m <<= 1) {
		size_t ht, u, v1;

		ht = t >> 1;
		for (u = 0, v1 = 0; u < m; u ++, v1 += t) {
			uint32_t s;
			size_t v;
			uint32_t *r1, *r2;

			s = gm[m + u];
			r1 = a + v1;
			r2 = r1 + ht;
			for (v = 0; v < ht; v ++, r1 ++, r2 ++) {
				uint32_t x, y;

				x = *r1;
				y = modp_montymul(*r2, s, p, p0i);
				*r1 = modp_add(x, y, p);
				*r2 = modp_sub(x, y, p);
			}
		}
		t = ht;
	}
}

/*
 * Compute the inverse NTT over a polynomial (binary case).
 */
static void
modp_iNTT2(uint32_t *a, const uint32_t *igm, unsigned logn,
	uint32_t p, uint32_t p0i)
{
	size_t t, m, n, k;
	uint32_t ni;
	uint32_t *r;

	if (logn == 0) {
		return;
	}
	n = MKN(logn);
	t = 1;
	for (m = n; m > 1; m >>= 1) {
		size_t hm, dt, u, v1;

		hm = m >> 1;
		dt = t << 1;
		for (u = 0, v1 = 0; u < hm; u ++, v1 += dt) {
			uint32_t s;
			size_t v;
			uint32_t *r1, *r2;

			s = igm[hm + u];
			r1 = a + v1;
			r2 = r1 + t;
			for (v = 0; v < t; v ++, r1 ++, r2 ++) {
				uint32_t x, y;

				x = *r1;
				y = *r2;
				*r1 = modp_add(x, y, p);
				*r2 = modp_montymul(
					modp_sub(x, y, p), s, p, p0i);
			}
		}
		t = dt;
	}

	/*
	 * We need 1/n in Montgomery representation, i.e. R/n. Since
	 * 1 <= logn <= 10, R/n is an integer; morever, R/n <= 2^30 < p,
	 * thus a simple shift will do.
	 */
	ni = (uint32_t)1 << (31 - logn);
	for (k = 0, r = a; k < n; k ++, r ++) {
		*r = modp_montymul(*r, ni, p, p0i);
	}
}

static void
int16_to_ntt(uint32_t *a, const int16_t *x, const uint32_t *gm,
	uint32_t p, uint32_t p0i, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u++) {
		a[u] = modp_set(x[u], p);
	}

	modp_NTT2(a, gm, logn, p, p0i);
}

/*
 * Set the polynomial a equal to h - 2 * s and converts it to Montgomery NTT
 * representation.
 */
static void
hash_to_ntt(uint32_t *a, const uint8_t *h, const int16_t *s,
	uint32_t *gm, uint32_t R2, uint32_t p, uint32_t p0i, unsigned logn)
{
	size_t n, u, v;
	int32_t hash;

	n = MKN(logn);
	if (logn <= 3) {
		hash = h[0];
		for (u = 0; u < n; u ++) {
			a[u] = modp_set(((hash >> u) & 1) - 2 * s[u], p);
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				a[u] = modp_set((hash & 1) - 2 * s[u], p);
				hash >>= 1;
			}
		}
	}

	modp_NTT2(a, gm, logn, p, p0i);
	for (u = 0; u < n; u++) {
		a[u] = modp_montymul(a[u], R2, p, p0i);
	}
}

/*
 * This a helper function for Zf(uncompressed_verify_NTT) and Zf(verify_NTT),
 * since these both calculate the norm of (s0, s1) modulo various primes with
 * respect to a quadratic form Q.
 * However, in Zf(verify_NTT), some of the computations are already done.
 */
static int
verify_norm_NTT(const uint8_t *restrict h,
	const int16_t *restrict s0, const int16_t *restrict s1,
	const int16_t *restrict q00, const int16_t *restrict q10,
	uint32_t p, uint32_t p0i, uint32_t R, uint32_t R2,
	unsigned logn, uint8_t *restrict tmp)
{
	/*
	 * TODO: since q00, q11 are selfadjoint, write a NTT, iNTT version for
	 * selfadjoints that only uses first half of coefficients, saving 4n bytes.
	 */

	uint32_t norm, norm0, term;
	uint32_t *s0_i32, *s1_i32, *r0_i32, *r1_i32, *gm, *igm;
	int32_t *q11;
	size_t n, hn, u, v;

	n = MKN(logn);
	hn = n >> 1;
	norm = 0;

	s0_i32 = (uint32_t *)tmp;
	s1_i32 = s0_i32 + n;
	r0_i32 = s1_i32 + n;
	r1_i32 = r0_i32 + n;
	q11 = (int32_t *)r1_i32;
	gm = r1_i32 + n;
	igm = gm + n;

	hash_to_ntt(s0_i32, h, s0, gm, R2, p, p0i, logn);
	hash_to_ntt(s1_i32, SECOND_HASH(h, logn), s1, gm, R2, p, p0i, logn);

	/*
	 * Store q00 in r0_i32 and assume q10 is in r1_i32 (NTT representation).
	 */
	int16_to_ntt(r0_i32, q00, gm, p, p0i, logn);

	// s0* q00 s0
	for (u = 0; u < hn; u++) {
		term = modp_montymul(s0_i32[u], s0_i32[n - 1 - u], p, p0i);
		norm = modp_add(norm, modp_montymul(term, r0_i32[u], p, p0i), p);
	}

	// s0* q01 s1 + s1* q10 s0
	for (u = 0; u < n; u++) {
		term = modp_montymul(s1_i32[u], s0_i32[n - 1 - u], p, p0i);
		norm = modp_add(norm, modp_montymul(term, r1_i32[u], p, p0i), p);
	}

	/*
	 * Now determine q11 = (1 + q10*q01) / q00.
	 */
	for (u = 0; u < hn; u++) {
		r1_i32[u] = modp_montymul(r1_i32[u], R2, p, p0i);
		r1_i32[u] = modp_montymul(r1_i32[u], r1_i32[n - 1 - u], p, p0i);
		r1_i32[u] = modp_add(r1_i32[u], 1u, p);
		r1_i32[u] = modp_div(r1_i32[u], r0_i32[u], p, p0i, R);

		// Needed for the iNTT:
		r1_i32[n - 1 - u] = r1_i32[u];
	}

	// s1* q11 s1
	for (u = 0; u < hn; u++) {
		term = modp_montymul(s1_i32[u], s1_i32[n - 1 - u], p, p0i);
		norm = modp_add(norm, modp_montymul(term, r1_i32[u], p, p0i), p);
	}

	/*
	 * Store q11 for later checks, since taking the inverse (modp_div) is
	 * expensive modulo a prime p.
	 */
	modp_iNTT2(r1_i32, igm, logn, p, p0i);
	for (u = 0; u < n; u++) {
		q11[u] = modp_norm(r1_i32[u], p);
	}

	norm0 = norm;
	for (v = 1; v < 2; v++) {
		norm = 0;
		p = PRIMES[v].p;
		p0i = modp_ninv31(p);
		R = modp_R(p);
		R2 = modp_R2(p, p0i);
		modp_mkgm2(gm, igm, logn, PRIMES[v].g, p, p0i);

		hash_to_ntt(s0_i32, h, s0, gm, R2, p, p0i, logn);
		hash_to_ntt(s1_i32, SECOND_HASH(h, logn), s1, gm, R2, p, p0i, logn);

		// s0* q00 s0
		int16_to_ntt(r0_i32, q00, gm, p, p0i, logn);
		for (u = 0; u < hn; u++) {
			term = modp_montymul(s0_i32[u], s0_i32[n - 1 - u], p, p0i);
			norm = modp_add(norm, modp_montymul(term, r0_i32[u], p, p0i), p);
		}

		// s0* q01 s1 + s1* q10 s0
		int16_to_ntt(r0_i32, q10, gm, p, p0i, logn);
		for (u = 0; u < n; u++) {
			term = modp_montymul(s1_i32[u], s0_i32[n - 1 - u], p, p0i);
			norm = modp_add(norm, modp_montymul(term, r0_i32[u], p, p0i), p);
		}

		// s1* q11 s1
		// int32_to_ntt(r0_i32, q11, gm, p, p0i, logn);
		for (u = 0; u < n; u++)
			r0_i32[u] = modp_set(q11[u], p);
		modp_NTT2(r0_i32, gm, logn, p, p0i);

		for (u = 0; u < hn; u++) {
			term = modp_montymul(s1_i32[u], s1_i32[n - 1 - u], p, p0i);
			norm = modp_add(norm, modp_montymul(term, r0_i32[u], p, p0i), p);
		}

		// norm0 := \infty, when norm0 != norm.
		norm0 |= -((norm0 - norm) >> 31);
		norm0 |= -((norm - norm0) >> 31);
	}

	return Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn)
		&& (norm0 & (hn - 1)) == 0
		&& (norm0 >> (logn - 1)) <= L2BOUND(logn);
}

/* see inner.h */
int
Zf(uncompressed_verify_NTT)(const uint8_t *restrict h,
	const int16_t *restrict s0, const int16_t *restrict s1,
	const int16_t *restrict q00, const int16_t *restrict q10,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n;
	uint32_t *gm, *igm, *q10_i32, p, p0i, R, R2;

	n = MKN(logn);
	p = PRIMES[0].p;
	p0i = modp_ninv31(p);
	R = modp_R(p);
	R2 = modp_R2(p, p0i);

	q10_i32 = ((uint32_t *)tmp) + 3*n;
	gm = q10_i32 + n;
	igm = gm + n;

	modp_mkgm2(gm, igm, logn, PRIMES[0].g, p, p0i);
	int16_to_ntt(q10_i32, q10, gm, p, p0i, logn);

	return verify_norm_NTT(h, s0, s1, q00, q10, p, p0i, R, R2, logn, tmp);
}

// ================================================================================
// NEW CODE
// ================================================================================
#define FT double

const FT vrfy_fpr_gm_tab[] = {
	0, 0, /* unused */
	-0.000000000000000000000000000, 1.000000000000000000000000000,
	 0.707106781186547524400844362, 0.707106781186547524400844362,
	-0.707106781186547524400844362, 0.707106781186547524400844362,
	 0.923879532511286756128183189, 0.382683432365089771728459984,
	-0.382683432365089771728459984, 0.923879532511286756128183189,
	 0.382683432365089771728459984, 0.923879532511286756128183189,
	-0.923879532511286756128183189, 0.382683432365089771728459984,
	 0.980785280403230449126182236, 0.195090322016128267848284868,
	-0.195090322016128267848284868, 0.980785280403230449126182236,
	 0.555570233019602224742830814, 0.831469612302545237078788378,
	-0.831469612302545237078788378, 0.555570233019602224742830814,
	 0.831469612302545237078788378, 0.555570233019602224742830814,
	-0.555570233019602224742830814, 0.831469612302545237078788378,
	 0.195090322016128267848284868, 0.980785280403230449126182236,
	-0.980785280403230449126182236, 0.195090322016128267848284868,
	 0.995184726672196886244836953, 0.098017140329560601994195564,
	-0.098017140329560601994195564, 0.995184726672196886244836953,
	 0.634393284163645498215171613, 0.773010453362736960810906610,
	-0.773010453362736960810906610, 0.634393284163645498215171613,
	 0.881921264348355029712756864, 0.471396736825997648556387626,
	-0.471396736825997648556387626, 0.881921264348355029712756864,
	 0.290284677254462367636192376, 0.956940335732208864935797887,
	-0.956940335732208864935797887, 0.290284677254462367636192376,
	 0.956940335732208864935797887, 0.290284677254462367636192376,
	-0.290284677254462367636192376, 0.956940335732208864935797887,
	 0.471396736825997648556387626, 0.881921264348355029712756864,
	-0.881921264348355029712756864, 0.471396736825997648556387626,
	 0.773010453362736960810906610, 0.634393284163645498215171613,
	-0.634393284163645498215171613, 0.773010453362736960810906610,
	 0.098017140329560601994195564, 0.995184726672196886244836953,
	-0.995184726672196886244836953, 0.098017140329560601994195564,
	 0.998795456205172392714771605, 0.049067674327418014254954977,
	-0.049067674327418014254954977, 0.998795456205172392714771605,
	 0.671558954847018400625376850, 0.740951125354959091175616897,
	-0.740951125354959091175616897, 0.671558954847018400625376850,
	 0.903989293123443331586200297, 0.427555093430282094320966857,
	-0.427555093430282094320966857, 0.903989293123443331586200297,
	 0.336889853392220050689253213, 0.941544065183020778412509403,
	-0.941544065183020778412509403, 0.336889853392220050689253213,
	 0.970031253194543992603984207, 0.242980179903263889948274162,
	-0.242980179903263889948274162, 0.970031253194543992603984207,
	 0.514102744193221726593693839, 0.857728610000272069902269984,
	-0.857728610000272069902269984, 0.514102744193221726593693839,
	 0.803207531480644909806676513, 0.595699304492433343467036529,
	-0.595699304492433343467036529, 0.803207531480644909806676513,
	 0.146730474455361751658850130, 0.989176509964780973451673738,
	-0.989176509964780973451673738, 0.146730474455361751658850130,
	 0.989176509964780973451673738, 0.146730474455361751658850130,
	-0.146730474455361751658850130, 0.989176509964780973451673738,
	 0.595699304492433343467036529, 0.803207531480644909806676513,
	-0.803207531480644909806676513, 0.595699304492433343467036529,
	 0.857728610000272069902269984, 0.514102744193221726593693839,
	-0.514102744193221726593693839, 0.857728610000272069902269984,
	 0.242980179903263889948274162, 0.970031253194543992603984207,
	-0.970031253194543992603984207, 0.242980179903263889948274162,
	 0.941544065183020778412509403, 0.336889853392220050689253213,
	-0.336889853392220050689253213, 0.941544065183020778412509403,
	 0.427555093430282094320966857, 0.903989293123443331586200297,
	-0.903989293123443331586200297, 0.427555093430282094320966857,
	 0.740951125354959091175616897, 0.671558954847018400625376850,
	-0.671558954847018400625376850, 0.740951125354959091175616897,
	 0.049067674327418014254954977, 0.998795456205172392714771605,
	-0.998795456205172392714771605, 0.049067674327418014254954977,
	 0.999698818696204220115765650, 0.024541228522912288031734529,
	-0.024541228522912288031734529, 0.999698818696204220115765650,
	 0.689540544737066924616730630, 0.724247082951466920941069243,
	-0.724247082951466920941069243, 0.689540544737066924616730630,
	 0.914209755703530654635014829, 0.405241314004989870908481306,
	-0.405241314004989870908481306, 0.914209755703530654635014829,
	 0.359895036534988148775104572, 0.932992798834738887711660256,
	-0.932992798834738887711660256, 0.359895036534988148775104572,
	 0.975702130038528544460395766, 0.219101240156869797227737547,
	-0.219101240156869797227737547, 0.975702130038528544460395766,
	 0.534997619887097210663076905, 0.844853565249707073259571205,
	-0.844853565249707073259571205, 0.534997619887097210663076905,
	 0.817584813151583696504920884, 0.575808191417845300745972454,
	-0.575808191417845300745972454, 0.817584813151583696504920884,
	 0.170961888760301226363642357, 0.985277642388941244774018433,
	-0.985277642388941244774018433, 0.170961888760301226363642357,
	 0.992479534598709998156767252, 0.122410675199216198498704474,
	-0.122410675199216198498704474, 0.992479534598709998156767252,
	 0.615231590580626845484913563, 0.788346427626606262009164705,
	-0.788346427626606262009164705, 0.615231590580626845484913563,
	 0.870086991108711418652292404, 0.492898192229784036873026689,
	-0.492898192229784036873026689, 0.870086991108711418652292404,
	 0.266712757474898386325286515, 0.963776065795439866686464356,
	-0.963776065795439866686464356, 0.266712757474898386325286515,
	 0.949528180593036667195936074, 0.313681740398891476656478846,
	-0.313681740398891476656478846, 0.949528180593036667195936074,
	 0.449611329654606600046294579, 0.893224301195515320342416447,
	-0.893224301195515320342416447, 0.449611329654606600046294579,
	 0.757208846506484547575464054, 0.653172842953776764084203014,
	-0.653172842953776764084203014, 0.757208846506484547575464054,
	 0.073564563599667423529465622, 0.997290456678690216135597140,
	-0.997290456678690216135597140, 0.073564563599667423529465622,
	 0.997290456678690216135597140, 0.073564563599667423529465622,
	-0.073564563599667423529465622, 0.997290456678690216135597140,
	 0.653172842953776764084203014, 0.757208846506484547575464054,
	-0.757208846506484547575464054, 0.653172842953776764084203014,
	 0.893224301195515320342416447, 0.449611329654606600046294579,
	-0.449611329654606600046294579, 0.893224301195515320342416447,
	 0.313681740398891476656478846, 0.949528180593036667195936074,
	-0.949528180593036667195936074, 0.313681740398891476656478846,
	 0.963776065795439866686464356, 0.266712757474898386325286515,
	-0.266712757474898386325286515, 0.963776065795439866686464356,
	 0.492898192229784036873026689, 0.870086991108711418652292404,
	-0.870086991108711418652292404, 0.492898192229784036873026689,
	 0.788346427626606262009164705, 0.615231590580626845484913563,
	-0.615231590580626845484913563, 0.788346427626606262009164705,
	 0.122410675199216198498704474, 0.992479534598709998156767252,
	-0.992479534598709998156767252, 0.122410675199216198498704474,
	 0.985277642388941244774018433, 0.170961888760301226363642357,
	-0.170961888760301226363642357, 0.985277642388941244774018433,
	 0.575808191417845300745972454, 0.817584813151583696504920884,
	-0.817584813151583696504920884, 0.575808191417845300745972454,
	 0.844853565249707073259571205, 0.534997619887097210663076905,
	-0.534997619887097210663076905, 0.844853565249707073259571205,
	 0.219101240156869797227737547, 0.975702130038528544460395766,
	-0.975702130038528544460395766, 0.219101240156869797227737547,
	 0.932992798834738887711660256, 0.359895036534988148775104572,
	-0.359895036534988148775104572, 0.932992798834738887711660256,
	 0.405241314004989870908481306, 0.914209755703530654635014829,
	-0.914209755703530654635014829, 0.405241314004989870908481306,
	 0.724247082951466920941069243, 0.689540544737066924616730630,
	-0.689540544737066924616730630, 0.724247082951466920941069243,
	 0.024541228522912288031734529, 0.999698818696204220115765650,
	-0.999698818696204220115765650, 0.024541228522912288031734529,
	 0.999924701839144540921646491, 0.012271538285719926079408262,
	-0.012271538285719926079408262, 0.999924701839144540921646491,
	 0.698376249408972853554813503, 0.715730825283818654125532623,
	-0.715730825283818654125532623, 0.698376249408972853554813503,
	 0.919113851690057743908477789, 0.393992040061048108596188661,
	-0.393992040061048108596188661, 0.919113851690057743908477789,
	 0.371317193951837543411934967, 0.928506080473215565937167396,
	-0.928506080473215565937167396, 0.371317193951837543411934967,
	 0.978317370719627633106240097, 0.207111376192218549708116020,
	-0.207111376192218549708116020, 0.978317370719627633106240097,
	 0.545324988422046422313987347, 0.838224705554838043186996856,
	-0.838224705554838043186996856, 0.545324988422046422313987347,
	 0.824589302785025264474803737, 0.565731810783613197389765011,
	-0.565731810783613197389765011, 0.824589302785025264474803737,
	 0.183039887955140958516532578, 0.983105487431216327180301155,
	-0.983105487431216327180301155, 0.183039887955140958516532578,
	 0.993906970002356041546922813, 0.110222207293883058807899140,
	-0.110222207293883058807899140, 0.993906970002356041546922813,
	 0.624859488142386377084072816, 0.780737228572094478301588484,
	-0.780737228572094478301588484, 0.624859488142386377084072816,
	 0.876070094195406607095844268, 0.482183772079122748517344481,
	-0.482183772079122748517344481, 0.876070094195406607095844268,
	 0.278519689385053105207848526, 0.960430519415565811199035138,
	-0.960430519415565811199035138, 0.278519689385053105207848526,
	 0.953306040354193836916740383, 0.302005949319228067003463232,
	-0.302005949319228067003463232, 0.953306040354193836916740383,
	 0.460538710958240023633181487, 0.887639620402853947760181617,
	-0.887639620402853947760181617, 0.460538710958240023633181487,
	 0.765167265622458925888815999, 0.643831542889791465068086063,
	-0.643831542889791465068086063, 0.765167265622458925888815999,
	 0.085797312344439890461556332, 0.996312612182778012627226190,
	-0.996312612182778012627226190, 0.085797312344439890461556332,
	 0.998118112900149207125155861, 0.061320736302208577782614593,
	-0.061320736302208577782614593, 0.998118112900149207125155861,
	 0.662415777590171761113069817, 0.749136394523459325469203257,
	-0.749136394523459325469203257, 0.662415777590171761113069817,
	 0.898674465693953843041976744, 0.438616238538527637647025738,
	-0.438616238538527637647025738, 0.898674465693953843041976744,
	 0.325310292162262934135954708, 0.945607325380521325730945387,
	-0.945607325380521325730945387, 0.325310292162262934135954708,
	 0.966976471044852109087220226, 0.254865659604514571553980779,
	-0.254865659604514571553980779, 0.966976471044852109087220226,
	 0.503538383725717558691867071, 0.863972856121586737918147054,
	-0.863972856121586737918147054, 0.503538383725717558691867071,
	 0.795836904608883536262791915, 0.605511041404325513920626941,
	-0.605511041404325513920626941, 0.795836904608883536262791915,
	 0.134580708507126186316358409, 0.990902635427780025108237011,
	-0.990902635427780025108237011, 0.134580708507126186316358409,
	 0.987301418157858382399815802, 0.158858143333861441684385360,
	-0.158858143333861441684385360, 0.987301418157858382399815802,
	 0.585797857456438860328080838, 0.810457198252594791726703434,
	-0.810457198252594791726703434, 0.585797857456438860328080838,
	 0.851355193105265142261290312, 0.524589682678468906215098464,
	-0.524589682678468906215098464, 0.851355193105265142261290312,
	 0.231058108280671119643236018, 0.972939952205560145467720114,
	-0.972939952205560145467720114, 0.231058108280671119643236018,
	 0.937339011912574923201899593, 0.348418680249434568419308588,
	-0.348418680249434568419308588, 0.937339011912574923201899593,
	 0.416429560097637182562598911, 0.909167983090522376563884788,
	-0.909167983090522376563884788, 0.416429560097637182562598911,
	 0.732654271672412834615546649, 0.680600997795453050594430464,
	-0.680600997795453050594430464, 0.732654271672412834615546649,
	 0.036807222941358832324332691, 0.999322384588349500896221011,
	-0.999322384588349500896221011, 0.036807222941358832324332691,
	 0.999322384588349500896221011, 0.036807222941358832324332691,
	-0.036807222941358832324332691, 0.999322384588349500896221011,
	 0.680600997795453050594430464, 0.732654271672412834615546649,
	-0.732654271672412834615546649, 0.680600997795453050594430464,
	 0.909167983090522376563884788, 0.416429560097637182562598911,
	-0.416429560097637182562598911, 0.909167983090522376563884788,
	 0.348418680249434568419308588, 0.937339011912574923201899593,
	-0.937339011912574923201899593, 0.348418680249434568419308588,
	 0.972939952205560145467720114, 0.231058108280671119643236018,
	-0.231058108280671119643236018, 0.972939952205560145467720114,
	 0.524589682678468906215098464, 0.851355193105265142261290312,
	-0.851355193105265142261290312, 0.524589682678468906215098464,
	 0.810457198252594791726703434, 0.585797857456438860328080838,
	-0.585797857456438860328080838, 0.810457198252594791726703434,
	 0.158858143333861441684385360, 0.987301418157858382399815802,
	-0.987301418157858382399815802, 0.158858143333861441684385360,
	 0.990902635427780025108237011, 0.134580708507126186316358409,
	-0.134580708507126186316358409, 0.990902635427780025108237011,
	 0.605511041404325513920626941, 0.795836904608883536262791915,
	-0.795836904608883536262791915, 0.605511041404325513920626941,
	 0.863972856121586737918147054, 0.503538383725717558691867071,
	-0.503538383725717558691867071, 0.863972856121586737918147054,
	 0.254865659604514571553980779, 0.966976471044852109087220226,
	-0.966976471044852109087220226, 0.254865659604514571553980779,
	 0.945607325380521325730945387, 0.325310292162262934135954708,
	-0.325310292162262934135954708, 0.945607325380521325730945387,
	 0.438616238538527637647025738, 0.898674465693953843041976744,
	-0.898674465693953843041976744, 0.438616238538527637647025738,
	 0.749136394523459325469203257, 0.662415777590171761113069817,
	-0.662415777590171761113069817, 0.749136394523459325469203257,
	 0.061320736302208577782614593, 0.998118112900149207125155861,
	-0.998118112900149207125155861, 0.061320736302208577782614593,
	 0.996312612182778012627226190, 0.085797312344439890461556332,
	-0.085797312344439890461556332, 0.996312612182778012627226190,
	 0.643831542889791465068086063, 0.765167265622458925888815999,
	-0.765167265622458925888815999, 0.643831542889791465068086063,
	 0.887639620402853947760181617, 0.460538710958240023633181487,
	-0.460538710958240023633181487, 0.887639620402853947760181617,
	 0.302005949319228067003463232, 0.953306040354193836916740383,
	-0.953306040354193836916740383, 0.302005949319228067003463232,
	 0.960430519415565811199035138, 0.278519689385053105207848526,
	-0.278519689385053105207848526, 0.960430519415565811199035138,
	 0.482183772079122748517344481, 0.876070094195406607095844268,
	-0.876070094195406607095844268, 0.482183772079122748517344481,
	 0.780737228572094478301588484, 0.624859488142386377084072816,
	-0.624859488142386377084072816, 0.780737228572094478301588484,
	 0.110222207293883058807899140, 0.993906970002356041546922813,
	-0.993906970002356041546922813, 0.110222207293883058807899140,
	 0.983105487431216327180301155, 0.183039887955140958516532578,
	-0.183039887955140958516532578, 0.983105487431216327180301155,
	 0.565731810783613197389765011, 0.824589302785025264474803737,
	-0.824589302785025264474803737, 0.565731810783613197389765011,
	 0.838224705554838043186996856, 0.545324988422046422313987347,
	-0.545324988422046422313987347, 0.838224705554838043186996856,
	 0.207111376192218549708116020, 0.978317370719627633106240097,
	-0.978317370719627633106240097, 0.207111376192218549708116020,
	 0.928506080473215565937167396, 0.371317193951837543411934967,
	-0.371317193951837543411934967, 0.928506080473215565937167396,
	 0.393992040061048108596188661, 0.919113851690057743908477789,
	-0.919113851690057743908477789, 0.393992040061048108596188661,
	 0.715730825283818654125532623, 0.698376249408972853554813503,
	-0.698376249408972853554813503, 0.715730825283818654125532623,
	 0.012271538285719926079408262, 0.999924701839144540921646491,
	-0.999924701839144540921646491, 0.012271538285719926079408262,
	 0.999981175282601142656990438, 0.006135884649154475359640235,
	-0.006135884649154475359640235, 0.999981175282601142656990438,
	 0.702754744457225302452914421, 0.711432195745216441522130290,
	-0.711432195745216441522130290, 0.702754744457225302452914421,
	 0.921514039342041943465396332, 0.388345046698826291624993541,
	-0.388345046698826291624993541, 0.921514039342041943465396332,
	 0.377007410216418256726567823, 0.926210242138311341974793388,
	-0.926210242138311341974793388, 0.377007410216418256726567823,
	 0.979569765685440534439326110, 0.201104634842091911558443546,
	-0.201104634842091911558443546, 0.979569765685440534439326110,
	 0.550457972936604802977289893, 0.834862874986380056304401383,
	-0.834862874986380056304401383, 0.550457972936604802977289893,
	 0.828045045257755752067527592, 0.560661576197336023839710223,
	-0.560661576197336023839710223, 0.828045045257755752067527592,
	 0.189068664149806212754997837, 0.981963869109555264072848154,
	-0.981963869109555264072848154, 0.189068664149806212754997837,
	 0.994564570734255452119106243, 0.104121633872054579120943880,
	-0.104121633872054579120943880, 0.994564570734255452119106243,
	 0.629638238914927025372981341, 0.776888465673232450040827983,
	-0.776888465673232450040827983, 0.629638238914927025372981341,
	 0.879012226428633477831323711, 0.476799230063322133342158117,
	-0.476799230063322133342158117, 0.879012226428633477831323711,
	 0.284407537211271843618310615, 0.958703474895871555374645792,
	-0.958703474895871555374645792, 0.284407537211271843618310615,
	 0.955141168305770721498157712, 0.296150888243623824121786128,
	-0.296150888243623824121786128, 0.955141168305770721498157712,
	 0.465976495767966177902756065, 0.884797098430937780104007041,
	-0.884797098430937780104007041, 0.465976495767966177902756065,
	 0.769103337645579639346626069, 0.639124444863775743801488193,
	-0.639124444863775743801488193, 0.769103337645579639346626069,
	 0.091908956497132728624990979, 0.995767414467659793982495643,
	-0.995767414467659793982495643, 0.091908956497132728624990979,
	 0.998475580573294752208559038, 0.055195244349689939809447526,
	-0.055195244349689939809447526, 0.998475580573294752208559038,
	 0.666999922303637506650154222, 0.745057785441465962407907310,
	-0.745057785441465962407907310, 0.666999922303637506650154222,
	 0.901348847046022014570746093, 0.433093818853151968484222638,
	-0.433093818853151968484222638, 0.901348847046022014570746093,
	 0.331106305759876401737190737, 0.943593458161960361495301445,
	-0.943593458161960361495301445, 0.331106305759876401737190737,
	 0.968522094274417316221088329, 0.248927605745720168110682816,
	-0.248927605745720168110682816, 0.968522094274417316221088329,
	 0.508830142543107036931749324, 0.860866938637767279344583877,
	-0.860866938637767279344583877, 0.508830142543107036931749324,
	 0.799537269107905033500246232, 0.600616479383868926653875896,
	-0.600616479383868926653875896, 0.799537269107905033500246232,
	 0.140658239332849230714788846, 0.990058210262297105505906464,
	-0.990058210262297105505906464, 0.140658239332849230714788846,
	 0.988257567730749491404792538, 0.152797185258443427720336613,
	-0.152797185258443427720336613, 0.988257567730749491404792538,
	 0.590759701858874228423887908, 0.806847553543799272206514313,
	-0.806847553543799272206514313, 0.590759701858874228423887908,
	 0.854557988365400520767862276, 0.519355990165589587361829932,
	-0.519355990165589587361829932, 0.854557988365400520767862276,
	 0.237023605994367206867735915, 0.971503890986251775537099622,
	-0.971503890986251775537099622, 0.237023605994367206867735915,
	 0.939459223602189911962669246, 0.342660717311994397592781983,
	-0.342660717311994397592781983, 0.939459223602189911962669246,
	 0.422000270799799685941287941, 0.906595704514915365332960588,
	-0.906595704514915365332960588, 0.422000270799799685941287941,
	 0.736816568877369875090132520, 0.676092703575315960360419228,
	-0.676092703575315960360419228, 0.736816568877369875090132520,
	 0.042938256934940823077124540, 0.999077727752645382888781997,
	-0.999077727752645382888781997, 0.042938256934940823077124540,
	 0.999529417501093163079703322, 0.030674803176636625934021028,
	-0.030674803176636625934021028, 0.999529417501093163079703322,
	 0.685083667772700381362052545, 0.728464390448225196492035438,
	-0.728464390448225196492035438, 0.685083667772700381362052545,
	 0.911706032005429851404397325, 0.410843171057903942183466675,
	-0.410843171057903942183466675, 0.911706032005429851404397325,
	 0.354163525420490382357395796, 0.935183509938947577642207480,
	-0.935183509938947577642207480, 0.354163525420490382357395796,
	 0.974339382785575860518721668, 0.225083911359792835991642120,
	-0.225083911359792835991642120, 0.974339382785575860518721668,
	 0.529803624686294668216054671, 0.848120344803297251279133563,
	-0.848120344803297251279133563, 0.529803624686294668216054671,
	 0.814036329705948361654516690, 0.580813958095764545075595272,
	-0.580813958095764545075595272, 0.814036329705948361654516690,
	 0.164913120489969921418189113, 0.986308097244598647863297524,
	-0.986308097244598647863297524, 0.164913120489969921418189113,
	 0.991709753669099522860049931, 0.128498110793793172624415589,
	-0.128498110793793172624415589, 0.991709753669099522860049931,
	 0.610382806276309452716352152, 0.792106577300212351782342879,
	-0.792106577300212351782342879, 0.610382806276309452716352152,
	 0.867046245515692651480195629, 0.498227666972781852410983869,
	-0.498227666972781852410983869, 0.867046245515692651480195629,
	 0.260794117915275518280186509, 0.965394441697689374550843858,
	-0.965394441697689374550843858, 0.260794117915275518280186509,
	 0.947585591017741134653387321, 0.319502030816015677901518272,
	-0.319502030816015677901518272, 0.947585591017741134653387321,
	 0.444122144570429231642069418, 0.895966249756185155914560282,
	-0.895966249756185155914560282, 0.444122144570429231642069418,
	 0.753186799043612482483430486, 0.657806693297078656931182264,
	-0.657806693297078656931182264, 0.753186799043612482483430486,
	 0.067443919563664057897972422, 0.997723066644191609848546728,
	-0.997723066644191609848546728, 0.067443919563664057897972422,
	 0.996820299291165714972629398, 0.079682437971430121147120656,
	-0.079682437971430121147120656, 0.996820299291165714972629398,
	 0.648514401022112445084560551, 0.761202385484261814029709836,
	-0.761202385484261814029709836, 0.648514401022112445084560551,
	 0.890448723244757889952150560, 0.455083587126343823535869268,
	-0.455083587126343823535869268, 0.890448723244757889952150560,
	 0.307849640041534893682063646, 0.951435020969008369549175569,
	-0.951435020969008369549175569, 0.307849640041534893682063646,
	 0.962121404269041595429604316, 0.272621355449948984493347477,
	-0.272621355449948984493347477, 0.962121404269041595429604316,
	 0.487550160148435954641485027, 0.873094978418290098636085973,
	-0.873094978418290098636085973, 0.487550160148435954641485027,
	 0.784556597155575233023892575, 0.620057211763289178646268191,
	-0.620057211763289178646268191, 0.784556597155575233023892575,
	 0.116318630911904767252544319, 0.993211949234794533104601012,
	-0.993211949234794533104601012, 0.116318630911904767252544319,
	 0.984210092386929073193874387, 0.177004220412148756196839844,
	-0.177004220412148756196839844, 0.984210092386929073193874387,
	 0.570780745886967280232652864, 0.821102514991104679060430820,
	-0.821102514991104679060430820, 0.570780745886967280232652864,
	 0.841554977436898409603499520, 0.540171472729892881297845480,
	-0.540171472729892881297845480, 0.841554977436898409603499520,
	 0.213110319916091373967757518, 0.977028142657754351485866211,
	-0.977028142657754351485866211, 0.213110319916091373967757518,
	 0.930766961078983731944872340, 0.365612997804773870011745909,
	-0.365612997804773870011745909, 0.930766961078983731944872340,
	 0.399624199845646828544117031, 0.916679059921042663116457013,
	-0.916679059921042663116457013, 0.399624199845646828544117031,
	 0.720002507961381629076682999, 0.693971460889654009003734389,
	-0.693971460889654009003734389, 0.720002507961381629076682999,
	 0.018406729905804820927366313, 0.999830581795823422015722275,
	-0.999830581795823422015722275, 0.018406729905804820927366313,
	 0.999830581795823422015722275, 0.018406729905804820927366313,
	-0.018406729905804820927366313, 0.999830581795823422015722275,
	 0.693971460889654009003734389, 0.720002507961381629076682999,
	-0.720002507961381629076682999, 0.693971460889654009003734389,
	 0.916679059921042663116457013, 0.399624199845646828544117031,
	-0.399624199845646828544117031, 0.916679059921042663116457013,
	 0.365612997804773870011745909, 0.930766961078983731944872340,
	-0.930766961078983731944872340, 0.365612997804773870011745909,
	 0.977028142657754351485866211, 0.213110319916091373967757518,
	-0.213110319916091373967757518, 0.977028142657754351485866211,
	 0.540171472729892881297845480, 0.841554977436898409603499520,
	-0.841554977436898409603499520, 0.540171472729892881297845480,
	 0.821102514991104679060430820, 0.570780745886967280232652864,
	-0.570780745886967280232652864, 0.821102514991104679060430820,
	 0.177004220412148756196839844, 0.984210092386929073193874387,
	-0.984210092386929073193874387, 0.177004220412148756196839844,
	 0.993211949234794533104601012, 0.116318630911904767252544319,
	-0.116318630911904767252544319, 0.993211949234794533104601012,
	 0.620057211763289178646268191, 0.784556597155575233023892575,
	-0.784556597155575233023892575, 0.620057211763289178646268191,
	 0.873094978418290098636085973, 0.487550160148435954641485027,
	-0.487550160148435954641485027, 0.873094978418290098636085973,
	 0.272621355449948984493347477, 0.962121404269041595429604316,
	-0.962121404269041595429604316, 0.272621355449948984493347477,
	 0.951435020969008369549175569, 0.307849640041534893682063646,
	-0.307849640041534893682063646, 0.951435020969008369549175569,
	 0.455083587126343823535869268, 0.890448723244757889952150560,
	-0.890448723244757889952150560, 0.455083587126343823535869268,
	 0.761202385484261814029709836, 0.648514401022112445084560551,
	-0.648514401022112445084560551, 0.761202385484261814029709836,
	 0.079682437971430121147120656, 0.996820299291165714972629398,
	-0.996820299291165714972629398, 0.079682437971430121147120656,
	 0.997723066644191609848546728, 0.067443919563664057897972422,
	-0.067443919563664057897972422, 0.997723066644191609848546728,
	 0.657806693297078656931182264, 0.753186799043612482483430486,
	-0.753186799043612482483430486, 0.657806693297078656931182264,
	 0.895966249756185155914560282, 0.444122144570429231642069418,
	-0.444122144570429231642069418, 0.895966249756185155914560282,
	 0.319502030816015677901518272, 0.947585591017741134653387321,
	-0.947585591017741134653387321, 0.319502030816015677901518272,
	 0.965394441697689374550843858, 0.260794117915275518280186509,
	-0.260794117915275518280186509, 0.965394441697689374550843858,
	 0.498227666972781852410983869, 0.867046245515692651480195629,
	-0.867046245515692651480195629, 0.498227666972781852410983869,
	 0.792106577300212351782342879, 0.610382806276309452716352152,
	-0.610382806276309452716352152, 0.792106577300212351782342879,
	 0.128498110793793172624415589, 0.991709753669099522860049931,
	-0.991709753669099522860049931, 0.128498110793793172624415589,
	 0.986308097244598647863297524, 0.164913120489969921418189113,
	-0.164913120489969921418189113, 0.986308097244598647863297524,
	 0.580813958095764545075595272, 0.814036329705948361654516690,
	-0.814036329705948361654516690, 0.580813958095764545075595272,
	 0.848120344803297251279133563, 0.529803624686294668216054671,
	-0.529803624686294668216054671, 0.848120344803297251279133563,
	 0.225083911359792835991642120, 0.974339382785575860518721668,
	-0.974339382785575860518721668, 0.225083911359792835991642120,
	 0.935183509938947577642207480, 0.354163525420490382357395796,
	-0.354163525420490382357395796, 0.935183509938947577642207480,
	 0.410843171057903942183466675, 0.911706032005429851404397325,
	-0.911706032005429851404397325, 0.410843171057903942183466675,
	 0.728464390448225196492035438, 0.685083667772700381362052545,
	-0.685083667772700381362052545, 0.728464390448225196492035438,
	 0.030674803176636625934021028, 0.999529417501093163079703322,
	-0.999529417501093163079703322, 0.030674803176636625934021028,
	 0.999077727752645382888781997, 0.042938256934940823077124540,
	-0.042938256934940823077124540, 0.999077727752645382888781997,
	 0.676092703575315960360419228, 0.736816568877369875090132520,
	-0.736816568877369875090132520, 0.676092703575315960360419228,
	 0.906595704514915365332960588, 0.422000270799799685941287941,
	-0.422000270799799685941287941, 0.906595704514915365332960588,
	 0.342660717311994397592781983, 0.939459223602189911962669246,
	-0.939459223602189911962669246, 0.342660717311994397592781983,
	 0.971503890986251775537099622, 0.237023605994367206867735915,
	-0.237023605994367206867735915, 0.971503890986251775537099622,
	 0.519355990165589587361829932, 0.854557988365400520767862276,
	-0.854557988365400520767862276, 0.519355990165589587361829932,
	 0.806847553543799272206514313, 0.590759701858874228423887908,
	-0.590759701858874228423887908, 0.806847553543799272206514313,
	 0.152797185258443427720336613, 0.988257567730749491404792538,
	-0.988257567730749491404792538, 0.152797185258443427720336613,
	 0.990058210262297105505906464, 0.140658239332849230714788846,
	-0.140658239332849230714788846, 0.990058210262297105505906464,
	 0.600616479383868926653875896, 0.799537269107905033500246232,
	-0.799537269107905033500246232, 0.600616479383868926653875896,
	 0.860866938637767279344583877, 0.508830142543107036931749324,
	-0.508830142543107036931749324, 0.860866938637767279344583877,
	 0.248927605745720168110682816, 0.968522094274417316221088329,
	-0.968522094274417316221088329, 0.248927605745720168110682816,
	 0.943593458161960361495301445, 0.331106305759876401737190737,
	-0.331106305759876401737190737, 0.943593458161960361495301445,
	 0.433093818853151968484222638, 0.901348847046022014570746093,
	-0.901348847046022014570746093, 0.433093818853151968484222638,
	 0.745057785441465962407907310, 0.666999922303637506650154222,
	-0.666999922303637506650154222, 0.745057785441465962407907310,
	 0.055195244349689939809447526, 0.998475580573294752208559038,
	-0.998475580573294752208559038, 0.055195244349689939809447526,
	 0.995767414467659793982495643, 0.091908956497132728624990979,
	-0.091908956497132728624990979, 0.995767414467659793982495643,
	 0.639124444863775743801488193, 0.769103337645579639346626069,
	-0.769103337645579639346626069, 0.639124444863775743801488193,
	 0.884797098430937780104007041, 0.465976495767966177902756065,
	-0.465976495767966177902756065, 0.884797098430937780104007041,
	 0.296150888243623824121786128, 0.955141168305770721498157712,
	-0.955141168305770721498157712, 0.296150888243623824121786128,
	 0.958703474895871555374645792, 0.284407537211271843618310615,
	-0.284407537211271843618310615, 0.958703474895871555374645792,
	 0.476799230063322133342158117, 0.879012226428633477831323711,
	-0.879012226428633477831323711, 0.476799230063322133342158117,
	 0.776888465673232450040827983, 0.629638238914927025372981341,
	-0.629638238914927025372981341, 0.776888465673232450040827983,
	 0.104121633872054579120943880, 0.994564570734255452119106243,
	-0.994564570734255452119106243, 0.104121633872054579120943880,
	 0.981963869109555264072848154, 0.189068664149806212754997837,
	-0.189068664149806212754997837, 0.981963869109555264072848154,
	 0.560661576197336023839710223, 0.828045045257755752067527592,
	-0.828045045257755752067527592, 0.560661576197336023839710223,
	 0.834862874986380056304401383, 0.550457972936604802977289893,
	-0.550457972936604802977289893, 0.834862874986380056304401383,
	 0.201104634842091911558443546, 0.979569765685440534439326110,
	-0.979569765685440534439326110, 0.201104634842091911558443546,
	 0.926210242138311341974793388, 0.377007410216418256726567823,
	-0.377007410216418256726567823, 0.926210242138311341974793388,
	 0.388345046698826291624993541, 0.921514039342041943465396332,
	-0.921514039342041943465396332, 0.388345046698826291624993541,
	 0.711432195745216441522130290, 0.702754744457225302452914421,
	-0.702754744457225302452914421, 0.711432195745216441522130290,
	 0.006135884649154475359640235, 0.999981175282601142656990438,
	-0.999981175282601142656990438, 0.006135884649154475359640235,
	 0.999995293809576171511580126, 0.003067956762965976270145365,
	-0.003067956762965976270145365, 0.999995293809576171511580126,
	 0.704934080375904908852523758, 0.709272826438865651316533772,
	-0.709272826438865651316533772, 0.704934080375904908852523758,
	 0.922701128333878570437264227, 0.385516053843918864075607949,
	-0.385516053843918864075607949, 0.922701128333878570437264227,
	 0.379847208924051170576281147, 0.925049240782677590302371869,
	-0.925049240782677590302371869, 0.379847208924051170576281147,
	 0.980182135968117392690210009, 0.198098410717953586179324918,
	-0.198098410717953586179324918, 0.980182135968117392690210009,
	 0.553016705580027531764226988, 0.833170164701913186439915922,
	-0.833170164701913186439915922, 0.553016705580027531764226988,
	 0.829761233794523042469023765, 0.558118531220556115693702964,
	-0.558118531220556115693702964, 0.829761233794523042469023765,
	 0.192080397049892441679288205, 0.981379193313754574318224190,
	-0.981379193313754574318224190, 0.192080397049892441679288205,
	 0.994879330794805620591166107, 0.101069862754827824987887585,
	-0.101069862754827824987887585, 0.994879330794805620591166107,
	 0.632018735939809021909403706, 0.774953106594873878359129282,
	-0.774953106594873878359129282, 0.632018735939809021909403706,
	 0.880470889052160770806542929, 0.474100214650550014398580015,
	-0.474100214650550014398580015, 0.880470889052160770806542929,
	 0.287347459544729526477331841, 0.957826413027532890321037029,
	-0.957826413027532890321037029, 0.287347459544729526477331841,
	 0.956045251349996443270479823, 0.293219162694258650606608599,
	-0.293219162694258650606608599, 0.956045251349996443270479823,
	 0.468688822035827933697617870, 0.883363338665731594736308015,
	-0.883363338665731594736308015, 0.468688822035827933697617870,
	 0.771060524261813773200605759, 0.636761861236284230413943435,
	-0.636761861236284230413943435, 0.771060524261813773200605759,
	 0.094963495329638998938034312, 0.995480755491926941769171600,
	-0.995480755491926941769171600, 0.094963495329638998938034312,
	 0.998640218180265222418199049, 0.052131704680283321236358216,
	-0.052131704680283321236358216, 0.998640218180265222418199049,
	 0.669282588346636065720696366, 0.743007952135121693517362293,
	-0.743007952135121693517362293, 0.669282588346636065720696366,
	 0.902673318237258806751502391, 0.430326481340082633908199031,
	-0.430326481340082633908199031, 0.902673318237258806751502391,
	 0.333999651442009404650865481, 0.942573197601446879280758735,
	-0.942573197601446879280758735, 0.333999651442009404650865481,
	 0.969281235356548486048290738, 0.245955050335794611599924709,
	-0.245955050335794611599924709, 0.969281235356548486048290738,
	 0.511468850437970399504391001, 0.859301818357008404783582139,
	-0.859301818357008404783582139, 0.511468850437970399504391001,
	 0.801376171723140219430247777, 0.598160706996342311724958652,
	-0.598160706996342311724958652, 0.801376171723140219430247777,
	 0.143695033150294454819773349, 0.989622017463200834623694454,
	-0.989622017463200834623694454, 0.143695033150294454819773349,
	 0.988721691960323767604516485, 0.149764534677321517229695737,
	-0.149764534677321517229695737, 0.988721691960323767604516485,
	 0.593232295039799808047809426, 0.805031331142963597922659282,
	-0.805031331142963597922659282, 0.593232295039799808047809426,
	 0.856147328375194481019630732, 0.516731799017649881508753876,
	-0.516731799017649881508753876, 0.856147328375194481019630732,
	 0.240003022448741486568922365, 0.970772140728950302138169611,
	-0.970772140728950302138169611, 0.240003022448741486568922365,
	 0.940506070593268323787291309, 0.339776884406826857828825803,
	-0.339776884406826857828825803, 0.940506070593268323787291309,
	 0.424779681209108833357226189, 0.905296759318118774354048329,
	-0.905296759318118774354048329, 0.424779681209108833357226189,
	 0.738887324460615147933116508, 0.673829000378756060917568372,
	-0.673829000378756060917568372, 0.738887324460615147933116508,
	 0.046003182130914628814301788, 0.998941293186856850633930266,
	-0.998941293186856850633930266, 0.046003182130914628814301788,
	 0.999618822495178597116830637, 0.027608145778965741612354872,
	-0.027608145778965741612354872, 0.999618822495178597116830637,
	 0.687315340891759108199186948, 0.726359155084345976817494315,
	-0.726359155084345976817494315, 0.687315340891759108199186948,
	 0.912962190428398164628018233, 0.408044162864978680820747499,
	-0.408044162864978680820747499, 0.912962190428398164628018233,
	 0.357030961233430032614954036, 0.934092550404258914729877883,
	-0.934092550404258914729877883, 0.357030961233430032614954036,
	 0.975025345066994146844913468, 0.222093620973203534094094721,
	-0.222093620973203534094094721, 0.975025345066994146844913468,
	 0.532403127877197971442805218, 0.846490938774052078300544488,
	-0.846490938774052078300544488, 0.532403127877197971442805218,
	 0.815814410806733789010772660, 0.578313796411655563342245019,
	-0.578313796411655563342245019, 0.815814410806733789010772660,
	 0.167938294974731178054745536, 0.985797509167567424700995000,
	-0.985797509167567424700995000, 0.167938294974731178054745536,
	 0.992099313142191757112085445, 0.125454983411546238542336453,
	-0.125454983411546238542336453, 0.992099313142191757112085445,
	 0.612810082429409703935211936, 0.790230221437310055030217152,
	-0.790230221437310055030217152, 0.612810082429409703935211936,
	 0.868570705971340895340449876, 0.495565261825772531150266670,
	-0.495565261825772531150266670, 0.868570705971340895340449876,
	 0.263754678974831383611349322, 0.964589793289812723836432159,
	-0.964589793289812723836432159, 0.263754678974831383611349322,
	 0.948561349915730288158494826, 0.316593375556165867243047035,
	-0.316593375556165867243047035, 0.948561349915730288158494826,
	 0.446868840162374195353044389, 0.894599485631382678433072126,
	-0.894599485631382678433072126, 0.446868840162374195353044389,
	 0.755201376896536527598710756, 0.655492852999615385312679701,
	-0.655492852999615385312679701, 0.755201376896536527598710756,
	 0.070504573389613863027351471, 0.997511456140303459699448390,
	-0.997511456140303459699448390, 0.070504573389613863027351471,
	 0.997060070339482978987989949, 0.076623861392031492278332463,
	-0.076623861392031492278332463, 0.997060070339482978987989949,
	 0.650846684996380915068975573, 0.759209188978388033485525443,
	-0.759209188978388033485525443, 0.650846684996380915068975573,
	 0.891840709392342727796478697, 0.452349587233770874133026703,
	-0.452349587233770874133026703, 0.891840709392342727796478697,
	 0.310767152749611495835997250, 0.950486073949481721759926101,
	-0.950486073949481721759926101, 0.310767152749611495835997250,
	 0.962953266873683886347921481, 0.269668325572915106525464462,
	-0.269668325572915106525464462, 0.962953266873683886347921481,
	 0.490226483288291154229598449, 0.871595086655951034842481435,
	-0.871595086655951034842481435, 0.490226483288291154229598449,
	 0.786455213599085757522319464, 0.617647307937803932403979402,
	-0.617647307937803932403979402, 0.786455213599085757522319464,
	 0.119365214810991364593637790, 0.992850414459865090793563344,
	-0.992850414459865090793563344, 0.119365214810991364593637790,
	 0.984748501801904218556553176, 0.173983873387463827950700807,
	-0.173983873387463827950700807, 0.984748501801904218556553176,
	 0.573297166698042212820171239, 0.819347520076796960824689637,
	-0.819347520076796960824689637, 0.573297166698042212820171239,
	 0.843208239641845437161743865, 0.537587076295645482502214932,
	-0.537587076295645482502214932, 0.843208239641845437161743865,
	 0.216106797076219509948385131, 0.976369731330021149312732194,
	-0.976369731330021149312732194, 0.216106797076219509948385131,
	 0.931884265581668106718557199, 0.362755724367397216204854462,
	-0.362755724367397216204854462, 0.931884265581668106718557199,
	 0.402434650859418441082533934, 0.915448716088267819566431292,
	-0.915448716088267819566431292, 0.402434650859418441082533934,
	 0.722128193929215321243607198, 0.691759258364157774906734132,
	-0.691759258364157774906734132, 0.722128193929215321243607198,
	 0.021474080275469507418374898, 0.999769405351215321657617036,
	-0.999769405351215321657617036, 0.021474080275469507418374898,
	 0.999882347454212525633049627, 0.015339206284988101044151868,
	-0.015339206284988101044151868, 0.999882347454212525633049627,
	 0.696177131491462944788582591, 0.717870045055731736211325329,
	-0.717870045055731736211325329, 0.696177131491462944788582591,
	 0.917900775621390457642276297, 0.396809987416710328595290911,
	-0.396809987416710328595290911, 0.917900775621390457642276297,
	 0.368466829953372331712746222, 0.929640895843181265457918066,
	-0.929640895843181265457918066, 0.368466829953372331712746222,
	 0.977677357824509979943404762, 0.210111836880469621717489972,
	-0.210111836880469621717489972, 0.977677357824509979943404762,
	 0.542750784864515906586768661, 0.839893794195999504583383987,
	-0.839893794195999504583383987, 0.542750784864515906586768661,
	 0.822849781375826332046780034, 0.568258952670131549790548489,
	-0.568258952670131549790548489, 0.822849781375826332046780034,
	 0.180022901405699522679906590, 0.983662419211730274396237776,
	-0.983662419211730274396237776, 0.180022901405699522679906590,
	 0.993564135520595333782021697, 0.113270952177564349018228733,
	-0.113270952177564349018228733, 0.993564135520595333782021697,
	 0.622461279374149972519166721, 0.782650596166575738458949301,
	-0.782650596166575738458949301, 0.622461279374149972519166721,
	 0.874586652278176112634431897, 0.484869248000791101822951699,
	-0.484869248000791101822951699, 0.874586652278176112634431897,
	 0.275571819310958163076425168, 0.961280485811320641748659653,
	-0.961280485811320641748659653, 0.275571819310958163076425168,
	 0.952375012719765858529893608, 0.304929229735402406490728633,
	-0.304929229735402406490728633, 0.952375012719765858529893608,
	 0.457813303598877221904961155, 0.889048355854664562540777729,
	-0.889048355854664562540777729, 0.457813303598877221904961155,
	 0.763188417263381271704838297, 0.646176012983316364832802220,
	-0.646176012983316364832802220, 0.763188417263381271704838297,
	 0.082740264549375693111987083, 0.996571145790554847093566910,
	-0.996571145790554847093566910, 0.082740264549375693111987083,
	 0.997925286198596012623025462, 0.064382630929857460819324537,
	-0.064382630929857460819324537, 0.997925286198596012623025462,
	 0.660114342067420478559490747, 0.751165131909686411205819422,
	-0.751165131909686411205819422, 0.660114342067420478559490747,
	 0.897324580705418281231391836, 0.441371268731716692879988968,
	-0.441371268731716692879988968, 0.897324580705418281231391836,
	 0.322407678801069848384807478, 0.946600913083283570044599823,
	-0.946600913083283570044599823, 0.322407678801069848384807478,
	 0.966190003445412555433832961, 0.257831102162159005614471295,
	-0.257831102162159005614471295, 0.966190003445412555433832961,
	 0.500885382611240786241285004, 0.865513624090569082825488358,
	-0.865513624090569082825488358, 0.500885382611240786241285004,
	 0.793975477554337164895083757, 0.607949784967773667243642671,
	-0.607949784967773667243642671, 0.793975477554337164895083757,
	 0.131540028702883111103387493, 0.991310859846115418957349799,
	-0.991310859846115418957349799, 0.131540028702883111103387493,
	 0.986809401814185476970235952, 0.161886393780111837641387995,
	-0.161886393780111837641387995, 0.986809401814185476970235952,
	 0.583308652937698294392830961, 0.812250586585203913049744181,
	-0.812250586585203913049744181, 0.583308652937698294392830961,
	 0.849741768000852489471268395, 0.527199134781901348464274575,
	-0.527199134781901348464274575, 0.849741768000852489471268395,
	 0.228072083170885739254457379, 0.973644249650811925318383912,
	-0.973644249650811925318383912, 0.228072083170885739254457379,
	 0.936265667170278246576310996, 0.351292756085567125601307623,
	-0.351292756085567125601307623, 0.936265667170278246576310996,
	 0.413638312238434547471944324, 0.910441292258067196934095369,
	-0.910441292258067196934095369, 0.413638312238434547471944324,
	 0.730562769227827561177758850, 0.682845546385248068164596123,
	-0.682845546385248068164596123, 0.730562769227827561177758850,
	 0.033741171851377584833716112, 0.999430604555461772019008327,
	-0.999430604555461772019008327, 0.033741171851377584833716112,
	 0.999204758618363895492950001, 0.039872927587739811128578738,
	-0.039872927587739811128578738, 0.999204758618363895492950001,
	 0.678350043129861486873655042, 0.734738878095963464563223604,
	-0.734738878095963464563223604, 0.678350043129861486873655042,
	 0.907886116487666212038681480, 0.419216888363223956433010020,
	-0.419216888363223956433010020, 0.907886116487666212038681480,
	 0.345541324963989065539191723, 0.938403534063108112192420774,
	-0.938403534063108112192420774, 0.345541324963989065539191723,
	 0.972226497078936305708321144, 0.234041958583543423191242045,
	-0.234041958583543423191242045, 0.972226497078936305708321144,
	 0.521975292937154342694258318, 0.852960604930363657746588082,
	-0.852960604930363657746588082, 0.521975292937154342694258318,
	 0.808656181588174991946968128, 0.588281548222645304786439813,
	-0.588281548222645304786439813, 0.808656181588174991946968128,
	 0.155828397654265235743101486, 0.987784141644572154230969032,
	-0.987784141644572154230969032, 0.155828397654265235743101486,
	 0.990485084256457037998682243, 0.137620121586486044948441663,
	-0.137620121586486044948441663, 0.990485084256457037998682243,
	 0.603066598540348201693430617, 0.797690840943391108362662755,
	-0.797690840943391108362662755, 0.603066598540348201693430617,
	 0.862423956111040538690933878, 0.506186645345155291048942344,
	-0.506186645345155291048942344, 0.862423956111040538690933878,
	 0.251897818154216950498106628, 0.967753837093475465243391912,
	-0.967753837093475465243391912, 0.251897818154216950498106628,
	 0.944604837261480265659265493, 0.328209843579092526107916817,
	-0.328209843579092526107916817, 0.944604837261480265659265493,
	 0.435857079922255491032544080, 0.900015892016160228714535267,
	-0.900015892016160228714535267, 0.435857079922255491032544080,
	 0.747100605980180144323078847, 0.664710978203344868130324985,
	-0.664710978203344868130324985, 0.747100605980180144323078847,
	 0.058258264500435759613979782, 0.998301544933892840738782163,
	-0.998301544933892840738782163, 0.058258264500435759613979782,
	 0.996044700901251989887944810, 0.088853552582524596561586535,
	-0.088853552582524596561586535, 0.996044700901251989887944810,
	 0.641481012808583151988739898, 0.767138911935820381181694573,
	-0.767138911935820381181694573, 0.641481012808583151988739898,
	 0.886222530148880631647990821, 0.463259783551860197390719637,
	-0.463259783551860197390719637, 0.886222530148880631647990821,
	 0.299079826308040476750336973, 0.954228095109105629780430732,
	-0.954228095109105629780430732, 0.299079826308040476750336973,
	 0.959571513081984528335528181, 0.281464937925757984095231007,
	-0.281464937925757984095231007, 0.959571513081984528335528181,
	 0.479493757660153026679839798, 0.877545290207261291668470750,
	-0.877545290207261291668470750, 0.479493757660153026679839798,
	 0.778816512381475953374724325, 0.627251815495144113509622565,
	-0.627251815495144113509622565, 0.778816512381475953374724325,
	 0.107172424956808849175529148, 0.994240449453187946358413442,
	-0.994240449453187946358413442, 0.107172424956808849175529148,
	 0.982539302287441255907040396, 0.186055151663446648105438304,
	-0.186055151663446648105438304, 0.982539302287441255907040396,
	 0.563199344013834115007363772, 0.826321062845663480311195452,
	-0.826321062845663480311195452, 0.563199344013834115007363772,
	 0.836547727223511984524285790, 0.547894059173100165608820571,
	-0.547894059173100165608820571, 0.836547727223511984524285790,
	 0.204108966092816874181696950, 0.978948175319062194715480124,
	-0.978948175319062194715480124, 0.204108966092816874181696950,
	 0.927362525650401087274536959, 0.374164062971457997104393020,
	-0.374164062971457997104393020, 0.927362525650401087274536959,
	 0.391170384302253888687512949, 0.920318276709110566440076541,
	-0.920318276709110566440076541, 0.391170384302253888687512949,
	 0.713584868780793592903125099, 0.700568793943248366792866380,
	-0.700568793943248366792866380, 0.713584868780793592903125099,
	 0.009203754782059819315102378, 0.999957644551963866333120920,
	-0.999957644551963866333120920, 0.009203754782059819315102378,
	 0.999957644551963866333120920, 0.009203754782059819315102378,
	-0.009203754782059819315102378, 0.999957644551963866333120920,
	 0.700568793943248366792866380, 0.713584868780793592903125099,
	-0.713584868780793592903125099, 0.700568793943248366792866380,
	 0.920318276709110566440076541, 0.391170384302253888687512949,
	-0.391170384302253888687512949, 0.920318276709110566440076541,
	 0.374164062971457997104393020, 0.927362525650401087274536959,
	-0.927362525650401087274536959, 0.374164062971457997104393020,
	 0.978948175319062194715480124, 0.204108966092816874181696950,
	-0.204108966092816874181696950, 0.978948175319062194715480124,
	 0.547894059173100165608820571, 0.836547727223511984524285790,
	-0.836547727223511984524285790, 0.547894059173100165608820571,
	 0.826321062845663480311195452, 0.563199344013834115007363772,
	-0.563199344013834115007363772, 0.826321062845663480311195452,
	 0.186055151663446648105438304, 0.982539302287441255907040396,
	-0.982539302287441255907040396, 0.186055151663446648105438304,
	 0.994240449453187946358413442, 0.107172424956808849175529148,
	-0.107172424956808849175529148, 0.994240449453187946358413442,
	 0.627251815495144113509622565, 0.778816512381475953374724325,
	-0.778816512381475953374724325, 0.627251815495144113509622565,
	 0.877545290207261291668470750, 0.479493757660153026679839798,
	-0.479493757660153026679839798, 0.877545290207261291668470750,
	 0.281464937925757984095231007, 0.959571513081984528335528181,
	-0.959571513081984528335528181, 0.281464937925757984095231007,
	 0.954228095109105629780430732, 0.299079826308040476750336973,
	-0.299079826308040476750336973, 0.954228095109105629780430732,
	 0.463259783551860197390719637, 0.886222530148880631647990821,
	-0.886222530148880631647990821, 0.463259783551860197390719637,
	 0.767138911935820381181694573, 0.641481012808583151988739898,
	-0.641481012808583151988739898, 0.767138911935820381181694573,
	 0.088853552582524596561586535, 0.996044700901251989887944810,
	-0.996044700901251989887944810, 0.088853552582524596561586535,
	 0.998301544933892840738782163, 0.058258264500435759613979782,
	-0.058258264500435759613979782, 0.998301544933892840738782163,
	 0.664710978203344868130324985, 0.747100605980180144323078847,
	-0.747100605980180144323078847, 0.664710978203344868130324985,
	 0.900015892016160228714535267, 0.435857079922255491032544080,
	-0.435857079922255491032544080, 0.900015892016160228714535267,
	 0.328209843579092526107916817, 0.944604837261480265659265493,
	-0.944604837261480265659265493, 0.328209843579092526107916817,
	 0.967753837093475465243391912, 0.251897818154216950498106628,
	-0.251897818154216950498106628, 0.967753837093475465243391912,
	 0.506186645345155291048942344, 0.862423956111040538690933878,
	-0.862423956111040538690933878, 0.506186645345155291048942344,
	 0.797690840943391108362662755, 0.603066598540348201693430617,
	-0.603066598540348201693430617, 0.797690840943391108362662755,
	 0.137620121586486044948441663, 0.990485084256457037998682243,
	-0.990485084256457037998682243, 0.137620121586486044948441663,
	 0.987784141644572154230969032, 0.155828397654265235743101486,
	-0.155828397654265235743101486, 0.987784141644572154230969032,
	 0.588281548222645304786439813, 0.808656181588174991946968128,
	-0.808656181588174991946968128, 0.588281548222645304786439813,
	 0.852960604930363657746588082, 0.521975292937154342694258318,
	-0.521975292937154342694258318, 0.852960604930363657746588082,
	 0.234041958583543423191242045, 0.972226497078936305708321144,
	-0.972226497078936305708321144, 0.234041958583543423191242045,
	 0.938403534063108112192420774, 0.345541324963989065539191723,
	-0.345541324963989065539191723, 0.938403534063108112192420774,
	 0.419216888363223956433010020, 0.907886116487666212038681480,
	-0.907886116487666212038681480, 0.419216888363223956433010020,
	 0.734738878095963464563223604, 0.678350043129861486873655042,
	-0.678350043129861486873655042, 0.734738878095963464563223604,
	 0.039872927587739811128578738, 0.999204758618363895492950001,
	-0.999204758618363895492950001, 0.039872927587739811128578738,
	 0.999430604555461772019008327, 0.033741171851377584833716112,
	-0.033741171851377584833716112, 0.999430604555461772019008327,
	 0.682845546385248068164596123, 0.730562769227827561177758850,
	-0.730562769227827561177758850, 0.682845546385248068164596123,
	 0.910441292258067196934095369, 0.413638312238434547471944324,
	-0.413638312238434547471944324, 0.910441292258067196934095369,
	 0.351292756085567125601307623, 0.936265667170278246576310996,
	-0.936265667170278246576310996, 0.351292756085567125601307623,
	 0.973644249650811925318383912, 0.228072083170885739254457379,
	-0.228072083170885739254457379, 0.973644249650811925318383912,
	 0.527199134781901348464274575, 0.849741768000852489471268395,
	-0.849741768000852489471268395, 0.527199134781901348464274575,
	 0.812250586585203913049744181, 0.583308652937698294392830961,
	-0.583308652937698294392830961, 0.812250586585203913049744181,
	 0.161886393780111837641387995, 0.986809401814185476970235952,
	-0.986809401814185476970235952, 0.161886393780111837641387995,
	 0.991310859846115418957349799, 0.131540028702883111103387493,
	-0.131540028702883111103387493, 0.991310859846115418957349799,
	 0.607949784967773667243642671, 0.793975477554337164895083757,
	-0.793975477554337164895083757, 0.607949784967773667243642671,
	 0.865513624090569082825488358, 0.500885382611240786241285004,
	-0.500885382611240786241285004, 0.865513624090569082825488358,
	 0.257831102162159005614471295, 0.966190003445412555433832961,
	-0.966190003445412555433832961, 0.257831102162159005614471295,
	 0.946600913083283570044599823, 0.322407678801069848384807478,
	-0.322407678801069848384807478, 0.946600913083283570044599823,
	 0.441371268731716692879988968, 0.897324580705418281231391836,
	-0.897324580705418281231391836, 0.441371268731716692879988968,
	 0.751165131909686411205819422, 0.660114342067420478559490747,
	-0.660114342067420478559490747, 0.751165131909686411205819422,
	 0.064382630929857460819324537, 0.997925286198596012623025462,
	-0.997925286198596012623025462, 0.064382630929857460819324537,
	 0.996571145790554847093566910, 0.082740264549375693111987083,
	-0.082740264549375693111987083, 0.996571145790554847093566910,
	 0.646176012983316364832802220, 0.763188417263381271704838297,
	-0.763188417263381271704838297, 0.646176012983316364832802220,
	 0.889048355854664562540777729, 0.457813303598877221904961155,
	-0.457813303598877221904961155, 0.889048355854664562540777729,
	 0.304929229735402406490728633, 0.952375012719765858529893608,
	-0.952375012719765858529893608, 0.304929229735402406490728633,
	 0.961280485811320641748659653, 0.275571819310958163076425168,
	-0.275571819310958163076425168, 0.961280485811320641748659653,
	 0.484869248000791101822951699, 0.874586652278176112634431897,
	-0.874586652278176112634431897, 0.484869248000791101822951699,
	 0.782650596166575738458949301, 0.622461279374149972519166721,
	-0.622461279374149972519166721, 0.782650596166575738458949301,
	 0.113270952177564349018228733, 0.993564135520595333782021697,
	-0.993564135520595333782021697, 0.113270952177564349018228733,
	 0.983662419211730274396237776, 0.180022901405699522679906590,
	-0.180022901405699522679906590, 0.983662419211730274396237776,
	 0.568258952670131549790548489, 0.822849781375826332046780034,
	-0.822849781375826332046780034, 0.568258952670131549790548489,
	 0.839893794195999504583383987, 0.542750784864515906586768661,
	-0.542750784864515906586768661, 0.839893794195999504583383987,
	 0.210111836880469621717489972, 0.977677357824509979943404762,
	-0.977677357824509979943404762, 0.210111836880469621717489972,
	 0.929640895843181265457918066, 0.368466829953372331712746222,
	-0.368466829953372331712746222, 0.929640895843181265457918066,
	 0.396809987416710328595290911, 0.917900775621390457642276297,
	-0.917900775621390457642276297, 0.396809987416710328595290911,
	 0.717870045055731736211325329, 0.696177131491462944788582591,
	-0.696177131491462944788582591, 0.717870045055731736211325329,
	 0.015339206284988101044151868, 0.999882347454212525633049627,
	-0.999882347454212525633049627, 0.015339206284988101044151868,
	 0.999769405351215321657617036, 0.021474080275469507418374898,
	-0.021474080275469507418374898, 0.999769405351215321657617036,
	 0.691759258364157774906734132, 0.722128193929215321243607198,
	-0.722128193929215321243607198, 0.691759258364157774906734132,
	 0.915448716088267819566431292, 0.402434650859418441082533934,
	-0.402434650859418441082533934, 0.915448716088267819566431292,
	 0.362755724367397216204854462, 0.931884265581668106718557199,
	-0.931884265581668106718557199, 0.362755724367397216204854462,
	 0.976369731330021149312732194, 0.216106797076219509948385131,
	-0.216106797076219509948385131, 0.976369731330021149312732194,
	 0.537587076295645482502214932, 0.843208239641845437161743865,
	-0.843208239641845437161743865, 0.537587076295645482502214932,
	 0.819347520076796960824689637, 0.573297166698042212820171239,
	-0.573297166698042212820171239, 0.819347520076796960824689637,
	 0.173983873387463827950700807, 0.984748501801904218556553176,
	-0.984748501801904218556553176, 0.173983873387463827950700807,
	 0.992850414459865090793563344, 0.119365214810991364593637790,
	-0.119365214810991364593637790, 0.992850414459865090793563344,
	 0.617647307937803932403979402, 0.786455213599085757522319464,
	-0.786455213599085757522319464, 0.617647307937803932403979402,
	 0.871595086655951034842481435, 0.490226483288291154229598449,
	-0.490226483288291154229598449, 0.871595086655951034842481435,
	 0.269668325572915106525464462, 0.962953266873683886347921481,
	-0.962953266873683886347921481, 0.269668325572915106525464462,
	 0.950486073949481721759926101, 0.310767152749611495835997250,
	-0.310767152749611495835997250, 0.950486073949481721759926101,
	 0.452349587233770874133026703, 0.891840709392342727796478697,
	-0.891840709392342727796478697, 0.452349587233770874133026703,
	 0.759209188978388033485525443, 0.650846684996380915068975573,
	-0.650846684996380915068975573, 0.759209188978388033485525443,
	 0.076623861392031492278332463, 0.997060070339482978987989949,
	-0.997060070339482978987989949, 0.076623861392031492278332463,
	 0.997511456140303459699448390, 0.070504573389613863027351471,
	-0.070504573389613863027351471, 0.997511456140303459699448390,
	 0.655492852999615385312679701, 0.755201376896536527598710756,
	-0.755201376896536527598710756, 0.655492852999615385312679701,
	 0.894599485631382678433072126, 0.446868840162374195353044389,
	-0.446868840162374195353044389, 0.894599485631382678433072126,
	 0.316593375556165867243047035, 0.948561349915730288158494826,
	-0.948561349915730288158494826, 0.316593375556165867243047035,
	 0.964589793289812723836432159, 0.263754678974831383611349322,
	-0.263754678974831383611349322, 0.964589793289812723836432159,
	 0.495565261825772531150266670, 0.868570705971340895340449876,
	-0.868570705971340895340449876, 0.495565261825772531150266670,
	 0.790230221437310055030217152, 0.612810082429409703935211936,
	-0.612810082429409703935211936, 0.790230221437310055030217152,
	 0.125454983411546238542336453, 0.992099313142191757112085445,
	-0.992099313142191757112085445, 0.125454983411546238542336453,
	 0.985797509167567424700995000, 0.167938294974731178054745536,
	-0.167938294974731178054745536, 0.985797509167567424700995000,
	 0.578313796411655563342245019, 0.815814410806733789010772660,
	-0.815814410806733789010772660, 0.578313796411655563342245019,
	 0.846490938774052078300544488, 0.532403127877197971442805218,
	-0.532403127877197971442805218, 0.846490938774052078300544488,
	 0.222093620973203534094094721, 0.975025345066994146844913468,
	-0.975025345066994146844913468, 0.222093620973203534094094721,
	 0.934092550404258914729877883, 0.357030961233430032614954036,
	-0.357030961233430032614954036, 0.934092550404258914729877883,
	 0.408044162864978680820747499, 0.912962190428398164628018233,
	-0.912962190428398164628018233, 0.408044162864978680820747499,
	 0.726359155084345976817494315, 0.687315340891759108199186948,
	-0.687315340891759108199186948, 0.726359155084345976817494315,
	 0.027608145778965741612354872, 0.999618822495178597116830637,
	-0.999618822495178597116830637, 0.027608145778965741612354872,
	 0.998941293186856850633930266, 0.046003182130914628814301788,
	-0.046003182130914628814301788, 0.998941293186856850633930266,
	 0.673829000378756060917568372, 0.738887324460615147933116508,
	-0.738887324460615147933116508, 0.673829000378756060917568372,
	 0.905296759318118774354048329, 0.424779681209108833357226189,
	-0.424779681209108833357226189, 0.905296759318118774354048329,
	 0.339776884406826857828825803, 0.940506070593268323787291309,
	-0.940506070593268323787291309, 0.339776884406826857828825803,
	 0.970772140728950302138169611, 0.240003022448741486568922365,
	-0.240003022448741486568922365, 0.970772140728950302138169611,
	 0.516731799017649881508753876, 0.856147328375194481019630732,
	-0.856147328375194481019630732, 0.516731799017649881508753876,
	 0.805031331142963597922659282, 0.593232295039799808047809426,
	-0.593232295039799808047809426, 0.805031331142963597922659282,
	 0.149764534677321517229695737, 0.988721691960323767604516485,
	-0.988721691960323767604516485, 0.149764534677321517229695737,
	 0.989622017463200834623694454, 0.143695033150294454819773349,
	-0.143695033150294454819773349, 0.989622017463200834623694454,
	 0.598160706996342311724958652, 0.801376171723140219430247777,
	-0.801376171723140219430247777, 0.598160706996342311724958652,
	 0.859301818357008404783582139, 0.511468850437970399504391001,
	-0.511468850437970399504391001, 0.859301818357008404783582139,
	 0.245955050335794611599924709, 0.969281235356548486048290738,
	-0.969281235356548486048290738, 0.245955050335794611599924709,
	 0.942573197601446879280758735, 0.333999651442009404650865481,
	-0.333999651442009404650865481, 0.942573197601446879280758735,
	 0.430326481340082633908199031, 0.902673318237258806751502391,
	-0.902673318237258806751502391, 0.430326481340082633908199031,
	 0.743007952135121693517362293, 0.669282588346636065720696366,
	-0.669282588346636065720696366, 0.743007952135121693517362293,
	 0.052131704680283321236358216, 0.998640218180265222418199049,
	-0.998640218180265222418199049, 0.052131704680283321236358216,
	 0.995480755491926941769171600, 0.094963495329638998938034312,
	-0.094963495329638998938034312, 0.995480755491926941769171600,
	 0.636761861236284230413943435, 0.771060524261813773200605759,
	-0.771060524261813773200605759, 0.636761861236284230413943435,
	 0.883363338665731594736308015, 0.468688822035827933697617870,
	-0.468688822035827933697617870, 0.883363338665731594736308015,
	 0.293219162694258650606608599, 0.956045251349996443270479823,
	-0.956045251349996443270479823, 0.293219162694258650606608599,
	 0.957826413027532890321037029, 0.287347459544729526477331841,
	-0.287347459544729526477331841, 0.957826413027532890321037029,
	 0.474100214650550014398580015, 0.880470889052160770806542929,
	-0.880470889052160770806542929, 0.474100214650550014398580015,
	 0.774953106594873878359129282, 0.632018735939809021909403706,
	-0.632018735939809021909403706, 0.774953106594873878359129282,
	 0.101069862754827824987887585, 0.994879330794805620591166107,
	-0.994879330794805620591166107, 0.101069862754827824987887585,
	 0.981379193313754574318224190, 0.192080397049892441679288205,
	-0.192080397049892441679288205, 0.981379193313754574318224190,
	 0.558118531220556115693702964, 0.829761233794523042469023765,
	-0.829761233794523042469023765, 0.558118531220556115693702964,
	 0.833170164701913186439915922, 0.553016705580027531764226988,
	-0.553016705580027531764226988, 0.833170164701913186439915922,
	 0.198098410717953586179324918, 0.980182135968117392690210009,
	-0.980182135968117392690210009, 0.198098410717953586179324918,
	 0.925049240782677590302371869, 0.379847208924051170576281147,
	-0.379847208924051170576281147, 0.925049240782677590302371869,
	 0.385516053843918864075607949, 0.922701128333878570437264227,
	-0.922701128333878570437264227, 0.385516053843918864075607949,
	 0.709272826438865651316533772, 0.704934080375904908852523758,
	-0.704934080375904908852523758, 0.709272826438865651316533772,
	 0.003067956762965976270145365, 0.999995293809576171511580126,
	-0.999995293809576171511580126, 0.003067956762965976270145365
};

void
float_FFT(FT *f, unsigned logn)
{
	/*
	 * FFT algorithm in bit-reversal order uses the following
	 * iterative algorithm:
	 *
	 *   t = N
	 *   for m = 1; m < N; m *= 2:
	 *       ht = t/2
	 *       for i1 = 0; i1 < m; i1 ++:
	 *           j1 = i1 * t
	 *           s = GM[m + i1]
	 *           for j = j1; j < (j1 + ht); j ++:
	 *               x = f[j]
	 *               y = s * f[j + ht]
	 *               f[j] = x + y
	 *               f[j + ht] = x - y
	 *       t = ht
	 *
	 * GM[k] contains w^rev(k) for primitive root w = exp(i*pi/N).
	 *
	 * In the description above, f[] is supposed to contain complex
	 * numbers. In our in-memory representation, the real and
	 * imaginary parts of f[k] are in array slots k and k+N/2.
	 *
	 * We only keep the first half of the complex numbers. We can
	 * see that after the first iteration, the first and second halves
	 * of the array of complex numbers have separate lives, so we
	 * simply ignore the second part.
	 */

	unsigned u;
	size_t t, n, hn, m;

	/*
	 * First iteration: compute f[j] + i * f[j+N/2] for all j < N/2
	 * (because GM[1] = w^rev(1) = w^(N/2) = i).
	 * In our chosen representation, this is a no-op: everything is
	 * already where it should be.
	 */

	/*
	 * Subsequent iterations are truncated to use only the first
	 * half of values.
	 */
	n = MKN(logn);
	hn = n >> 1;
	t = hn;
	for (u = 1, m = 2; u < logn; u ++, m <<= 1) {
		size_t ht, hm, i1, j1;

		ht = t >> 1;
		hm = m >> 1;
		for (i1 = 0, j1 = 0; i1 < hm; i1 ++, j1 += t) {
			size_t j, j2;

			j2 = j1 + ht;
			FT s_re, s_im;

			s_re = vrfy_fpr_gm_tab[((m + i1) << 1) + 0];
			s_im = vrfy_fpr_gm_tab[((m + i1) << 1) + 1];
			for (j = j1; j < j2; j ++) {
				FT x_re, x_im, y_re, y_im, tmp;

				x_re = f[j];
				x_im = f[j + hn];
				y_re = f[j + ht];
				y_im = f[j + ht + hn];

				tmp = y_re * s_re - y_im * s_im;
				y_im = y_re * s_im + y_im * s_re;
				y_re = tmp;

				f[j      ] = x_re + y_re;
				f[j   +hn] = x_im + y_im;
				f[j+ht   ] = x_re - y_re;
				f[j+ht+hn] = x_im - y_im;
			}
		}
		t = ht;
	}
}

void
float_iFFT(FT *f, unsigned logn)
{
	/*
	 * Inverse FFT algorithm in bit-reversal order uses the following
	 * iterative algorithm:
	 *
	 *   t = 1
	 *   for m = N; m > 1; m /= 2:
	 *       hm = m/2
	 *       dt = t*2
	 *       for i1 = 0; i1 < hm; i1 ++:
	 *           j1 = i1 * dt
	 *           s = iGM[hm + i1]
	 *           for j = j1; j < (j1 + t); j ++:
	 *               x = f[j]
	 *               y = f[j + t]
	 *               f[j] = x + y
	 *               f[j + t] = s * (x - y)
	 *       t = dt
	 *   for i1 = 0; i1 < N; i1 ++:
	 *       f[i1] = f[i1] / N
	 *
	 * iGM[k] contains (1/w)^rev(k) for primitive root w = exp(i*pi/N)
	 * (actually, iGM[k] = 1/GM[k] = conj(GM[k])).
	 *
	 * In the main loop (not counting the final division loop), in
	 * all iterations except the last, the first and second half of f[]
	 * (as an array of complex numbers) are separate. In our chosen
	 * representation, we do not keep the second half.
	 *
	 * The last iteration recombines the recomputed half with the
	 * implicit half, and should yield only real numbers since the
	 * target polynomial is real; moreover, s = i at that step.
	 * Thus, when considering x and y:
	 *    y = conj(x) since the final f[j] must be real
	 *    Therefore, f[j] is filled with 2*Re(x), and f[j + t] is
	 *    filled with 2*Im(x).
	 * But we already have Re(x) and Im(x) in array slots j and j+t
	 * in our chosen representation. That last iteration is thus a
	 * simple doubling of the values in all the array.
	 *
	 * We make the last iteration a no-op by tweaking the final
	 * division into a division by N/2, not N.
	 */
	size_t u, n, hn, t, m;

	n = MKN(logn);
	t = 1;
	m = n;
	hn = n >> 1;
	for (u = logn; u > 1; u --) {
		size_t hm, dt, i1, j1;

		hm = m >> 1;
		dt = t << 1;
		for (i1 = 0, j1 = 0; j1 < hn; i1 ++, j1 += dt) {
			size_t j, j2;

			j2 = j1 + t;
			FT s_re, s_im;

			s_re = vrfy_fpr_gm_tab[((hm + i1) << 1) + 0];
			s_im = -vrfy_fpr_gm_tab[((hm + i1) << 1) + 1];
			for (j = j1; j < j2; j ++) {
				FT x_re, x_im, y_re, y_im, tmp;

				x_re = f[j];
				x_im = f[j + hn];
				y_re = f[j + t];
				y_im = f[j + t + hn];

				f[j     ] = x_re + y_re;
				f[j + hn] = x_im + y_im;

				x_re -= y_re;
				x_im -= y_im;

				tmp = x_re * s_re - x_im * s_im;
				f[j + t + hn] = x_re * s_im + x_im * s_re;
				f[j + t     ] = tmp;
			}
		}
		t = dt;
		m = hm;
	}

	/*
	 * Last iteration is a no-op, provided that we divide by N/2
	 * instead of N. We need to make a special case for logn = 0.
	 */
	if (logn > 0) {
		FT ni = 2.0 / n;
		for (u = 0; u < n; u ++) {
			f[u] *= ni;
		}
	}
}


/* see inner.h */
int Zf(verify_NTT)(const uint8_t *restrict h,
	const int16_t *restrict s1,
	const int16_t *restrict q00, const int16_t *restrict q10,
	unsigned logn, uint8_t *restrict tmp)
{
	/*
	 * Recover s0 with s0 = round(h0 / 2 + (h1 / 2 - s1) q10 / q00)
	 * = round((h0 q00 + (h1 - 2 s1) q10) / (2q00)).
	 * Put (t0, t1) = (h0 - 2 * s0, h1 - 2 * s1) in FFT representation.
	 *
	 * Use floats only here...
	 */

	size_t n, u, v, w;
	int16_t *s0;
	uint32_t p, p0i, R, R2, *gm, *igm, *q10_i32, *s1_i32;
	FT *t0, *t1;
	uint8_t h0;

	n = MKN(logn);

	t0 = (FT *)tmp;
	t1 = (FT *)(tmp + 8*n);
	// t1 = t0 + 2 * n; // To have same alignment as if it would be an fpr (double)

	gm = (uint32_t *)(tmp + 16*n);
	igm = gm + n;
	s0 = (int16_t *)(igm + n);
	q10_i32 = igm + n;
	s1_i32 = ((uint32_t *)t1) + n;

	p = PRIMES[0].p;
	p0i = modp_ninv31(p);
	R = modp_R(p);
	R2 = modp_R2(p, p0i);
	modp_mkgm2(gm, igm, logn, PRIMES[0].g, p, p0i);

	int16_to_ntt(q10_i32, q10, gm, p, p0i, logn);
	hash_to_ntt(s1_i32, SECOND_HASH(h, logn), s1, gm, R2, p, p0i, logn);
	for (u = 0; u < n; u++) {
		s1_i32[u] = modp_montymul(s1_i32[u], q10_i32[u], p, p0i);
	}

	/*
	 * Convert (h1 - 2 s1) q10, calculated in NTT, to FFT format, so we can do
	 * division and rounding.
	 */
	modp_iNTT2(s1_i32, igm, logn, p, p0i);
	for (u = 0; u < n; u++) {
		t0[u] = (FT) modp_norm(s1_i32[u], p);
		t1[u] = (FT) q00[u];
	}
	float_FFT(t0, logn);
	float_FFT(t1, logn);

	for (u = 0; u < n/2; u++) {
		t0[u] /= t1[u];
		t0[u + n/2] /= t1[u];
	}

	float_iFFT(t0, logn);

	/*
	 * verify_norm_NTT expects that q10 is in NTT representation in the second
	 * half of t1.
	 */
	memcpy(s1_i32, q10_i32, n * sizeof *q10_i32);

	/*
	 * Recover s0 with s0 = round(h0 / 2 + (h1 / 2 - s1) q10 / q00).
	 * Put (t0, t1) = (h0 - 2 * s0, h1 - 2 * s1) in FFT representation.
	 */
	if (logn <= 3) {
		h0 = h[0];
		for (u = 0; u < n; u ++) {
			s0[u] = lround(((h0 & 1) + t0[u]) / 2.0);
			h0 >>= 1;
		}
	} else {
		for (u = 0, w = 0; w < n; u ++) {
			h0 = h[u];
			for (v = 0; v < 8; v ++, w ++) {
				s0[w] = lround(((h0 & 1) + t0[w]) / 2.0);
				h0 >>= 1;
			}
		}
	}

	return verify_norm_NTT(h, s0, s1, q00, q10, p, p0i, R, R2, logn, tmp);
}
