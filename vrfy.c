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

/* fpr.h for 32 bit */
/*
 * Custom floating-point implementation with integer arithmetics. We
 * use IEEE-754 "binary32" format (float), with some simplifications:
 *
 *   - Top bit is s = 1 for negative, 0 for positive.
 *
 *   - Exponent e uses the next 8 bits (bits 23 to 30, inclusive).
 *
 *   - Mantissa m uses the 23 low bits.
 *
 * Encoded value is, in general: (-1)^s * 2^(e-127) * (1 + m*2^(-23))
 * i.e. the mantissa really is a 24-bit number (less than 2.0, but not
 * less than 1.0), but the top bit (equal to 1 by definition) is omitted
 * in the encoding.
 */
typedef uint32_t fpr32;

/*
 * Normalize a provided unsigned integer to the 2^31..2^32-1 range by
 * left-shifting it if necessary. The exponent e is adjusted accordingly
 * (i.e. if the value was left-shifted by n bits, then n is subtracted
 * from e). If source m is 0, then it remains 0, but e is altered.
 * Both m and e must be simple variables (no expressions allowed).
 */
#define FPR_NORM32(m, e)   do { \
		uint32_t nt; \
 \
		(e) -= 31; \
 \
		nt = (uint32_t)((m) >> 16); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) << 16)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 4); \
 \
		nt = (uint32_t)((m) >> 24); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) <<  8)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 3); \
 \
		nt = (uint32_t)((m) >> 28); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) <<  4)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 2); \
 \
		nt = (uint32_t)((m) >> 30); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) <<  2)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 1); \
 \
		nt = (uint32_t)((m) >> 31); \
		(m) ^= ((m) ^ ((m) <<  1)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt); \
	} while (0)

/*
 * For computations, we split values into an integral mantissa in the
 * 2^25..2^26 range, and an (adjusted) exponent. The lowest bit is
 * "sticky" (it is set to 1 if any of the bits below it is 1); when
 * re-encoding, the low two bits are dropped, but may induce an
 * increment in the value for proper rounding.
 */

/*
 * Expectations:
 *   s = 0 or 1
 *   exponent e is "arbitrary" and unbiased
 *   2^25 <= m < 2^26
 * Numerical value is (-1)^2 * m * 2^e
 *
 * Exponents which are too low lead to value zero. If the exponent is
 * too large, the returned value is indeterminate.
 *
 * If m = 0, then a zero is returned (using the provided sign).
 * If e < -151, then a zero is returned (regardless of the value of m).
 * If e >= -151 and e != 0, m must be within the expected range
 * (2^25 to 2^26-1).
 */
static inline fpr32
FPR32(int s, int e, uint32_t m)
{
	fpr32 x;
	uint32_t t;
	unsigned f;

	/*
	 * If e >= -151, then the value is "normal"; otherwise, it
	 * should be a subnormal, which we clamp down to zero.
	 */
	e += 151;
	t = (uint32_t)e >> 31;
	m &= (uint32_t)t - 1;

	/*
	 * If m = 0 then we want a zero; make e = 0 too, but conserve
	 * the sign.
	 */
	t = (uint32_t)(m >> 25);
	e &= -(int)t;

	/*
	 * The 23 mantissa bits come from m. Value m has its top bit set
	 * (unless it is a zero); we leave it "as is": the top bit will
	 * increment the exponent by 1, except when m = 0, which is
	 * exactly what we want.
	 */
	x = (((uint32_t)s << 31) | (m >> 2)) + ((uint32_t)e << 23);

	/*
	 * Rounding: if the low three bits of m are 011, 110 or 111,
	 * then the value should be incremented to get the next
	 * representable value. This implements the usual
	 * round-to-nearest rule (with preference to even values in case
	 * of a tie). Note that the increment may make a carry spill
	 * into the exponent field, which is again exactly what we want
	 * in that case.
	 */
	f = (unsigned)m & 7U;
	x += (0xC8U >> f) & 1;
	return x;
}

// #define fpr32_scaled   Zf(fpr32_scaled)
fpr32 fpr32_scaled(int32_t i, int sc);

static inline fpr32
fpr32_of(int64_t i)
{
	return fpr32_scaled(i, 0);
}

/*
static const fpr32 fpr32_zero = 0;
static const fpr32 fpr32_one = 1065353216; // 127 << 23
static const fpr32 fpr32_two = 1073741824; // 1 << 30
static const fpr32 fpr32_onehalf = 1056964608; // 126 << 23
*/

static inline int32_t
fpr32_rint(fpr32 x)
{
	uint32_t m, d, dd, s, f;
	int e;

	/*
	 * We assume that the value fits in -(2^31-1)..+(2^31-1). We can
	 * thus extract the mantissa as a 31-bit integer, then right-shift
	 * it as needed.
	 */
	m = ((x << 7) | ((uint32_t)1 << 30)) & (((uint32_t)1 << 31) - 1);
	e = 157 - ((int)(x >> 23) & 0xFF);

	/*
	 * If a shift of more than 31 bits is needed, then simply set m
	 * to zero. This also covers the case of an input operand equal
	 * to zero.
	 */
	m &= -((uint32_t)(e - 32) >> 31);
	e &= 31;

	/*
	 * Right-shift m as needed. Shift count is e. Proper rounding
	 * mandates that:
	 *   - If the highest dropped bit is zero, then round low.
	 *   - If the highest dropped bit is one, and at least one of the
	 *     other dropped bits is one, then round up.
	 *   - If the highest dropped bit is one, and all other dropped
	 *     bits are zero, then round up if the lowest kept bit is 1,
	 *     or low otherwise (i.e. ties are broken by "rounding to even").
	 *
	 * We thus first extract a word consisting of all the dropped bit
	 * AND the lowest kept bit; then we shrink it down to three bits,
	 * the lowest being "sticky".
	 */
	d = m << (31 - e);
	/*
	 * dd != 0 iff there is one of the lowest 29 bits of d that is a '1'.
	 * In this, set f = f | 1;
	 */
	dd = d & 0x1FFFFFF;
	f = (uint32_t)(d >> 29) | ((dd | -dd) >> 31);
	m = (m >> e) + (uint32_t)((0xC8U >> f) & 1U);

	/*
	 * Apply the sign bit.
	 */
	s = (uint32_t)(x >> 31);
	return ((int32_t)m ^ -(int32_t)s) + (int32_t)s;
}

#define fpr32_add   Zf(fpr32_add)
fpr32 fpr32_add(fpr32 x, fpr32 y);

static inline fpr32
fpr32_neg(fpr32 x)
{
	return x ^ (uint32_t)1 << 31;
}

static inline fpr32
fpr32_sub(fpr32 x, fpr32 y)
{
	return fpr32_add(x, fpr32_neg(y));
}

static inline fpr32
fpr32_half(fpr32 x)
{
	/*
	 * To divide a value by 2, we just have to subtract 1 from its
	 * exponent, but we have to take care of zero.
	 */
	uint32_t t;

	x -= (uint32_t)1 << 23;
	t = (((uint32_t)(x >> 23) & 0xFF) + 1) >> 8;
	x &= (uint32_t)t - 1;
	return x;
}

#define fpr32_div   Zf(fpr32_div)
fpr32 fpr32_div(fpr32 x, fpr32 y);

/* fpr.c */
/*
 * Normalize a provided unsigned integer to the 2^31..2^32-1 range by
 * left-shifting it if necessary. The exponent e is adjusted accordingly
 * (i.e. if the value was left-shifted by n bits, then n is subtracted
 * from e). If source m is 0, then it remains 0, but e is altered.
 * Both m and e must be simple variables (no expressions allowed).
 */
#define FPR_NORM32(m, e)   do { \
		uint32_t nt; \
 \
		(e) -= 31; \
 \
		nt = (uint32_t)((m) >> 16); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) << 16)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 4); \
 \
		nt = (uint32_t)((m) >> 24); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) <<  8)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 3); \
 \
		nt = (uint32_t)((m) >> 28); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) <<  4)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 2); \
 \
		nt = (uint32_t)((m) >> 30); \
		nt = (nt | -nt) >> 31; \
		(m) ^= ((m) ^ ((m) <<  2)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt << 1); \
 \
		nt = (uint32_t)((m) >> 31); \
		(m) ^= ((m) ^ ((m) <<  1)) & ((uint64_t)nt - 1); \
		(e) += (int)(nt); \
	} while (0)


fpr32
fpr32_scaled(int32_t i, int sc)
{
	/*
	 * To convert from int to float, we have to do the following:
	 *  1. Get the absolute value of the input, and its sign
	 *  2. Shift right or left the value as appropriate
	 *  3. Pack the result
	 *
	 * We can assume that the source integer is not -2^31.
	 */
	int s, e;
	uint32_t t, m;

	/*
	 * Extract sign bit.
	 * We have: -i = 1 + ~i
	 */
	s = (int)((uint32_t)i >> 31);
	i ^= -(int32_t)s;
	i += s;

	/*
	 * For now we suppose that i != 0.
	 * Otherwise, we set m to i and left-shift it as much as needed
	 * to get a 1 in the top bit. We can do that in a logarithmic
	 * number of conditional shifts.
	 */
	m = (uint32_t)i;
	e = 6 + sc;
	FPR_NORM32(m, e);

	/*
	 * Now m is in the 2^31..2^32-1 range. We must divide it by 64;
	 * if one of the dropped bits is a 1, this should go into the
	 * "sticky bit".
	 */
	m |= ((uint32_t)m & 0x3F) + 0x3F;
	m >>= 6;

	/*
	 * Corrective action: if i = 0 then all of the above was
	 * incorrect, and we clamp e and m down to zero.
	 */
	t = (uint32_t)((uint32_t)(i | -i) >> 31);
	m &= -t;
	e &= -(int)t;

	/*
	 * Assemble back everything. The FPR32() function will handle cases
	 * where e is too low.
	 */
	return FPR32(s, e, m);
}

fpr32
fpr32_add(fpr32 x, fpr32 y)
{
	uint32_t m, xu, yu, za, cs;
	int ex, ey, sx, sy, cc;

	/*
	 * Make sure that the first operand (x) has the larger absolute
	 * value. This guarantees that the exponent of y is less than
	 * or equal to the exponent of x, and, if they are equal, then
	 * the mantissa of y will not be greater than the mantissa of x.
	 *
	 * After this swap, the result will have the sign x, except in
	 * the following edge case: abs(x) = abs(y), and x and y have
	 * opposite sign bits; in that case, the result shall be +0
	 * even if the sign bit of x is 1. To handle this case properly,
	 * we do the swap if abs(x) = abs(y) AND the sign of x is 1.
	 */
	m = ((uint32_t)1 << 31) - 1;
	za = (x & m) - (y & m);
	cs = (uint32_t)(za >> 31)
		| ((1U - (uint32_t)(-za >> 31)) & (uint32_t)(x >> 31));
	m = (x ^ y) & -(uint32_t)cs;
	x ^= m;
	y ^= m;

	/*
	 * Extract sign bits, exponents and mantissas. The mantissas are
	 * scaled up to 2^26..2^27-1, and the exponent is unbiased. If
	 * an operand is zero, its mantissa is set to 0 at this step, and
	 * its exponent will be -153.
	 */
	ex = (int)(x >> 23);
	sx = ex >> 8;
	ex &= 0xFF;
	m = (uint32_t)((ex + 0xFF) >> 8) << 23;
	xu = ((x & (((uint32_t)1 << 23) - 1)) | m) << 3;
	ex -= 153;
	ey = (int)(y >> 23);
	sy = ey >> 8;
	ey &= 0xFF;
	m = (uint32_t)((ey + 0xFF) >> 8) << 23;
	yu = ((y & (((uint32_t)1 << 23) - 1)) | m) << 3;
	ey -= 153;

	/*
	 * x has the larger exponent; hence, we only need to right-shift y.
	 * If the shift count is larger than 27 bits then we clamp the
	 * value to zero.
	 */
	cc = ex - ey;
	yu &= -((uint32_t)(cc - 28) >> 31);
	cc &= 31;

	/*
	 * The lowest bit of yu is "sticky".
	 */
	m = (1u << cc) - 1;
	yu |= (yu & m) + m;
	yu = yu >> cc;

	/*
	 * If the operands have the same sign, then we add the mantissas;
	 * otherwise, we subtract the mantissas.
	 */
	xu += yu - ((yu << 1) & -(uint32_t)(sx ^ sy));

	/*
	 * The result may be smaller, or slightly larger. We normalize
	 * it to the 2^31..2^32-1 range (if xu is zero, then it stays
	 * at zero).
	 */
	FPR_NORM32(xu, ex);

	/*
	 * Scale down the value to 2^25..2^26-1, handling the last bit
	 * as sticky.
	 */
	xu |= ((uint32_t)xu & 0x3F) + 0x3F;
	xu >>= 6;
	ex += 6;

	/*
	 * In general, the result has the sign of x. However, if the
	 * result is exactly zero, then the following situations may
	 * be encountered:
	 *   x > 0, y = -x   -> result should be +0
	 *   x < 0, y = -x   -> result should be +0
	 *   x = +0, y = +0  -> result should be +0
	 *   x = -0, y = +0  -> result should be +0
	 *   x = +0, y = -0  -> result should be +0
	 *   x = -0, y = -0  -> result should be -0
	 *
	 * But at the conditional swap step at the start of the
	 * function, we ensured that if abs(x) = abs(y) and the
	 * sign of x was 1, then x and y were swapped. Thus, the
	 * two following cases cannot actually happen:
	 *   x < 0, y = -x
	 *   x = -0, y = +0
	 * In all other cases, the sign bit of x is conserved, which
	 * is what the FPR32() function does. The FPR32() function also
	 * properly clamps values to zero when the exponent is too
	 * low, but does not alter the sign in that case.
	 */
	return FPR32(sx, ex, xu);
}



fpr32
fpr32_mul(fpr32 x, fpr32 y)
{
	uint32_t xu, yu, zu, zv, w;
	uint64_t z;
	int ex, ey, d, e, s;

	/*
	 * Extract absolute values as scaled unsigned integers. We
	 * don't extract exponents yet.
	 */
	xu = (x & (((uint32_t)1 << 23) - 1)) | ((uint32_t)1 << 23);
	yu = (y & (((uint32_t)1 << 23) - 1)) | ((uint32_t)1 << 23);

	/*
	 * Multiply the two 24-bit integers.
	 */
	z = (uint64_t)xu * (uint64_t)yu;

	/*
	 * Since xu and yu are both in the 2^23..2^24-1 range, the
	 * product is in the 2^46..2^48-1 range. We first reassemble
	 * it and round it into the 2^25..2^27-1 range; the bottom bit
	 * is made "sticky".
	 */
	// zu = (uint32_t)((z + ((uint32_t)1 << 21) - 1)  >> 21);
	zu = (uint32_t)(z >> 21);
	zv = (uint32_t)z & 0x1FFFFF;
	zu |= (zv + 0x1FFFFF) >> 21;

	/*
	 * We normalize zu to the 2^25..2^26-1 range: it could be one
	 * bit too large at this point. This is done with a conditional
	 * right-shift that takes into account the sticky bit.
	 */
	zv = (zu >> 1) | (zu & 1);
	w = zu >> 26;
	zu ^= (zu ^ zv) & -w;

	/*
	 * Get the aggregate scaling factor:
	 *
	 *   - Each exponent is biased by 127.
	 *
	 *   - Integral mantissas are scaled by 2^23, hence an
	 *     extra 23 bias for each exponent.
	 *
	 *   - However, we right-shifted z by 21 bits, and then
	 *     by 0 or 1 extra bit (depending on the value of w).
	 *
	 * In total, we must add the exponents, then subtract
	 * 2 * (127 + 23), then add 21 + w.
	 */
	ex = (int)((x >> 23) & 0xFF);
	ey = (int)((y >> 23) & 0xFF);
	e = ex + ey - 279 + (int)w;

	/*
	 * Sign bit is the XOR of the operand sign bits.
	 */
	s = (int)((x ^ y) >> 31);

	/*
	 * Corrective actions for zeros: if either of the operands is
	 * zero, then the computations above were wrong. Test for zero
	 * is whether ex or ey is zero. We just have to set the mantissa
	 * (zu) to zero, the FPR32() function will normalize e.
	 */
	d = ((ex + 0xFF) & (ey + 0xFF)) >> 8;
	zu &= (uint32_t) -d;

	/*
	 * FPR32() packs the result and applies proper rounding.
	 */
	return FPR32(s, e, zu);
}

fpr32
fpr32_div(fpr32 x, fpr32 y)
{
	uint32_t xu, yu, q, q2, w;
	int i, ex, ey, e, d, s;

	/*
	 * Extract mantissas of x and y (unsigned).
	 */
	xu = (x & (((uint32_t)1 << 23) - 1)) | ((uint32_t)1 << 23);
	yu = (y & (((uint32_t)1 << 23) - 1)) | ((uint32_t)1 << 23);

	/*
	 * Perform bit-by-bit division of xu by yu. We run it for 26 bits.
	 */
	q = 0;
	for (i = 0; i < 26; i ++) {
		/*
		 * If yu is less than or equal xu, then subtract it and
		 * push a 1 in the quotient; otherwise, leave xu unchanged
		 * and push a 0.
		 */
		uint32_t b;

		b = ((xu - yu) >> 31) - 1;
		xu -= b & yu;
		q |= b & 1;
		xu <<= 1;
		q <<= 1;
	}

	/*
	 * We got 26 bits in the quotient, followed by an extra zero. We
	 * want that 27th bit to be "sticky": it should be a 1 if and
	 * only if the remainder (xu) is non-zero.
	 */
	q |= (xu | -xu) >> 31;

	/*
	 * Quotient is at most 2^27-1. Its top bit may be zero, but in
	 * that case the next-to-top bit will be a one, since the
	 * initial xu and yu were both in the 2^23..2^24-1 range.
	 * We perform a conditional shift to normalize q to the
	 * 2^25..2^26-1 range (with the bottom bit being sticky).
	 */
	q2 = (q >> 1) | (q & 1);
	w = q >> 26;
	q ^= (q ^ q2) & -w;

	/*
	 * Extract exponents to compute the scaling factor:
	 *
	 *   - Each exponent is biased and we scaled them up by
	 *     23 bits; but these biases will cancel out.
	 *
	 *   - The division loop produced a 26-bit shifted result,
	 *     so we must scale it down by 26 bits.
	 *
	 *   - If w = 1, we right-shifted the integer by 1 bit,
	 *     hence we must add 1 to the scaling.
	 */
	ex = (int)((x >> 23) & 0xFF);
	ey = (int)((y >> 23) & 0xFF);
	e = ex - ey - 26 + (int)w;

	/*
	 * Sign is the XOR of the signs of the operands.
	 */
	s = (int)((x ^ y) >> 31);

	/*
	 * Corrective actions for zeros: if x = 0, then the computation
	 * is wrong, and we must clamp e and q to 0. We do not care
	 * about the case y = 0 (as per assumptions in this module,
	 * the caller does not perform divisions by zero).
	 */
	d = (ex + 0xFF) >> 8;
	s &= d;
	e &= -d;
	q &= -(uint32_t)d;

	/*
	 * FPR32() packs the result and applies proper rounding.
	 */
	return FPR32(s, e, q);
}

const fpr32 vrfy_fpr_gm_tab[] = {
         0u,          0u,/*unused*/2147483648u, 1065353216u,
		 1060439283u, 1060439283u, 3207922931u, 1060439283u,
		 1064076126u, 1053028117u, 3200511765u, 1064076126u,
		 1053028117u, 1064076126u, 3211559774u, 1053028117u,
		 1065030846u, 1044891074u, 3192374722u, 1065030846u,
		 1057896922u, 1062525745u, 3210009393u, 1057896922u,
		 1062525745u, 1057896922u, 3205380570u, 1062525745u,
		 1044891074u, 1065030846u, 3212514494u, 1044891074u,
		 1065272429u, 1036565814u, 3184049462u, 1065272429u,
		 1059219353u, 1061544963u, 3209028611u, 1059219353u,
		 1063372184u, 1056004842u, 3203488490u, 1063372184u,
		 1049927729u, 1064630795u, 3212114443u, 1049927729u,
		 1064630795u, 1049927729u, 3197411377u, 1064630795u,
		 1056004842u, 1063372184u, 3210855832u, 1056004842u,
		 1061544963u, 1059219353u, 3206703001u, 1061544963u,
		 1036565814u, 1065272429u, 3212756077u, 1036565814u,
		 1065333007u, 1028193072u, 3175676720u, 1065333007u,
		 1059842890u, 1061007097u, 3208490745u, 1059842890u,
		 1063742424u, 1054533760u, 3202017408u, 1063742424u,
		 1051491540u, 1064372488u, 3211856136u, 1051491540u,
		 1064850424u, 1048104908u, 3195588556u, 1064850424u,
		 1057201213u, 1062966298u, 3210449946u, 1057201213u,
		 1062051586u, 1058570176u, 3206053824u, 1062051586u,
		 1041645699u, 1065171628u, 3212655276u, 1041645699u,
		 1065171628u, 1041645699u, 3189129347u, 1065171628u,
		 1058570176u, 1062051586u, 3209535234u, 1058570176u,
		 1062966298u, 1057201213u, 3204684861u, 1062966298u,
		 1048104908u, 1064850424u, 3212334072u, 1048104908u,
		 1064372488u, 1051491540u, 3198975188u, 1064372488u,
		 1054533760u, 1063742424u, 3211226072u, 1054533760u,
		 1061007097u, 1059842890u, 3207326538u, 1061007097u,
		 1028193072u, 1065333007u, 3212816655u, 1028193072u,
		 1065348163u, 1019808432u, 3167292080u, 1065348163u,
		 1060144571u, 1060726850u, 3208210498u, 1060144571u,
		 1063913895u, 1053785034u, 3201268682u, 1063913895u,
		 1052263466u, 1064229022u, 3211712670u, 1052263466u,
		 1064945565u, 1046502419u, 3193986067u, 1064945565u,
		 1057551771u, 1062750291u, 3210233939u, 1057551771u,
		 1062292797u, 1058236458u, 3205720106u, 1062292797u,
		 1043271842u, 1065106216u, 3212589864u, 1043271842u,
		 1065227044u, 1039839859u, 3187323507u, 1065227044u,
		 1058897873u, 1061802258u, 3209285906u, 1058897873u,
		 1063173637u, 1056726311u, 3204209959u, 1063173637u,
		 1049136787u, 1064745479u, 3212229127u, 1049136787u,
		 1064506439u, 1050712805u, 3198196453u, 1064506439u,
		 1055273845u, 1063561817u, 3211045465u, 1055273845u,
		 1061279856u, 1059534422u, 3207018070u, 1061279856u,
		 1033283845u, 1065307757u, 3212791405u, 1033283845u,
		 1065307757u, 1033283845u, 3180767493u, 1065307757u,
		 1059534422u, 1061279856u, 3208763504u, 1059534422u,
		 1063561817u, 1055273845u, 3202757493u, 1063561817u,
		 1050712805u, 1064506439u, 3211990087u, 1050712805u,
		 1064745479u, 1049136787u, 3196620435u, 1064745479u,
		 1056726311u, 1063173637u, 3210657285u, 1056726311u,
		 1061802258u, 1058897873u, 3206381521u, 1061802258u,
		 1039839859u, 1065227044u, 3212710692u, 1039839859u,
		 1065106216u, 1043271842u, 3190755490u, 1065106216u,
		 1058236458u, 1062292797u, 3209776445u, 1058236458u,
		 1062750291u, 1057551771u, 3205035419u, 1062750291u,
		 1046502419u, 1064945565u, 3212429213u, 1046502419u,
		 1064229022u, 1052263466u, 3199747114u, 1064229022u,
		 1053785034u, 1063913895u, 3211397543u, 1053785034u,
		 1060726850u, 1060144571u, 3207628219u, 1060726850u,
		 1019808432u, 1065348163u, 3212831811u, 1019808432u,
		 1065351953u, 1011420816u, 3158904464u, 1065351953u,
		 1060292809u, 1060583971u, 3208067619u, 1060292809u,
		 1063996172u, 1053407571u, 3200891219u, 1063996172u,
		 1052646730u, 1064153747u, 3211637395u, 1052646730u,
		 1064989442u, 1045697793u, 3193181441u, 1064989442u,
		 1057725035u, 1062639077u, 3210122725u, 1057725035u,
		 1062410313u, 1058067405u, 3205551053u, 1062410313u,
		 1044082383u, 1065069773u, 3212553421u, 1044082383u,
		 1065250992u, 1038203950u, 3185687598u, 1065250992u,
		 1059059403u, 1061674597u, 3209158245u, 1059059403u,
		 1063274017u, 1056366795u, 3203850443u, 1063274017u,
		 1049532962u, 1064689350u, 3212172998u, 1049532962u,
		 1064569821u, 1050321030u, 3197804678u, 1064569821u,
		 1055640507u, 1063468122u, 3210951770u, 1055640507u,
		 1061413376u, 1059377701u, 3206861349u, 1061413376u,
		 1034925696u, 1065291352u, 3212775000u, 1034925696u,
		 1065321643u, 1031482228u, 3178965876u, 1065321643u,
		 1059689493u, 1061144423u, 3208628071u, 1059689493u,
		 1063653256u, 1054904911u, 3202388559u, 1063653256u,
		 1051102994u, 1064440658u, 3211924306u, 1051102994u,
		 1064799173u, 1048739264u, 3196222912u, 1064799173u,
		 1057023972u, 1063071059u, 3210554707u, 1057023972u,
		 1061927928u, 1058734790u, 3206218438u, 1061927928u,
		 1040830342u, 1065200588u, 3212684236u, 1040830342u,
		 1065140169u, 1042459574u, 3189943222u, 1065140169u,
		 1058404057u, 1062173215u, 3209656863u, 1058404057u,
		 1062859370u, 1057377154u, 3204860802u, 1062859370u,
		 1047304831u, 1064899224u, 3212382872u, 1047304831u,
		 1064301939u, 1051878383u, 3199362031u, 1064301939u,
		 1054160449u, 1063829308u, 3211312956u, 1054160449u,
		 1060867899u, 1059994590u, 3207478238u, 1060867899u,
		 1024901932u, 1065341847u, 3212825495u, 1024901932u,
		 1065341847u, 1024901932u, 3172385580u, 1065341847u,
		 1059994590u, 1060867899u, 3208351547u, 1059994590u,
		 1063829308u, 1054160449u, 3201644097u, 1063829308u,
		 1051878383u, 1064301939u, 3211785587u, 1051878383u,
		 1064899224u, 1047304831u, 3194788479u, 1064899224u,
		 1057377154u, 1062859370u, 3210343018u, 1057377154u,
		 1062173215u, 1058404057u, 3205887705u, 1062173215u,
		 1042459574u, 1065140169u, 3212623817u, 1042459574u,
		 1065200588u, 1040830342u, 3188313990u, 1065200588u,
		 1058734790u, 1061927928u, 3209411576u, 1058734790u,
		 1063071059u, 1057023972u, 3204507620u, 1063071059u,
		 1048739264u, 1064799173u, 3212282821u, 1048739264u,
		 1064440658u, 1051102994u, 3198586642u, 1064440658u,
		 1054904911u, 1063653256u, 3211136904u, 1054904911u,
		 1061144423u, 1059689493u, 3207173141u, 1061144423u,
		 1031482228u, 1065321643u, 3212805291u, 1031482228u,
		 1065291352u, 1034925696u, 3182409344u, 1065291352u,
		 1059377701u, 1061413376u, 3208897024u, 1059377701u,
		 1063468122u, 1055640507u, 3203124155u, 1063468122u,
		 1050321030u, 1064569821u, 3212053469u, 1050321030u,
		 1064689350u, 1049532962u, 3197016610u, 1064689350u,
		 1056366795u, 1063274017u, 3210757665u, 1056366795u,
		 1061674597u, 1059059403u, 3206543051u, 1061674597u,
		 1038203950u, 1065250992u, 3212734640u, 1038203950u,
		 1065069773u, 1044082383u, 3191566031u, 1065069773u,
		 1058067405u, 1062410313u, 3209893961u, 1058067405u,
		 1062639077u, 1057725035u, 3205208683u, 1062639077u,
		 1045697793u, 1064989442u, 3212473090u, 1045697793u,
		 1064153747u, 1052646730u, 3200130378u, 1064153747u,
		 1053407571u, 1063996172u, 3211479820u, 1053407571u,
		 1060583971u, 1060292809u, 3207776457u, 1060583971u,
		 1011420816u, 1065351953u, 3212835601u, 1011420816u,
		 1065352900u, 1003032456u, 3150516104u, 1065352900u,
		 1060366268u, 1060511852u, 3207995500u, 1060366268u,
		 1064036440u, 1053218089u, 3200701737u, 1064036440u,
		 1052837662u, 1064115229u, 3211598877u, 1052837662u,
		 1065010454u, 1045294688u, 3192778336u, 1065010454u,
		 1057811152u, 1062582675u, 3210066323u, 1057811152u,
		 1062468291u, 1057982340u, 3205465988u, 1062468291u,
		 1044486967u, 1065050620u, 3212534268u, 1044486967u,
		 1065262025u, 1037385145u, 3184868793u, 1065262025u,
		 1059139577u, 1061610026u, 3209093674u, 1059139577u,
		 1063323378u, 1056186119u, 3203669767u, 1063323378u,
		 1049730525u, 1064660375u, 3212144023u, 1049730525u,
		 1064600610u, 1050124567u, 3197608215u, 1064600610u,
		 1055822969u, 1063420432u, 3210904080u, 1055822969u,
		 1061479413u, 1059298729u, 3206782377u, 1061479413u,
		 1035745987u, 1065282205u, 3212765853u, 1035745987u,
		 1065327640u, 1029837929u, 3177321577u, 1065327640u,
		 1059766402u, 1061075995u, 3208559643u, 1059766402u,
		 1063698124u, 1054719609u, 3202203257u, 1063698124u,
		 1051297476u, 1064406871u, 3211890519u, 1051297476u,
		 1064825104u, 1048504033u, 3195987681u, 1064825104u,
		 1057112753u, 1063018951u, 3210502599u, 1057112753u,
		 1061990009u, 1058652672u, 3206136320u, 1061990009u,
		 1041238199u, 1065186420u, 3212670068u, 1041238199u,
		 1065156211u, 1042052830u, 3189536478u, 1065156211u,
		 1058487303u, 1062112656u, 3209596304u, 1058487303u,
		 1062913104u, 1057289348u, 3204772996u, 1062913104u,
		 1047705169u, 1064875131u, 3212358779u, 1047705169u,
		 1064337510u, 1051685178u, 3199168826u, 1064337510u,
		 1054347371u, 1063786152u, 3211269800u, 1054347371u,
		 1060937731u, 1059918953u, 3207402601u, 1060937731u,
		 1026547719u, 1065337743u, 3212821391u, 1026547719u,
		 1065345321u, 1023101370u, 3170585018u, 1065345321u,
		 1060069797u, 1060797604u, 3208281252u, 1060069797u,
		 1063871889u, 1053973001u, 3201456649u, 1063871889u,
		 1052071148u, 1064265776u, 3211749424u, 1052071148u,
		 1064922702u, 1046903910u, 3194387558u, 1064922702u,
		 1057464630u, 1062805098u, 3210288746u, 1057464630u,
		 1062233263u, 1058320441u, 3205804089u, 1062233263u,
		 1042865916u, 1065123504u, 3212607152u, 1042865916u,
		 1065214129u, 1040422146u, 3187905794u, 1065214129u,
		 1058816524u, 1061865343u, 3209348991u, 1058816524u,
		 1063122622u, 1056905138u, 3204388786u, 1063122622u,
		 1048938190u, 1064772631u, 3212256279u, 1048938190u,
		 1064473848u, 1050908101u, 3198391749u, 1064473848u,
		 1055089658u, 1063607819u, 3211091467u, 1055089658u,
		 1061212378u, 1059612165u, 3207095813u, 1061212378u,
		 1032462346u, 1065315015u, 3212798663u, 1032462346u,
		 1065299869u, 1034104972u, 3181588620u, 1065299869u,
		 1059456266u, 1061346857u, 3208830505u, 1059456266u,
		 1063515251u, 1055457463u, 3202941111u, 1063515251u,
		 1050517112u, 1064538431u, 3212022079u, 1050517112u,
		 1064717719u, 1049335047u, 3196818695u, 1064717719u,
		 1056546861u, 1063224103u, 3210707751u, 1056546861u,
		 1061738675u, 1058978834u, 3206462482u, 1061738675u,
		 1039022198u, 1065239331u, 3212722979u, 1039022198u,
		 1065088305u, 1043677336u, 3191160984u, 1065088305u,
		 1058152112u, 1062351814u, 3209835462u, 1058152112u,
		 1062694950u, 1057638573u, 3205122221u, 1062694950u,
		 1046100375u, 1064967812u, 3212451460u, 1046100375u,
		 1064191678u, 1052455328u, 3199938976u, 1064191678u,
		 1053596555u, 1063955323u, 3211438971u, 1053596555u,
		 1060655638u, 1060218909u, 3207702557u, 1060655638u,
		 1016514998u, 1065350374u, 3212834022u, 1016514998u,
		 1065350374u, 1016514998u, 3163998646u, 1065350374u,
		 1060218909u, 1060655638u, 3208139286u, 1060218909u,
		 1063955323u, 1053596555u, 3201080203u, 1063955323u,
		 1052455328u, 1064191678u, 3211675326u, 1052455328u,
		 1064967812u, 1046100375u, 3193584023u, 1064967812u,
		 1057638573u, 1062694950u, 3210178598u, 1057638573u,
		 1062351814u, 1058152112u, 3205635760u, 1062351814u,
		 1043677336u, 1065088305u, 3212571953u, 1043677336u,
		 1065239331u, 1039022198u, 3186505846u, 1065239331u,
		 1058978834u, 1061738675u, 3209222323u, 1058978834u,
		 1063224103u, 1056546861u, 3204030509u, 1063224103u,
		 1049335047u, 1064717719u, 3212201367u, 1049335047u,
		 1064538431u, 1050517112u, 3198000760u, 1064538431u,
		 1055457463u, 1063515251u, 3210998899u, 1055457463u,
		 1061346857u, 1059456266u, 3206939914u, 1061346857u,
		 1034104972u, 1065299869u, 3212783517u, 1034104972u,
		 1065315015u, 1032462346u, 3179945994u, 1065315015u,
		 1059612165u, 1061212378u, 3208696026u, 1059612165u,
		 1063607819u, 1055089658u, 3202573306u, 1063607819u,
		 1050908101u, 1064473848u, 3211957496u, 1050908101u,
		 1064772631u, 1048938190u, 3196421838u, 1064772631u,
		 1056905138u, 1063122622u, 3210606270u, 1056905138u,
		 1061865343u, 1058816524u, 3206300172u, 1061865343u,
		 1040422146u, 1065214129u, 3212697777u, 1040422146u,
		 1065123504u, 1042865916u, 3190349564u, 1065123504u,
		 1058320441u, 1062233263u, 3209716911u, 1058320441u,
		 1062805098u, 1057464630u, 3204948278u, 1062805098u,
		 1046903910u, 1064922702u, 3212406350u, 1046903910u,
		 1064265776u, 1052071148u, 3199554796u, 1064265776u,
		 1053973001u, 1063871889u, 3211355537u, 1053973001u,
		 1060797604u, 1060069797u, 3207553445u, 1060797604u,
		 1023101370u, 1065345321u, 3212828969u, 1023101370u,
		 1065337743u, 1026547719u, 3174031367u, 1065337743u,
		 1059918953u, 1060937731u, 3208421379u, 1059918953u,
		 1063786152u, 1054347371u, 3201831019u, 1063786152u,
		 1051685178u, 1064337510u, 3211821158u, 1051685178u,
		 1064875131u, 1047705169u, 3195188817u, 1064875131u,
		 1057289348u, 1062913104u, 3210396752u, 1057289348u,
		 1062112656u, 1058487303u, 3205970951u, 1062112656u,
		 1042052830u, 1065156211u, 3212639859u, 1042052830u,
		 1065186420u, 1041238199u, 3188721847u, 1065186420u,
		 1058652672u, 1061990009u, 3209473657u, 1058652672u,
		 1063018951u, 1057112753u, 3204596401u, 1063018951u,
		 1048504033u, 1064825104u, 3212308752u, 1048504033u,
		 1064406871u, 1051297476u, 3198781124u, 1064406871u,
		 1054719609u, 1063698124u, 3211181772u, 1054719609u,
		 1061075995u, 1059766402u, 3207250050u, 1061075995u,
		 1029837929u, 1065327640u, 3212811288u, 1029837929u,
		 1065282205u, 1035745987u, 3183229635u, 1065282205u,
		 1059298729u, 1061479413u, 3208963061u, 1059298729u,
		 1063420432u, 1055822969u, 3203306617u, 1063420432u,
		 1050124567u, 1064600610u, 3212084258u, 1050124567u,
		 1064660375u, 1049730525u, 3197214173u, 1064660375u,
		 1056186119u, 1063323378u, 3210807026u, 1056186119u,
		 1061610026u, 1059139577u, 3206623225u, 1061610026u,
		 1037385145u, 1065262025u, 3212745673u, 1037385145u,
		 1065050620u, 1044486967u, 3191970615u, 1065050620u,
		 1057982340u, 1062468291u, 3209951939u, 1057982340u,
		 1062582675u, 1057811152u, 3205294800u, 1062582675u,
		 1045294688u, 1065010454u, 3212494102u, 1045294688u,
		 1064115229u, 1052837662u, 3200321310u, 1064115229u,
		 1053218089u, 1064036440u, 3211520088u, 1053218089u,
		 1060511852u, 1060366268u, 3207849916u, 1060511852u,
		 1003032456u, 1065352900u, 3212836548u, 1003032456u,
		 1065353137u,  994643910u, 3142127558u, 1065353137u,
		 1060402831u, 1060475623u, 3207959271u, 1060402831u,
		 1064056356u, 1053123164u, 3200606812u, 1064056356u,
		 1052932949u, 1064095751u, 3211579399u, 1052932949u,
		 1065020727u, 1045092943u, 3192576591u, 1065020727u,
		 1057854081u, 1062554276u, 3210037924u, 1057854081u,
		 1062497083u, 1057939675u, 3205423323u, 1062497083u,
		 1044689081u, 1065040811u, 3212524459u, 1044689081u,
		 1065267305u, 1036975543u, 3184459191u, 1065267305u,
		 1059179515u, 1061577556u, 3209061204u, 1059179515u,
		 1063347850u, 1056095555u, 3203579203u, 1063347850u,
		 1049829173u, 1064645661u, 3212129309u, 1049829173u,
		 1064615778u, 1050026194u, 3197509842u, 1064615778u,
		 1055913979u, 1063396378u, 3210880026u, 1055913979u,
		 1061512249u, 1059259091u, 3206742739u, 1061512249u,
		 1036155961u, 1065277396u, 3212761044u, 1036155961u,
		 1065330403u, 1029015566u, 3176499214u, 1065330403u,
		 1059804699u, 1061041605u, 3208525253u, 1059804699u,
		 1063720345u, 1054626753u, 3202110401u, 1063720345u,
		 1051394561u, 1064389754u, 3211873402u, 1051394561u,
		 1064837841u, 1048304548u, 3195788196u, 1064837841u,
		 1057157023u, 1062992692u, 3210476340u, 1057157023u,
		 1062020861u, 1058611471u, 3206095119u, 1062020861u,
		 1041441994u, 1065179102u, 3212662750u, 1041441994u,
		 1065163997u, 1041849312u, 3189332960u, 1065163997u,
		 1058528786u, 1062082185u, 3209565833u, 1058528786u,
		 1062939769u, 1057245321u, 3204728969u, 1062939769u,
		 1047905114u, 1064862854u, 3212346502u, 1047905114u,
		 1064355073u, 1051588412u, 3199072060u, 1064355073u,
		 1054440633u, 1063764359u, 3211248007u, 1054440633u,
		 1060972472u, 1059880975u, 3207364623u, 1060972472u,
		 1027370453u, 1065335454u, 3212819102u, 1027370453u,
		 1065346821u, 1021454970u, 3168938618u, 1065346821u,
		 1060107238u, 1060762284u, 3208245932u, 1060107238u,
		 1063892964u, 1053879082u, 3201362730u, 1063892964u,
		 1052167363u, 1064247472u, 3211731120u, 1052167363u,
		 1064934211u, 1046703235u, 3194186883u, 1064934211u,
		 1057508242u, 1062777761u, 3210261409u, 1057508242u,
		 1062263095u, 1058278495u, 3205762143u, 1062263095u,
		 1043068932u, 1065114938u, 3212598586u, 1043068932u,
		 1065220664u, 1040217925u, 3187701573u, 1065220664u,
		 1058857247u, 1061833863u, 3209317511u, 1058857247u,
		 1063148198u, 1056815803u, 3204299451u, 1063148198u,
		 1049037530u, 1064759131u, 3212242779u, 1049037530u,
		 1064490219u, 1050810503u, 3198294151u, 1064490219u,
		 1055181822u, 1063584889u, 3211068537u, 1055181822u,
		 1061246177u, 1059573345u, 3207056993u, 1061246177u,
		 1032873140u, 1065311465u, 3212795113u, 1032873140u,
		 1065303892u, 1033694457u, 3181178105u, 1065303892u,
		 1059495395u, 1061313417u, 3208797065u, 1059495395u,
		 1063538604u, 1055365725u, 3202849373u, 1063538604u,
		 1050615007u, 1064522510u, 3212006158u, 1050615007u,
		 1064731675u, 1049235959u, 3196719607u, 1064731675u,
		 1056636663u, 1063198939u, 3210682587u, 1056636663u,
		 1061770529u, 1058938402u, 3206422050u, 1061770529u,
		 1039431104u, 1065233266u, 3212716914u, 1039431104u,
		 1065097338u, 1043474644u, 3190958292u, 1065097338u,
		 1058194330u, 1062322370u, 3209806018u, 1058194330u,
		 1062722687u, 1057595214u, 3205078862u, 1062722687u,
		 1046301466u, 1064956766u, 3212440414u, 1046301466u,
		 1064210424u, 1052359454u, 3199843102u, 1064210424u,
		 1053690858u, 1063934681u, 3211418329u, 1053690858u,
		 1060691301u, 1060181794u, 3207665442u, 1060691301u,
		 1018161769u, 1065349347u, 3212832995u, 1018161769u,
		 1065351242u, 1014714699u, 3162198347u, 1065351242u,
		 1060255914u, 1060619861u, 3208103509u, 1060255914u,
		 1063975820u, 1053502126u, 3200985774u, 1063975820u,
		 1052551087u, 1064172786u, 3211656434u, 1052551087u,
		 1064978704u, 1045899151u, 3193382799u, 1064978704u,
		 1057681847u, 1062667080u, 3210150728u, 1057681847u,
		 1062381129u, 1058109803u, 3205593451u, 1062381129u,
		 1043879916u, 1065079117u, 3212562765u, 1043879916u,
		 1065245240u, 1038613146u, 3186096794u, 1065245240u,
		 1059019167u, 1061706698u, 3209190346u, 1059019167u,
		 1063249129u, 1056456904u, 3203940552u, 1063249129u,
		 1049434048u, 1064703610u, 3212187258u, 1049434048u,
		 1064554201u, 1050419119u, 3197902767u, 1064554201u,
		 1055549057u, 1063491756u, 3210975404u, 1055549057u,
		 1061380177u, 1059417035u, 3206900683u, 1061380177u,
		 1034515386u, 1065295689u, 3212779337u, 1034515386u,
		 1065318408u, 1032051466u, 3179535114u, 1065318408u,
		 1059650881u, 1061178460u, 3208662108u, 1059650881u,
		 1063630608u, 1054997354u, 3202481002u, 1063630608u,
		 1051005599u, 1064457328u, 3211940976u, 1051005599u,
		 1064785978u, 1048838768u, 3196322416u, 1064785978u,
		 1056979462u, 1063096909u, 3210580557u, 1056979462u,
		 1061896698u, 1058775705u, 3206259353u, 1061896698u,
		 1040626286u, 1065207436u, 3212691084u, 1040626286u,
		 1065131914u, 1042662796u, 3190146444u, 1065131914u,
		 1058362295u, 1062203304u, 3209686952u, 1058362295u,
		 1062832301u, 1057420934u, 3204904582u, 1062832301u,
		 1047104442u, 1064911040u, 3212394688u, 1047104442u,
		 1064283931u, 1051974821u, 3199458469u, 1064283931u,
		 1054066791u, 1063850670u, 3211334318u, 1054066791u,
		 1060832809u, 1060032247u, 3207515895u, 1060832809u,
		 1024078895u, 1065343663u, 3212827311u, 1024078895u,
		 1065339874u, 1025724875u, 3173208523u, 1065339874u,
		 1059956825u, 1060902873u, 3208386521u, 1059956825u,
		 1063807801u, 1054253977u, 3201737625u, 1063807801u,
		 1051781835u, 1064319799u, 3211803447u, 1051781835u,
		 1064887254u, 1047505074u, 3194988722u, 1064887254u,
		 1057333292u, 1062886304u, 3210369952u, 1057333292u,
		 1062142999u, 1058445727u, 3205929375u, 1062142999u,
		 1042256251u, 1065148268u, 3212631916u, 1042256251u,
		 1065193582u, 1041034314u, 3188517962u, 1065193582u,
		 1058693779u, 1061959032u, 3209442680u, 1058693779u,
		 1063045073u, 1057068403u, 3204552051u, 1063045073u,
		 1048639680u, 1064812215u, 3212295863u, 1048639680u,
		 1064423839u, 1051200287u, 3198683935u, 1064423839u,
		 1054812329u, 1063675761u, 3211159409u, 1054812329u,
		 1061110268u, 1059728000u, 3207211648u, 1061110268u,
		 1030660152u, 1065324721u, 3212808369u, 1030660152u,
		 1065286857u, 1035335898u, 3182819546u, 1065286857u,
		 1059338266u, 1061446455u, 3208930103u, 1059338266u,
		 1063444347u, 1055731811u, 3203215459u, 1063444347u,
		 1050222846u, 1064585291u, 3212068939u, 1050222846u,
		 1064674939u, 1049631788u, 3197115436u, 1064674939u,
		 1056276533u, 1063298767u, 3210782415u, 1056276533u,
		 1061642373u, 1059099539u, 3206583187u, 1061642373u,
		 1037794615u, 1065256587u, 3212740235u, 1037794615u,
		 1065060274u, 1044284734u, 3191768382u, 1065060274u,
		 1058024917u, 1062439367u, 3209923015u, 1058024917u,
		 1062610942u, 1057768137u, 3205251785u, 1062610942u,
		 1045496305u, 1065000025u, 3212483673u, 1045496305u,
		 1064134561u, 1052742255u, 3200225903u, 1064134561u,
		 1053312892u, 1064016379u, 3211500027u, 1053312892u,
		 1060547967u, 1060329594u, 3207813242u, 1060547967u,
		 1008126808u, 1065352505u, 3212836153u, 1008126808u,
		 1065352505u, 1008126808u, 3155610456u, 1065352505u,
		 1060329594u, 1060547967u, 3208031615u, 1060329594u,
		 1064016379u, 1053312892u, 3200796540u, 1064016379u,
		 1052742255u, 1064134561u, 3211618209u, 1052742255u,
		 1065000025u, 1045496305u, 3192979953u, 1065000025u,
		 1057768137u, 1062610942u, 3210094590u, 1057768137u,
		 1062439367u, 1058024917u, 3205508565u, 1062439367u,
		 1044284734u, 1065060274u, 3212543922u, 1044284734u,
		 1065256587u, 1037794615u, 3185278263u, 1065256587u,
		 1059099539u, 1061642373u, 3209126021u, 1059099539u,
		 1063298767u, 1056276533u, 3203760181u, 1063298767u,
		 1049631788u, 1064674939u, 3212158587u, 1049631788u,
		 1064585291u, 1050222846u, 3197706494u, 1064585291u,
		 1055731811u, 1063444347u, 3210927995u, 1055731811u,
		 1061446455u, 1059338266u, 3206821914u, 1061446455u,
		 1035335898u, 1065286857u, 3212770505u, 1035335898u,
		 1065324721u, 1030660152u, 3178143800u, 1065324721u,
		 1059728000u, 1061110268u, 3208593916u, 1059728000u,
		 1063675761u, 1054812329u, 3202295977u, 1063675761u,
		 1051200287u, 1064423839u, 3211907487u, 1051200287u,
		 1064812215u, 1048639680u, 3196123328u, 1064812215u,
		 1057068403u, 1063045073u, 3210528721u, 1057068403u,
		 1061959032u, 1058693779u, 3206177427u, 1061959032u,
		 1041034314u, 1065193582u, 3212677230u, 1041034314u,
		 1065148268u, 1042256251u, 3189739899u, 1065148268u,
		 1058445727u, 1062142999u, 3209626647u, 1058445727u,
		 1062886304u, 1057333292u, 3204816940u, 1062886304u,
		 1047505074u, 1064887254u, 3212370902u, 1047505074u,
		 1064319799u, 1051781835u, 3199265483u, 1064319799u,
		 1054253977u, 1063807801u, 3211291449u, 1054253977u,
		 1060902873u, 1059956825u, 3207440473u, 1060902873u,
		 1025724875u, 1065339874u, 3212823522u, 1025724875u,
		 1065343663u, 1024078895u, 3171562543u, 1065343663u,
		 1060032247u, 1060832809u, 3208316457u, 1060032247u,
		 1063850670u, 1054066791u, 3201550439u, 1063850670u,
		 1051974821u, 1064283931u, 3211767579u, 1051974821u,
		 1064911040u, 1047104442u, 3194588090u, 1064911040u,
		 1057420934u, 1062832301u, 3210315949u, 1057420934u,
		 1062203304u, 1058362295u, 3205845943u, 1062203304u,
		 1042662796u, 1065131914u, 3212615562u, 1042662796u,
		 1065207436u, 1040626286u, 3188109934u, 1065207436u,
		 1058775705u, 1061896698u, 3209380346u, 1058775705u,
		 1063096909u, 1056979462u, 3204463110u, 1063096909u,
		 1048838768u, 1064785978u, 3212269626u, 1048838768u,
		 1064457328u, 1051005599u, 3198489247u, 1064457328u,
		 1054997354u, 1063630608u, 3211114256u, 1054997354u,
		 1061178460u, 1059650881u, 3207134529u, 1061178460u,
		 1032051466u, 1065318408u, 3212802056u, 1032051466u,
		 1065295689u, 1034515386u, 3181999034u, 1065295689u,
		 1059417035u, 1061380177u, 3208863825u, 1059417035u,
		 1063491756u, 1055549057u, 3203032705u, 1063491756u,
		 1050419119u, 1064554201u, 3212037849u, 1050419119u,
		 1064703610u, 1049434048u, 3196917696u, 1064703610u,
		 1056456904u, 1063249129u, 3210732777u, 1056456904u,
		 1061706698u, 1059019167u, 3206502815u, 1061706698u,
		 1038613146u, 1065245240u, 3212728888u, 1038613146u,
		 1065079117u, 1043879916u, 3191363564u, 1065079117u,
		 1058109803u, 1062381129u, 3209864777u, 1058109803u,
		 1062667080u, 1057681847u, 3205165495u, 1062667080u,
		 1045899151u, 1064978704u, 3212462352u, 1045899151u,
		 1064172786u, 1052551087u, 3200034735u, 1064172786u,
		 1053502126u, 1063975820u, 3211459468u, 1053502126u,
		 1060619861u, 1060255914u, 3207739562u, 1060619861u,
		 1014714699u, 1065351242u, 3212834890u, 1014714699u,
		 1065349347u, 1018161769u, 3165645417u, 1065349347u,
		 1060181794u, 1060691301u, 3208174949u, 1060181794u,
		 1063934681u, 1053690858u, 3201174506u, 1063934681u,
		 1052359454u, 1064210424u, 3211694072u, 1052359454u,
		 1064956766u, 1046301466u, 3193785114u, 1064956766u,
		 1057595214u, 1062722687u, 3210206335u, 1057595214u,
		 1062322370u, 1058194330u, 3205677978u, 1062322370u,
		 1043474644u, 1065097338u, 3212580986u, 1043474644u,
		 1065233266u, 1039431104u, 3186914752u, 1065233266u,
		 1058938402u, 1061770529u, 3209254177u, 1058938402u,
		 1063198939u, 1056636663u, 3204120311u, 1063198939u,
		 1049235959u, 1064731675u, 3212215323u, 1049235959u,
		 1064522510u, 1050615007u, 3198098655u, 1064522510u,
		 1055365725u, 1063538604u, 3211022252u, 1055365725u,
		 1061313417u, 1059495395u, 3206979043u, 1061313417u,
		 1033694457u, 1065303892u, 3212787540u, 1033694457u,
		 1065311465u, 1032873140u, 3180356788u, 1065311465u,
		 1059573345u, 1061246177u, 3208729825u, 1059573345u,
		 1063584889u, 1055181822u, 3202665470u, 1063584889u,
		 1050810503u, 1064490219u, 3211973867u, 1050810503u,
		 1064759131u, 1049037530u, 3196521178u, 1064759131u,
		 1056815803u, 1063148198u, 3210631846u, 1056815803u,
		 1061833863u, 1058857247u, 3206340895u, 1061833863u,
		 1040217925u, 1065220664u, 3212704312u, 1040217925u,
		 1065114938u, 1043068932u, 3190552580u, 1065114938u,
		 1058278495u, 1062263095u, 3209746743u, 1058278495u,
		 1062777761u, 1057508242u, 3204991890u, 1062777761u,
		 1046703235u, 1064934211u, 3212417859u, 1046703235u,
		 1064247472u, 1052167363u, 3199651011u, 1064247472u,
		 1053879082u, 1063892964u, 3211376612u, 1053879082u,
		 1060762284u, 1060107238u, 3207590886u, 1060762284u,
		 1021454970u, 1065346821u, 3212830469u, 1021454970u,
		 1065335454u, 1027370453u, 3174854101u, 1065335454u,
		 1059880975u, 1060972472u, 3208456120u, 1059880975u,
		 1063764359u, 1054440633u, 3201924281u, 1063764359u,
		 1051588412u, 1064355073u, 3211838721u, 1051588412u,
		 1064862854u, 1047905114u, 3195388762u, 1064862854u,
		 1057245321u, 1062939769u, 3210423417u, 1057245321u,
		 1062082185u, 1058528786u, 3206012434u, 1062082185u,
		 1041849312u, 1065163997u, 3212647645u, 1041849312u,
		 1065179102u, 1041441994u, 3188925642u, 1065179102u,
		 1058611471u, 1062020861u, 3209504509u, 1058611471u,
		 1062992692u, 1057157023u, 3204640671u, 1062992692u,
		 1048304548u, 1064837841u, 3212321489u, 1048304548u,
		 1064389754u, 1051394561u, 3198878209u, 1064389754u,
		 1054626753u, 1063720345u, 3211203993u, 1054626753u,
		 1061041605u, 1059804699u, 3207288347u, 1061041605u,
		 1029015566u, 1065330403u, 3212814051u, 1029015566u,
		 1065277396u, 1036155961u, 3183639609u, 1065277396u,
		 1059259091u, 1061512249u, 3208995897u, 1059259091u,
		 1063396378u, 1055913979u, 3203397627u, 1063396378u,
		 1050026194u, 1064615778u, 3212099426u, 1050026194u,
		 1064645661u, 1049829173u, 3197312821u, 1064645661u,
		 1056095555u, 1063347850u, 3210831498u, 1056095555u,
		 1061577556u, 1059179515u, 3206663163u, 1061577556u,
		 1036975543u, 1065267305u, 3212750953u, 1036975543u,
		 1065040811u, 1044689081u, 3192172729u, 1065040811u,
		 1057939675u, 1062497083u, 3209980731u, 1057939675u,
		 1062554276u, 1057854081u, 3205337729u, 1062554276u,
		 1045092943u, 1065020727u, 3212504375u, 1045092943u,
		 1064095751u, 1052932949u, 3200416597u, 1064095751u,
		 1053123164u, 1064056356u, 3211540004u, 1053123164u,
		 1060475623u, 1060402831u, 3207886479u, 1060475623u,
		  994643910u, 1065353137u, 3212836785u,  994643910u,
};

void
float_FFT(fpr32 *f, unsigned logn)
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
			fpr32 s_re, s_im;

			s_re = vrfy_fpr_gm_tab[((m + i1) << 1) + 0];
			s_im = vrfy_fpr_gm_tab[((m + i1) << 1) + 1];
			for (j = j1; j < j2; j ++) {
				fpr32 x_re, x_im, y_re, y_im, tmp;

				x_re = f[j];
				x_im = f[j + hn];
				y_re = f[j + ht];
				y_im = f[j + ht + hn];

				tmp  = fpr32_sub(fpr32_mul(y_re, s_re), fpr32_mul(y_im, s_im));
				y_im = fpr32_add(fpr32_mul(y_re, s_im), fpr32_mul(y_im, s_re));
				y_re = tmp;

				f[j      ] = fpr32_add(x_re, y_re);
				f[j   +hn] = fpr32_add(x_im, y_im);
				f[j+ht   ] = fpr32_sub(x_re, y_re);
				f[j+ht+hn] = fpr32_sub(x_im, y_im);
			}
		}
		t = ht;
	}
}

void
float_iFFT(fpr32 *f, unsigned logn)
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
			fpr32 s_re, s_im;

			s_re = vrfy_fpr_gm_tab[((hm + i1) << 1) + 0];
			s_im = fpr32_neg(vrfy_fpr_gm_tab[((hm + i1) << 1) + 1]);
			for (j = j1; j < j2; j ++) {
				fpr32 x_re, x_im, y_re, y_im, tmp;

				x_re = f[j];
				x_im = f[j + hn];
				y_re = f[j + t];
				y_im = f[j + t + hn];

				f[j     ] = fpr32_add(x_re, y_re);
				f[j + hn] = fpr32_add(x_im, y_im);

				x_re = fpr32_sub(x_re, y_re);
				x_im = fpr32_sub(x_im, y_im);

				tmp = fpr32_sub(
					fpr32_mul(x_re, s_re), fpr32_mul(x_im, s_im));
				f[j + t + hn] = fpr32_add(
					fpr32_mul(x_re, s_im), fpr32_mul(x_im, s_re));
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
		fpr32 ni = fpr32_scaled(1, 1 - logn);
		for (u = 0; u < n; u ++) {
			f[u] = fpr32_mul(f[u], ni);
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

	size_t n, hn, u, v, w;
	int16_t *s0;
	uint32_t p, p0i, R, R2, *gm, *igm, *q10_i32, *s1_i32;
	fpr32 *t0, *t1;
	uint8_t h0;

	n = MKN(logn);
	hn = n >> 1;

	t0 = (fpr32 *)tmp;
	t1 = t0 + n;

	gm = (uint32_t *)(tmp + 16*n);
	igm = gm + n;
	s0 = (int16_t *)(igm + n);
	q10_i32 = igm + n;
	s1_i32 = ((uint32_t *)tmp) + 3*n;

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
		t0[u] = fpr32_of(modp_norm(s1_i32[u], p));
		t1[u] = fpr32_of(q00[u]);
	}
	float_FFT(t0, logn);
	float_FFT(t1, logn);

	for (u = 0; u < hn; u++) {
		t0[u] = fpr32_div(t0[u], t1[u]);
		t0[u + hn] = fpr32_div(t0[u + hn], t1[u]);
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
			s0[u] = fpr32_rint(fpr32_half(
				fpr32_add(fpr32_of(h0 & 1), t0[u])));
			h0 >>= 1;
		}
	} else {
		for (u = 0, w = 0; w < n; u ++) {
			h0 = h[u];
			for (v = 0; v < 8; v ++, w ++) {
				s0[w] = fpr32_rint(fpr32_half(
					fpr32_add(fpr32_of(h0 & 1), t0[w])));
				h0 >>= 1;
			}
		}
	}

	return verify_norm_NTT(h, s0, s1, q00, q10, p, p0i, R, R2, logn, tmp);
}
