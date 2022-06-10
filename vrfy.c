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
	size_t u, n, hn;

	n = MKN(logn);
	hn = n >> 1;

	// doing this in reverse, allows iq00, iq10 to overlap with the begin
	// of q00.
	for (u = n; u -- > 0; ) {
		q10[u] = fpr_of(iq10[u]);
	}
	for (u = n; u -- > 0; ) {
		q00[u] = fpr_of(iq00[u]);
	}

	Zf(FFT)(q00, logn);
	Zf(FFT)(q10, logn);

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
Zf(verify_simple_rounding)(const uint8_t *restrict h,
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
Zf(verify_simple_rounding_fft)(const uint8_t *restrict h,
	const int16_t *restrict s1, const fpr *restrict q00,
	const fpr *restrict q10, const fpr *restrict q11,
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
	{ 2147377153, 1977035326,  968223422 },
	{ 2147358721, 1067163706,  132460015 },
	{ 2147352577, 1606082042,  598693809 },
	{ 2147346433, 2033915641, 1056257184 },
	{ 2147338241, 1653770625,  421286710 },
	{ 2147309569,  631200819, 1111201074 },
	{ 2147297281, 2038364663, 1042003613 },
	{ 0, 0, 0 }
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

	n = (size_t)1 << logn;

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

/*
 * Set the polynomial a equal to h - 2 * s and converts it to Montgomery NTT
 * representation.
 */
static void
hash_to_i32(uint32_t *a, const uint8_t *h, const int16_t *s,
	uint32_t *gm, uint32_t R2, uint32_t p, uint32_t p0i, unsigned logn)
{
	int32_t *b;
	size_t n, u, v;
	uint8_t hash;

	n = MKN(logn);
	b = (int32_t *)a;

	if (logn <= 3) {
		for (u = 0; u < n; u ++) {
			a[u] = ((h[0] >> u) & 1) - 2 * s[u];
		}
	} else {
		for (u = 0; u < n; ) {
			hash = *h++;
			for (v = 0; v < 8; v ++, u ++) {
				a[u] = (hash & 1) - 2 * s[u];
				hash >>= 1;
			}
		}
	}

	for (u = 0; u < n; u++) {
		a[u] = modp_set(b[u], p);
	}
	modp_NTT2(a, gm, logn, p, p0i);
	for (u = 0; u < n; u++) {
		a[u] = modp_montymul(a[u], R2, p, p0i);
	}
}

/* see inner.h */
int
Zf(uncompressed_verify_NTT)(const uint8_t *restrict h,
	const int16_t *restrict s0, const int16_t *restrict s1,
	const int16_t *restrict q00, const int16_t *restrict q10,
	unsigned logn, uint8_t *restrict tmp)
{
	uint32_t norm, norm0, term, p, p0i, R, R2;
	uint32_t *gm, *igm, *q00_i32, *q10_i32, *q11_i32, *s0_i32, *s1_i32;
	int32_t *q11;
	size_t n, hn, u, v;

	n = MKN(logn);
	norm = 0;
	hn = n >> 1;
	q11 = (int32_t *)tmp;
	gm = (uint32_t *)(q11 + n);
	igm = gm + n;
	q00_i32 = igm + n;
	q10_i32 = q00_i32 + n;
	q11_i32 = q10_i32 + n; 
	s0_i32 = q11_i32 + n;
	s1_i32 = s0_i32 + n;

	p = PRIMES[0].p;
	p0i = modp_ninv31(p);
	R = modp_R(p);
	R2 = modp_R2(p, p0i);
	modp_mkgm2(gm, igm, logn, PRIMES[0].g, p, p0i);

	for (u = 0; u < n; u++) {
		q00_i32[u] = modp_set(q00[u], p);
		q10_i32[u] = modp_set(q10[u], p);
	}

	modp_NTT2(q00_i32, gm, logn, p, p0i);
	modp_NTT2(q10_i32, gm, logn, p, p0i);

	for (u = 0; u < hn; u++) {
		q11_i32[u] = modp_montymul(q10_i32[u], R2, p, p0i);
		q11_i32[u] = modp_montymul(q11_i32[u], q10_i32[n - 1 - u], p, p0i);
		q11_i32[u] = modp_add(q11_i32[u], 1u, p);
		q11_i32[u] = modp_div(q11_i32[u], q00_i32[u], p, p0i, R);

		// Needed for the iNTT:
		q11_i32[n - 1 - u] = q11_i32[u];
	}

	hash_to_i32(s0_i32, h, s0, gm, R2, p, p0i, logn);
	hash_to_i32(s1_i32, SECOND_HASH(h, logn), s1, gm, R2, p, p0i, logn);

	// s0* q00 s0
	for (u = 0; u < hn; u++) {
		term = modp_montymul(s0_i32[u], s0_i32[n - 1 - u], p, p0i);
		norm = modp_add(norm, modp_montymul(term, q00_i32[u], p, p0i), p);
	}

	// s1* q11 s1
	for (u = 0; u < hn; u++) {
		term = modp_montymul(s1_i32[u], s1_i32[n - 1 - u], p, p0i);
		norm = modp_add(norm, modp_montymul(term, q11_i32[u], p, p0i), p);
	}

	// s0* q01 s1 + s1* q10 s0
	for (u = 0; u < n; u++) {
		term = modp_montymul(s1_i32[u], s0_i32[n - 1 - u], p, p0i);
		norm = modp_add(norm, modp_montymul(term, q10_i32[u], p, p0i), p);
	}

	modp_iNTT2(q11_i32, igm, logn, p, p0i);
	for (u = 0; u < n; u++)
		q11[u] = modp_norm(q11_i32[u], p);

	norm0 = norm;
	for (v = 1; v < 3; v++) {
		norm = 0;
		p = PRIMES[v].p;
		p0i = modp_ninv31(p);
		R = modp_R(p);
		R2 = modp_R2(p, p0i);
		modp_mkgm2(gm, igm, logn, PRIMES[v].g, p, p0i);

		// TODO: we can optimize RAM by putting the q00, q10, q11 resp. in NTT
		// version when needed below
		for (u = 0; u < n; u++) {
			q00_i32[u] = modp_set(q00[u], p);
			q10_i32[u] = modp_set(q10[u], p);
			q11_i32[u] = modp_set(q11[u], p);
		}

		modp_NTT2(q00_i32, gm, logn, p, p0i);
		modp_NTT2(q10_i32, gm, logn, p, p0i);
		modp_NTT2(q11_i32, gm, logn, p, p0i);

		hash_to_i32(s0_i32, h, s0, gm, R2, p, p0i, logn);
		hash_to_i32(s1_i32, SECOND_HASH(h, logn), s1, gm, R2, p, p0i, logn);

		// s0* q00 s0
		for (u = 0; u < hn; u++) {
			term = modp_montymul(s0_i32[u], s0_i32[n - 1 - u], p, p0i);
			norm = modp_add(norm, modp_montymul(term, q00_i32[u], p, p0i), p);
		}

		// s1* q11 s1
		for (u = 0; u < hn; u++) {
			term = modp_montymul(s1_i32[u], s1_i32[n - 1 - u], p, p0i);
			norm = modp_add(norm, modp_montymul(term, q11_i32[u], p, p0i), p);
		}

		// s0* q01 s1 + s1* q10 s0
		for (u = 0; u < n; u++) {
			term = modp_montymul(s1_i32[u], s0_i32[n - 1 - u], p, p0i);
			norm = modp_add(norm, modp_montymul(term, q10_i32[u], p, p0i), p);
		}

		// norm0 := \infty, when norm0 != norm.
		norm0 |= -((norm0 - norm) >> 31);
		norm0 |= -((norm - norm0) >> 31);
	}

	return Zf(in_positive_half)(s1, SECOND_HASH(h, logn), logn)
		&& (norm0 & (hn - 1)) == 0
		&& (norm0 >> (logn - 1)) <= L2BOUND(logn);
}

#define hashbit(h, idx) ((h[idx >> 3] >> (idx & 7)) & 1u)

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
	// uint32_t p, p0i, R2;
	// uint32_t *gm, *igm, *q00_i32, *q10_i32, *e0_i32, *e1_i32;
	// fpr *numerator, *denominator;
	fpr *t0, *t1;
	uint8_t h0;

	n = MKN(logn);
	s0 = (int16_t *)tmp;

	t0 = (fpr *)tmp;
	t1 = t0 + n;

	hash_to_fft(t0, SECOND_HASH(h, logn), s1, logn);

	for (u = 0; u < n; u++) t1[u] = fpr_of(q10[u]);
	Zf(FFT)(t1, logn);
	Zf(poly_mul_fft)(t0, t1, logn);

	for (u = 0; u < n; u++) t1[u] = fpr_of(q00[u]);
	Zf(FFT)(t1, logn);
	Zf(poly_div_autoadj_fft)(t0, t1, logn);

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

/*
	gm = (uint32_t *)(s0 + n);
	igm = gm + n;
	q00_i32 = igm + n;
	q10_i32 = q00_i32 + n;
	e0_i32 = q10_i32 + n;
	e1_i32 = e0_i32 + n;
	numerator = (fpr *)(e1_i32 + n);
	denominator = (fpr *)(numerator + n);

	p = PRIMES[0].p;
	p0i = modp_ninv31(p);
	R2 = modp_R2(p, p0i);
	modp_mkgm2(gm, igm, logn, PRIMES[0].g, p, p0i);

	for (u = 0; u < n; u++) {
		q00_i32[u] = modp_set(q00[u], p);
		q10_i32[u] = modp_set(q10[u], p);
		e0_i32[u] = modp_set(hashbit(h, u), p);
		e1_i32[u] = modp_set(hashbit(SECOND_HASH(h, logn), u) - 2 * s1[u], p);
	}
	modp_NTT2(q00_i32, gm, logn, p, p0i);
	modp_NTT2(q10_i32, gm, logn, p, p0i);
	modp_NTT2(e0_i32, gm, logn, p, p0i);
	modp_NTT2(e1_i32, gm, logn, p, p0i);

	for (u = 0; u < n; u++) {
		e0_i32[u] = modp_montymul(q00_i32[u],
			modp_montymul(e0_i32[u], R2, p, p0i), p, p0i);
		e1_i32[u] = modp_montymul(q10_i32[u],
			modp_montymul(e1_i32[u], R2, p, p0i), p, p0i);

		e0_i32[u] = modp_sub(e0_i32[u], e1_i32[u], p);
	}

	modp_iNTT2(e0_i32, igm, logn, p, p0i);
	printf("\n");
	for (u = 0; u < n; u++) {
		printf("%d ", modp_norm(e0_i32[u], p));
	}
	printf("\n");
	for (u = 0; u < n; u ++) {
		numerator[u] = fpr_half(fpr_of(modp_norm(e0_i32[u], p)));
		denominator[u] = fpr_of(q00[u]);
	}
	Zf(FFT)(numerator, logn);
	Zf(FFT)(denominator, logn);

	Zf(poly_div_autoadj_fft)(numerator, denominator, logn);
	Zf(iFFT)(numerator, logn);

	for (u = 0; u < n; u ++) {
		s0[u] = fpr_rint(numerator[u]);
		printf("%d ", s0[u]);
	}
	printf("\n");
*/

	// return Zf(uncompressed_verify_NTT)(h, s0, s1, q00, q10, logn, (uint8_t *) gm);
	// return 1;
	return Zf(uncompressed_verify_NTT)(h, s0, s1, q00, q10, logn, (uint8_t *)(s0 + n));
}
