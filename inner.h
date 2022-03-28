#ifndef HAWK_INNER_H__
#define HAWK_INNER_H__

/*
 * Internal functions for Falcon. This is not the API intended to be
 * used by applications; instead, this internal API provides all the
 * primitives on which wrappers build to provide external APIs.
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

/*
 * IMPORTANT API RULES
 * -------------------
 *
 * This API has some non-trivial usage rules:
 *
 *
 *  - All public functions (i.e. the non-static ones) must be referenced
 *    with the Zf() macro (e.g. Zf(keygen) for the keygen()
 *    function). That macro adds a prefix to the name, which is
 *    configurable with the HAWK_PREFIX macro. This allows compiling
 *    the code into a specific "namespace" and potentially including
 *    several versions of this code into a single application (e.g. to
 *    have an AVX2 and a non-AVX2 variants and select the one to use at
 *    runtime based on availability of AVX2 opcodes).
 *
 *  - Functions that need temporary buffers expects them as a final
 *    tmp[] array of type uint8_t*, with a size which is documented for
 *    each function. However, most have some alignment requirements,
 *    because they will use the array to store 16-bit, 32-bit or 64-bit
 *    values (e.g. uint64_t or double). The caller must ensure proper
 *    alignment. What happens on unaligned access depends on the
 *    underlying architecture, ranging from a slight time penalty
 *    to immediate termination of the process.
 *
 *  - Some functions rely on specific rounding rules and precision for
 *    floating-point numbers. On some systems (in particular 32-bit x86
 *    with the 387 FPU), this requires setting an hardware control
 *    word. The caller MUST use set_fpu_cw() to ensure proper precision:
 *
 *      oldcw = set_fpu_cw(2);
 *      Zf(sign_dyn)(...);
 *      set_fpu_cw(oldcw);
 *
 *    On systems where the native floating-point precision is already
 *    proper, or integer-based emulation is used, the set_fpu_cw()
 *    function does nothing, so it can be called systematically.
 */


#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* =================================================================== */

/*
 * Compute degree N from logarithm 'logn'.
 */
#define MKN(logn)   ((size_t)1 << (logn))

/* ==================================================================== */

/*
 * "Naming" macro used to apply a consistent prefix over all global
 * symbols.
 */
#ifndef HAWK_PREFIX
#define HAWK_PREFIX   hawk_inner
#endif
#define Zf(name)             Zf_(HAWK_PREFIX, name)
#define Zf_(prefix, name)    Zf__(prefix, name)
#define Zf__(prefix, name)   prefix ## _ ## name


/*
 * Some computations with floating-point elements, in particular
 * rounding to the nearest integer, rely on operations using _exactly_
 * the precision of IEEE-754 binary64 type (i.e. 52 bits). On 32-bit
 * x86, the 387 FPU may be used (depending on the target OS) and, in
 * that case, may use more precision bits (i.e. 64 bits, for an 80-bit
 * total type length); to prevent miscomputations, we define an explicit
 * function that modifies the precision in the FPU control word.
 *
 * set_fpu_cw() sets the precision to the provided value, and returns
 * the previously set precision; callers are supposed to restore the
 * previous precision on exit. The correct (52-bit) precision is
 * configured with the value "2". On unsupported compilers, or on
 * targets other than 32-bit x86, or when the native 'double' type is
 * not used, the set_fpu_cw() function does nothing at all.
 */
static inline unsigned
set_fpu_cw(unsigned x)
{
	return x;
}



/*
 * MSVC 2015 does not know the C99 keyword 'restrict'.
 */
#if defined _MSC_VER && _MSC_VER
#ifndef restrict
#define restrict   __restrict
#endif
#endif

/* ==================================================================== */
/*
 * SHAKE256 implementation (shake.c).
 *
 * API is defined to be easily replaced with the fips202.h API defined
 * as part of PQClean.
 */

typedef struct {
	union {
		uint64_t A[25];
		uint8_t dbuf[200];
	} st;
	uint64_t dptr;
} inner_shake256_context;

#define inner_shake256_init      Zf(i_shake256_init)
#define inner_shake256_inject    Zf(i_shake256_inject)
#define inner_shake256_flip      Zf(i_shake256_flip)
#define inner_shake256_extract   Zf(i_shake256_extract)

void Zf(i_shake256_init)(inner_shake256_context *sc);
void Zf(i_shake256_inject)(inner_shake256_context *sc,
	const uint8_t *in, size_t len);
void Zf(i_shake256_flip)(inner_shake256_context *sc);
void Zf(i_shake256_extract)(inner_shake256_context *sc,
	uint8_t *out, size_t len);

/* ==================================================================== */
/*
 * RNG (rng.c).
 *
 * A PRNG based on ChaCha20 is implemented; it is seeded from a SHAKE256
 * context (flipped) and is used for bulk pseudorandom generation.
 * A system-dependent seed generator is also provided.
 */

/*
 * Obtain a random seed from the system RNG.
 *
 * Returned value is 1 on success, 0 on error.
 */
int Zf(get_seed)(void *seed, size_t seed_len);

/*
 * Structure for a PRNG. This includes a large buffer so that values
 * get generated in advance. The 'state' is used to keep the current
 * PRNG algorithm state (contents depend on the selected algorithm).
 *
 * The unions with 'dummy_u64' are there to ensure proper alignment for
 * 64-bit direct access.
 */
typedef struct {
	union {
		uint8_t d[512]; /* MUST be 512, exactly */
		uint64_t dummy_u64;
	} buf;
	size_t ptr;
	union {
		uint8_t d[256];
		uint64_t dummy_u64;
	} state;
	int type;
} prng;

/*
 * Instantiate a PRNG. That PRNG will feed over the provided SHAKE256
 * context (in "flipped" state) to obtain its initial state.
 */
void Zf(prng_init)(prng *p, inner_shake256_context *src);

/*
 * Refill the PRNG buffer. This is normally invoked automatically, and
 * is declared here only so that prng_get_u64() may be inlined.
 */
void Zf(prng_refill)(prng *p);

/*
 * Get some bytes from a PRNG.
 */
void Zf(prng_get_bytes)(prng *p, void *dst, size_t len);

/*
 * Get a 64-bit random value from a PRNG.
 */
static inline uint64_t
prng_get_u64(prng *p)
{
	size_t u;

	/*
	 * If there are less than 9 bytes in the buffer, we refill it.
	 * This means that we may drop the last few bytes, but this allows
	 * for faster extraction code. Also, it means that we never leave
	 * an empty buffer.
	 */
	u = p->ptr;
	if (u >= (sizeof p->buf.d) - 9) {
		Zf(prng_refill)(p);
		u = 0;
	}
	p->ptr = u + 8;

	/*
	 * On systems that use little-endian encoding and allow
	 * unaligned accesses, we can simply read the data where it is.
	 */
	return (uint64_t)p->buf.d[u + 0]
		| ((uint64_t)p->buf.d[u + 1] << 8)
		| ((uint64_t)p->buf.d[u + 2] << 16)
		| ((uint64_t)p->buf.d[u + 3] << 24)
		| ((uint64_t)p->buf.d[u + 4] << 32)
		| ((uint64_t)p->buf.d[u + 5] << 40)
		| ((uint64_t)p->buf.d[u + 6] << 48)
		| ((uint64_t)p->buf.d[u + 7] << 56);
}

/*
 * Get an 8-bit random value from a PRNG.
 */
static inline unsigned
prng_get_u8(prng *p)
{
	unsigned v;

	v = p->buf.d[p->ptr ++];
	if (p->ptr == sizeof p->buf.d) {
		Zf(prng_refill)(p);
	}
	return v;
}

/* ==================================================================== */
/*
 * Implementation of floating-point real numbers (fpr.h, fpr.c).
 */

/*
 * Real numbers are implemented by an extra header file, included below.
 * This is meant to support pluggable implementations. The default
 * implementation relies on the C type 'double'.
 *
 * The included file must define the following types, functions and
 * constants:
 *
 *   fpr
 *         type for a real number
 *
 *   fpr fpr_of(int64_t i)
 *         cast an integer into a real number; source must be in the
 *         -(2^63-1)..+(2^63-1) range
 *
 *   fpr fpr_scaled(int64_t i, int sc)
 *         compute i*2^sc as a real number; source 'i' must be in the
 *         -(2^63-1)..+(2^63-1) range
 *
 *   fpr fpr_ldexp(fpr x, int e)
 *         compute x*2^e
 *
 *   int64_t fpr_rint(fpr x)
 *         round x to the nearest integer; x must be in the -(2^63-1)
 *         to +(2^63-1) range
 *
 *   int64_t fpr_trunc(fpr x)
 *         round to an integer; this rounds towards zero; value must
 *         be in the -(2^63-1) to +(2^63-1) range
 *
 *   fpr fpr_add(fpr x, fpr y)
 *         compute x + y
 *
 *   fpr fpr_sub(fpr x, fpr y)
 *         compute x - y
 *
 *   fpr fpr_neg(fpr x)
 *         compute -x
 *
 *   fpr fpr_half(fpr x)
 *         compute x/2
 *
 *   fpr fpr_double(fpr x)
 *         compute x*2
 *
 *   fpr fpr_mul(fpr x, fpr y)
 *         compute x * y
 *
 *   fpr fpr_sqr(fpr x)
 *         compute x * x
 *
 *   fpr fpr_inv(fpr x)
 *         compute 1/x
 *
 *   fpr fpr_div(fpr x, fpr y)
 *         compute x/y
 *
 *   fpr fpr_sqrt(fpr x)
 *         compute the square root of x
 *
 *   int fpr_lt(fpr x, fpr y)
 *         return 1 if x < y, 0 otherwise
 *
 *   uint64_t fpr_expm_p63(fpr x)
 *         return exp(x), assuming that 0 <= x < log(2). Returned value
 *         is scaled to 63 bits (i.e. it really returns 2^63*exp(-x),
 *         rounded to the nearest integer). Computation should have a
 *         precision of at least 45 bits.
 *
 *   const fpr fpr_gm_tab[]
 *         array of constants for FFT / iFFT
 *
 *   const fpr fpr_p2_tab[]
 *         precomputed powers of 2 (by index, 0 to 10)
 *
 * Constants of type 'fpr':
 *
 *   fpr fpr_q                 12289
 *   fpr fpr_inverse_of_q      1/12289
 *   fpr fpr_inv_2sqrsigma0    1/(2*(1.8205^2))
 *   fpr fpr_inv_sigma[]       1/sigma (indexed by logn, 1 to 10)
 *   fpr fpr_sigma_min[]       1/sigma_min (indexed by logn, 1 to 10)
 *   fpr fpr_log2              log(2)
 *   fpr fpr_inv_log2          1/log(2)
 *   fpr fpr_bnorm_max         16822.4121
 *   fpr fpr_zero              0
 *   fpr fpr_one               1
 *   fpr fpr_two               2
 *   fpr fpr_onehalf           0.5
 *   fpr fpr_ptwo31            2^31
 *   fpr fpr_ptwo31m1          2^31-1
 *   fpr fpr_mtwo31m1          -(2^31-1)
 *   fpr fpr_ptwo63m1          2^63-1
 *   fpr fpr_mtwo63m1          -(2^63-1)
 *   fpr fpr_ptwo63            2^63
 */
#include "fpr.h"

/* ==================================================================== */
/*
 * Discrete gaussian sampling (sampler.c).
 */

/*
 * Internal sampler engine. Exported for tests.
 *
 * sampler_context wraps around a source of random numbers (PRNG) and
 * the sigma_min value (nominally dependent on the degree).
 *
 * sampler() takes as parameters:
 *   ctx      pointer to the sampler_context structure
 *   mu       center for the distribution
 *   isigma   inverse of the distribution standard deviation
 * It returns an integer sampled along the Gaussian distribution centered
 * on mu and of standard deviation sigma = 1/isigma.
 *
 * gaussian0_sampler() takes as parameter a pointer to a PRNG, and
 * returns an integer sampled along a half-Gaussian with standard
 * deviation sigma0 = 1.8205 (center is 0, returned value is
 * nonnegative).
 */

typedef struct {
	prng p;
	fpr sigma_min;
} sampler_context;

int Zf(sampler)(void *ctx, fpr mu, fpr isigma);

int Zf(gaussian0_sampler)(prng *p);

/* ==================================================================== */
/*
 * Coding & decoding (compress.c)
 */

/*
 * Encodes the NTRU-basis which is the private key of the signature scheme,
 * with an efficient encoding. The smallest private key is obtained by storing
 * the seed for the key generation, while the largest is storing the basis as
 * 16bit integers. Here, a middle-way is chosen.
 */
size_t
Zf(encode_seckey)(void *out, size_t max_out_len,
	const int8_t *f, const int8_t *g, const int8_t *F, unsigned logn);

/*
 * Inverse of Zf(encode_seckey). To calculate G, use Zf(expand_seckey).
 */
size_t
Zf(decode_seckey)(int8_t *f, int8_t *g, int8_t *F,
	const void *in, size_t max_in_len,  unsigned logn);

/*
 * Encode the Gram matrix of the NTRU-basis with determinant 1.
 */
size_t
Zf(encode_pubkey)(void *out, size_t max_out_len,
	const int16_t *q00, const int16_t *q10, unsigned logn);

/*
 * Decode the Gram matrix of the NTRU-basis with determinant 1.
 * The q11 can be recovered by q11 = (1 + q10 q10^*) / q00.
 */
size_t
Zf(decode_pubkey)(int16_t *q00, int16_t *q10,
	const void *in, size_t max_in_len, unsigned logn);

/*
 * Encode a signature s1 by outputting 'lo_bits' bits of the lowest signficant
 * bits of x[i] and using unary for the other most significant bits.
 */
size_t
Zf(encode_sig)(void *out, size_t max_out_len, const int16_t *x, unsigned logn,
	size_t lo_bits);

/*
 * Decode a signature s1.
 */
size_t
Zf(decode_sig)(int16_t *x, const void *in, size_t max_in_len, unsigned logn,
	size_t lo_bits);

/*
 * TODO: remove Huffman-table compression
 */
size_t
Zf(encode_pubkey_huffman)(void *out, size_t max_out_len,
	const int16_t *q00, const int16_t *q10, unsigned logn);

size_t
Zf(decode_pubkey_huffman)(int16_t *q00, int16_t *q10,
	const void *in, size_t max_in_len, unsigned logn);

size_t
Zf(encode_sig_huffman)(void *out, size_t max_out_len, const int16_t *x,
	unsigned logn);

/* ==================================================================== */
/*
 * FFT (fft.c).
 *
 * A real polynomial is represented as an array of N 'fpr' elements. The
 * FFT representation of a real polynomial contains N/2 complex elements;
 * each is stored as two real numbers, for the real and imaginary parts,
 * respectively. See fft.c for details on the internal representation.
 */

/*
 * Compute FFT in-place: the source array should contain a real
 * polynomial (N coefficients); its storage area is reused to store
 * the FFT representation of that polynomial (N/2 complex numbers).
 *
 * 'logn' MUST lie between 1 and 10 (inclusive).
 */
void Zf(FFT)(fpr *f, unsigned logn);

/*
 * Compute the inverse FFT in-place: the source array should contain the
 * FFT representation of a real polynomial (N/2 elements); the resulting
 * real polynomial (N coefficients of type 'fpr') is written over the
 * array.
 *
 * 'logn' MUST lie between 1 and 10 (inclusive).
 */
void Zf(iFFT)(fpr *f, unsigned logn);

/*
 * Add polynomial b to polynomial a. a and b MUST NOT overlap. This
 * function works in both normal and FFT representations.
 */
void Zf(poly_add)(fpr *restrict a, const fpr *restrict b, unsigned logn);

/*
 * Subtract polynomial b from polynomial a. a and b MUST NOT overlap. This
 * function works in both normal and FFT representations.
 */
void Zf(poly_sub)(fpr *restrict a, const fpr *restrict b, unsigned logn);

/*
 * Negate polynomial a. This function works in both normal and FFT
 * representations.
 */
void Zf(poly_neg)(fpr *a, unsigned logn);

/*
 * Compute adjoint of polynomial a. This function works only in FFT
 * representation.
 */
void Zf(poly_adj_fft)(fpr *a, unsigned logn);

/*
 * Multiply polynomial a with polynomial b. a and b MUST NOT overlap.
 * This function works only in FFT representation.
 */
void Zf(poly_mul_fft)(fpr *restrict a, const fpr *restrict b, unsigned logn);

/*
 * Store product of polynomial a and polynomial b in d. d, a and b MUST NOT
 * overlap. This function works only in FFT representation.
 */
void Zf(poly_prod_fft)(fpr *restrict d, const fpr *restrict a,
	const fpr *restrict b, unsigned logn);

/*
 * Multiply polynomial a with the adjoint of polynomial b. a and b MUST NOT
 * overlap. This function works only in FFT representation.
 */
void Zf(poly_muladj_fft)(fpr *restrict a, const fpr *restrict b, unsigned logn);

/*
 * Multiply polynomial with its own adjoint. This function works only in FFT
 * representation.
 */
void Zf(poly_mulselfadj_fft)(fpr *a, unsigned logn);

/*
 * Multiply polynomial with a real constant. This function works in both
 * normal and FFT representations.
 */
void Zf(poly_mulconst)(fpr *a, fpr x, unsigned logn);

/*
 * Divide polynomial a by polynomial b, modulo X^N+1 (FFT representation).
 * a and b MUST NOT overlap.
 */
void Zf(poly_div_fft)(fpr *restrict a, const fpr *restrict b, unsigned logn);

/*
 * Given f and g (in FFT representation), compute 1/(f*adj(f)+g*adj(g))
 * (also in FFT representation). Since the result is auto-adjoint, all its
 * coordinates in FFT representation are real; as such, only the first N/2
 * values of d[] are filled (the imaginary parts are skipped).
 *
 * Array d MUST NOT overlap with either a or b.
 */
void Zf(poly_invnorm2_fft)(fpr *restrict d,
	const fpr *restrict a, const fpr *restrict b, unsigned logn);

/*
 * Given F, G, f and g (in FFT representation), compute F*adj(f)+G*adj(g)
 * (also in FFT representation). Destination d MUST NOT overlap with
 * any of the source arrays.
 */
void Zf(poly_add_muladj_fft)(fpr *restrict d,
	const fpr *restrict F, const fpr *restrict G,
	const fpr *restrict f, const fpr *restrict g, unsigned logn);

/*
 * Given F, G, f and g (in FFT representation), compute F*f+G*g
 * (also in FFT representation). Destination d MUST NOT overlap with
 * any of the source arrays.
 */
void Zf(poly_add_mul_fft)(fpr *restrict d,
	const fpr *restrict F, const fpr *restrict G,
	const fpr *restrict f, const fpr *restrict g, unsigned logn);

/*
 * Multiply polynomial a by polynomial b, where b is autoadjoint. Both
 * a and b are in FFT representation. Since b is autoadjoint, all its
 * FFT coefficients are real, and the array b contains only N/2 elements.
 * a and b MUST NOT overlap.
 */
void Zf(poly_mul_autoadj_fft)(fpr *restrict a,
	const fpr *restrict b, unsigned logn);

/*
 * Divide polynomial a by polynomial b, where b is autoadjoint. Both
 * a and b are in FFT representation. Since b is autoadjoint, all its
 * FFT coefficients are real, and the array b contains only N/2 elements.
 * a and b MUST NOT overlap.
 */
void Zf(poly_div_autoadj_fft)(fpr *restrict a,
	const fpr *restrict b, unsigned logn);

/*
 * Perform an LDL decomposition of an auto-adjoint matrix G, in FFT
 * representation. On input, g00, g01 and g11 are provided (where the
 * matrix G = [[g00, g01], [adj(g01), g11]]). On output, the d00, l10
 * and d11 values are written in g00, g01 and g11, respectively
 * (with D = [[d00, 0], [0, d11]] and L = [[1, 0], [l10, 1]]).
 * (In fact, d00 = g00, so the g00 operand is left unmodified.)
 */
void Zf(poly_LDL_fft)(const fpr *restrict g00,
	fpr *restrict g01, fpr *restrict g11, unsigned logn);

/*
 * Perform an LDL decomposition of an auto-adjoint matrix G, in FFT
 * representation. This is identical to poly_LDL_fft() except that
 * g00, g01 and g11 are unmodified; the outputs d11 and l10 are written
 * in two other separate buffers provided as extra parameters.
 */
void Zf(poly_LDLmv_fft)(fpr *restrict d11, fpr *restrict l10,
	const fpr *restrict g00, const fpr *restrict g01,
	const fpr *restrict g11, unsigned logn);

/*
 * Apply "split" operation on a polynomial in FFT representation:
 * f = f0(x^2) + x*f1(x^2), for half-size polynomials f0 and f1
 * (polynomials modulo X^(N/2)+1). f0, f1 and f MUST NOT overlap.
 */
void Zf(poly_split_fft)(fpr *restrict f0, fpr *restrict f1,
	const fpr *restrict f, unsigned logn);

/*
 * Apply "merge" operation on two polynomials in FFT representation:
 * given f0 and f1, polynomials moduo X^(N/2)+1, this function computes
 * f = f0(x^2) + x*f1(x^2), in FFT representation modulo X^N+1.
 * f MUST NOT overlap with either f0 or f1.
 */
void Zf(poly_merge_fft)(fpr *restrict f,
	const fpr *restrict f0, const fpr *restrict f1, unsigned logn);


/*
 * Add to a polynomial its own adjoint. This function works only in FFT
 * representation.
 */
void
Zf(poly_addselfadj_fft)(fpr *a, unsigned logn);

/*
 * Add polynomial b to polynomial a, where b is autoadjoint. Both a and b
 * are in FFT representation. Since b is autoadjoint, all its FFT
 * coefficients are real, and the array b contains only N/2 elements.
 */
void
Zf(poly_add_autoadj_fft)(fpr *a, fpr *b, unsigned logn);

/* ==================================================================== */
/*
 * Fast Fourier Orthogonalization.
 */

/*
 * Number of elements (fpr) in the LDL tree for an input with polynomials
 * of size 2^logn. Note that leaves do not store any information, nor does
 * the root of the whole tree. As on any level i=1,..., logn-1 we have
 * 2^logn elements, the total size is (logn-1) * 2^logn.
 */
#define LDL_TREESIZE(logn) (((logn) - 1) << (logn))

/*
 * Compute the ffLDL tree of an auto-adjoint q00 (in FFT representation),
 * so one can use Zf(ffNearestPlane_tree)
 *
 * The "tree" array is filled with the computed tree, of size logn 2^logn
 * elements (see ffLDL_treesize()).
 *
 * Input arrays MUST NOT overlap. tmp[] should have room for at least
 * three polynomials of 2^logn elements each.
 */
void
Zf(ffLDL_fft)(fpr *restrict tree, const fpr *restrict q00, unsigned logn,
	fpr *restrict tmp);

/*
 * Perform Fast Fourier Nearest Plane for target vector t and LDL tree T.
 * Result is stored in z. Note: tmp[] must have space for at least two
 * polynomials of size 2^logn.
 */
void
Zf(ffNearestPlane_tree)(fpr *restrict z, const fpr *restrict tree,
	const fpr *restrict t, unsigned logn, fpr *restrict tmp);

/*
 * Dynamic version of Zf(ffNearestPlane_tree). The tree of Gram matrix g
 * is calculated on the fly, and the return value is stored in t. Note: t
 * and g get modified in this function. tmp[] must have space for at least
 * 2 polynomials of size 2^logn.
 */
void
Zf(ffNearestPlane_dyn)(fpr *restrict t, fpr *restrict g, unsigned logn,
	fpr *restrict tmp);

/*
 * Babai reduce (F, G) w.r.t (f, g) (all in FFT representation) and adjust
 * the result in Fn, Gn as well (in coefficient representation). Note:
 * tmp[] must have space for at least 4 polynomials of size 2^logn.
 */
void
Zf(ffBabai_reduce)(const fpr *restrict f, const fpr *restrict g,
	fpr *restrict F, fpr *restrict G, int8_t *restrict Fn,
	int8_t *restrict Gn, unsigned logn, fpr *tmp);

/* ==================================================================== */
/*
 * Key pair generation.
 */

/*
 * Completes the NTRU basis, by finding (F, G) in Babai's fundamental domain
 * generated by (rotations of) (f, g) such that f * G - g * F = 1. Moreover, it
 * computes the public key which is the Gram matrix of the basis
 *
 *     [ (f, g), (F, G) ].
 *
 * This function returns 1 on success and 0 on failure which may happen during
 * NTRU solving.
 *
 * Note: tmp[] must have space for at least 48 * 2^logn bytes.
 * If Babai reduction is not done, at least 28 2^logn bytes are needed.
 */
int
Zf(complete_private)(const int8_t *restrict f, const int8_t *restrict g,
	int8_t *restrict F, int8_t *restrict G,
	fpr *restrict q00, fpr *restrict q10, fpr *restrict q11,
	unsigned logn, uint8_t *restrict tmp);

/*
 * Compute a secret key and a public key belonging to the signature
 * scheme. The secret key is a tuple (f, g, F, G) where each belongs to
 * the ring Z[X] / (X^n + 1) with small coefficients (abs. value << 127),
 * such that
 *
 *     f * G - g * F = 1 (mod X^n+1),
 *
 * holds. The public is given by the Gram matrix of the basis generated by
 * the basis vectors (f, g) and (F, G). A discrete gaussian with sigma = 1.425
 * is used to generate coefficients of f and g.
 *
 * Note: tmp[] must have space for at least 48 * 2^logn bytes.
 * If Babai reduction is not done, at least 28 2^logn bytes are needed.
 */
void
Zf(keygen)(inner_shake256_context *rng,
	int8_t *restrict f, int8_t *restrict g, // secret key
	int8_t *restrict F, int8_t *restrict G, // secret key
	fpr *restrict q00, fpr *restrict q10, fpr *restrict q11, // public key
	unsigned logn, uint8_t *restrict tmp);

/* ==================================================================== */
/*
 * Signature generation.
 */

/*
 * Compute a signature of hm: the signature is a vector (s0, s1) that is close
 * to (hm/2, 0) with respect to the quadratic form Q.
 * Here 'close' means that the squared l2-norm of (2 s0 - hm, 2 s1) wrt Q
 * is at most `bound`.
 *
 * If during generation it is not short enough, s0 and s1 are untouched and 0
 * is returned; the caller should then try again.
 * Otherwise, 1 is returned and s0, s1 contain a valid signature for hm.
 *
 * Note: tmp[] must have space for at least 40 * 2^logn bytes.
 */
void
Zf(complete_sign)(inner_shake256_context *rng,
	int16_t *restrict s0, int16_t *restrict s1,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const int8_t *restrict h0, const int8_t *restrict h1, uint32_t bound,
	unsigned logn, uint8_t *restrict tmp);

/*
 * Similar to Zf(complete_sign), except that s0 is not returned.
 * If Zf(sign) returns 1, a signature is only valid with high probability.
 *
 * Note: tmp[] must have space for at least 40 * 2^logn bytes.
 */
void
Zf(sign)(inner_shake256_context *rng, int16_t *restrict sig,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const int8_t *restrict h0, const int8_t *restrict h1, uint32_t bound,
	unsigned logn, uint8_t *restrict tmp);

/*
 * Similar to Zf(sign), except that the signature is always valid when 1 is
 * returned.
 *
 * Note: tmp[] must have space for at least 24 * 2^logn bytes.
 */
void
Zf(guaranteed_sign)(inner_shake256_context *rng, int16_t *restrict sig,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const fpr *restrict q00, const int8_t *restrict h0,
	const int8_t *restrict h1, uint32_t bound, unsigned logn,
	uint8_t *restrict tmp);

/*
 * Required space in terms of number of fpr values for an expanded secret key,
 * which can be used in Zf(fft_sign).
 */
#define EXPANDED_SECKEY_SIZE(logn) (9 << ((logn) - 1))

/*
 * Expands a secret key given by the key generation, to produce the secret key
 * basis [[f,g], [F,G]] in FFT-representation, as well as 1/(f^* f + g^* g).
 *
 * Note: expanded_seckey[] must have space for at least 4.5 2^logn fprs, since
 * there are 4 polynomials in the basis and 1/(f* f + g* g) is self-adjoint.
 */
void
Zf(expand_seckey)(fpr *restrict expanded_seckey,
	const int8_t *f, const int8_t *g, const int8_t *F, unsigned logn);

/*
 * Generates a signature of a message with hash h[] (of n / 8 bytes long),
 * which is guaranteed to be a valid signature.
 *
 * Note: tmp[] must have space for at least 26 * 2^logn bytes.
 */
void
Zf(fft_sign)(inner_shake256_context *rng, int16_t *restrict sig,
	const fpr *restrict expanded_seckey, const uint8_t *restrict h,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp);

/* ==================================================================== */
/*
 * Signature verification functions (vrfy.c).
 */

/*
 * Given q00_num, q10_num in integer representation, compute the full public
 * key, which is q00, q10 and q11 in FFT representation.
 * Note here that q11 is reconstructed using the rule q00 q11 - q10 q01 = 1.
 */
void
Zf(complete_pubkey)(const int16_t *q00_num, const int16_t *q10_num,
	fpr *q00, fpr *q10, fpr *q11, unsigned logn);

/*
 * Verify if a signature (s0, s1) is valid for a message hm.
 * The signature is accepted iff the squared l2-norm of (2s0 - hm, 2s1) is at
 * most bound.
 *
 * Note: tmp[] must have space for at least 16 * 2^logn bytes.
 */
int
Zf(complete_verify)(const int8_t *restrict h0, const int8_t *restrict h1,
	const int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp);

/*
 * Similar to Zf(complete_verify), except that it takes s1 as signature,
 * and reconstructs s0 with simple rounding:
 *     s0 = round((hm%2) / 2 - s1 q10 / q00).
 */
int
Zf(verify_simple_rounding)(const int8_t *restrict h0, const int8_t *restrict h1,
	int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp);

/*
 * Similar to Zf(complete_verify), except that it takes s1 as signature,
 * and reconstructs s0 with Babai's Nearest Plane Algorithm:
 *     s0 ~ (hm%2) / 2 - s1 q10 / q00,
 * such that s0 (f, g) + s1 (F, G) is in the parallelepiped generated by
 * the GSO'ed basis generated by (f, g) and (F, G).
 */
int
Zf(verify_nearest_plane)(const int8_t *restrict h0, const int8_t *restrict h1,
	int16_t *restrict s0, const int16_t *restrict s1,
	const fpr *restrict q00, const fpr *restrict q10, const fpr *restrict q11,
	uint32_t bound, unsigned logn, uint8_t *restrict tmp);

/*
 * Verify if a signature (s0, s1) is valid for a hashed message h of length
 * n / 4 bytes, where s0 is reconstructed with simple rounding.
 * The signature is accepted iff the squared l2-norm of (h0 - 2s0, h1 - 2s1) is
 * at most bound wrt quadratic form Q.
 *
 * Note: tmp[] must have space for at least 16 * 2^logn bytes.
 */
int
Zf(verify_simple_rounding_fft)(const uint8_t *restrict h,
	const int16_t *restrict s1, const fpr *restrict q00,
	const fpr *restrict q10, const fpr *restrict q11, uint32_t bound,
	unsigned logn, uint8_t *restrict tmp);

#endif
