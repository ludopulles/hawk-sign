/*
 * External Hawk API.
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2022 Hawk Project
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
 * @author   Ludo Pulles <ludo.pulles@cwi.nl>
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

#ifndef HAWK_H__
#define HAWK_H__

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================= */
/*
 * Hawk API Notes
 * ----------------
 *
 *
 * HAWK DEGREE
 *
 * Hawk is parameterized by a degree, which is a power of two. Formally, two
 * values are possible: 512 for Hawk-512 and 1024 for Hawk-1024. This
 * implementation also supports lower degrees, from 2 to 256; these reduced
 * variants do NOT provide adequate security and should be used for research
 * purposes only.
 *
 * In all functions and macros defined below, the degree is provided
 * logarithmically as the 'logn' parameter: logn ranges from 1 to 10,
 * and represents the degree 2^logn.
 *
 *
 * ERROR REPORTING
 *
 * All functions that may fail for some reason return an 'int' value. A
 * returned value of zero is a success; all error conditions are
 * reported as an error code. Error codes are negative. Macros are
 * defined for some error codes; in the interest of forward
 * compatiblity, applications that use this implementation should be
 * prepared to receive other error codes not listed in the macros below.
 *
 *
 * TEMPORARY BUFFERS
 *
 * Many functions expect temporary areas, provided as the parameter
 * 'tmp'. The caller is responsible for allocating these areas with the
 * proper size; the HAWK_TMPSIZE_* macros evaluate to constant
 * expressions that yield the proper size (in bytes) and can be used to
 * allocate the temporaries on the stack or elsewhere. There are no
 * alignment requirements on temporaries (the functions handle alignment
 * internally).
 *
 * The caller is responsible for clearing temporary buffer contents,
 * if such memory scrubbing is deemed relevant in the context where this
 * implementation is used.
 *
 * The same temporary buffer can be reused for several operations,
 * possibly distinct from each other. For all degrees (logn = 1 to 10), the
 * following sizes are in descending order:
 *
 *   HAWK_TMPSIZE_SIGNDYN
 *   HAWK_TMPSIZE_KEYGEN
 *   HAWK_TMPSIZE_VERIFY
 *   HAWK_TMPSIZE_VERIFY_NTT
 *   HAWK_TMPSIZE_SIGNDYN_NTT
 *   HAWK_TMPSIZE_SIGN
 *   HAWK_TMPSIZE_MAKEPUB
 *
 * i.e. a temporary buffer large enough for dynamic floating-point signing
 * ("SIGNDYN") will also be large enough for key pair generation ("KEYGEN"),
 * dynamic integer-only signing ("SIGNDYN_NTT") and any signature verification
 * ("VERIFY" / "VERIFY_NTT").
 *
 * Here are the actual values for the temporary buffer sizes (in bytes):
 *
 *  deg  sigD  kgen   ver  verN sigDN   sig mkpub
 *    2   115   103   101    77    63    61    51
 *    4   221   199   193   145   121   113    99
 *    8   433   391   377   281   237   217   195
 *   16   859   775   747   555   471   427   387
 *   32  1711  1543  1487  1103   939   847   771
 *   64  3415  3079  2967  2199  1875  1687  1539
 *  128  6823  6151  5927  4391  3747  3367  3075
 *  256 13639 12295 11847  8775  7491  6727  6147
 *  512 27271 24583 23687 17543 14979 13447 12291
 * 1024 54535 49159 47367 35079 29955 26887 24579
 *
 * Here, first the scheme Hawk is listed. The suffix 'D' stands for dynamic
 * signing as opposed to batched signing where you need to do some
 * preprocessing on the secret key to get the expanded key. The suffix 'N'
 * stands for number theoretic transform and these functions will be used by
 * the reference implementation as these functions do not make use of floating
 * point numbers while the AVX2 optimized version uses floats for faster
 * signing/verifying.
 *
 *
 * KEY SIZES
 *
 * Public and secret keys are exchanged as serialized sequences of bytes.
 * The secret key size is fixed for a given degree and the HAWK_SECKEY_SIZE
 * macro returns this value as constant expressions.
 * For a given degree, the public key size follows a distribution that is close
 * to a normal distribution. The average number of bytes of a public key
 * together with the standard deviation is reported below. The HAWK_PUBKEY_SIZE
 * array contains for all degrees (logn = 1 to 10) a number that is six
 * standard deviations above the average, so a generated public key will be
 * larger than this number only once in a million times. The public key is
 * padded with zeros to the size in HAWK_PUBKEY_SIZE.
 *
 * Below is the average size in bytes of the public key together with standard
 * deviation and the minimum and maximum size encountered measured over 100,000
 * key generations. The column 6sigma gives avg + 6 stddev which when rounded
 * to the nearest integer are the values found in the HAWK_PUBKEY_SIZE array.
 * Note: these numbers include the header byte which contains the degree.
 *
 *  deg     avg+/-stddev min  max  6sigma
 *    2    4.67+/- 0.49    4    8    7.62
 *    4    6.59+/- 0.56    6   12    9.96
 *    8   11.40+/- 0.67   10   18   15.45
 *   16   22.05+/- 0.83   20   28   27.02
 *   32   46.01+/- 1.35   43   55   54.09
 *   64   98.61+/- 1.61   94  112  108.25
 *  128  214.12+/- 2.78  206  232  230.82
 *  256  464.55+/- 3.43  454  484  485.14
 *  512 1006.51+/- 6.09  987 1035 1043.02
 * 1024 2329.17+/-11.18 2291 2392 2396.25
 *
 *
 * FORMATS
 *
 * There are two formats for signatures:
 *
 *   - COMPACT: this is the default format, which yields the shortest
 *     signatures on average. However, the size is variable though within a
 *     limited range (see below).
 *
 *   - PADDED: this is the compact format, but with extra padding bytes
 *     to obtain a fixed size known at compile-time. The size depends only
 *     on the degree; the HAWK_SIG_PADDED_SIZE macro computes it. The
 *     signature process enforces that size by signing with a fresh salt until
 *     an appropriate size is obtained (such restarts are uncommon enough that
 *     the computational overhead is negligible).
 *
 * The signature format is selected by the 'sig_type' parameter to
 * the signature generation and verification functions.
 *
 * Actual signature size has been measured over 1,000,000 signatures for
 * degrees 512 and 1024 (10,000 random keys, 100 signatures per key):
 *
 *   deg  Gol.Rice enc   huffman enc
 *   512   541.7+/-4.4   539.1+/-4.4
 *  1024  1194.7+/-5.7  1189.4+/-5.9
 *
 * Here the Golomb-Rice encodings output 5 and 6 of the lowest significant bits
 * in binary for Hawk-512 and Hawk-1024 respectively. All lengths are in bytes
 * and include the header byte and salt.
 *
 * A secret key, in its encoded format, can be used as parameter to
 * hawk_sign_dyn(). An "expanded secret key" is computed with
 * hawk_expand_seckey(), to be used with hawk_sign(). The
 * expanded secret key is much larger than the encoded secret key, and
 * its format is not portable. Its size (in bytes) is provided by
 * HAWK_EXPANDEDKEY_SIZE. There are no specific alignment requirements
 * on expanded keys, except that the alignment of a given expanded key
 * must not change, i.e. if an expanded key is moved from address addr1
 * to address addr2 then it must hold that addr1 = addr2 (mod 8).
 * Expanded secret keys are meant to be used when several signatures are
 * to be computed with the same secret key: amortized cost per signature
 * is about halved when using expanded secret keys (for short messages,
 * and depending on underlying architecture and implementation choices).
 *
 *
 * USE OF SHAKE256
 *
 * SHAKE256 is used in two places:
 *
 *  - As a PRNG: all functions that require randomness (key pair
 *    generation, signature generation) receive as parameter a SHAKE256
 *    object, in output mode, from which pseudorandom data is obtained.
 *
 *    A SHAKE256 instance, to be used as a RNG, can be initialized
 *    from an explicit 48-byte seed, or from an OS-provided RNG. Using
 *    an explicit seed is meant for reproducibility of test vectors,
 *    or to be used in cases where no OS-provided RNG is available and
 *    supported.
 *
 *  - As the hashing mechanism for the message which should be signed.
 *    The streamed signature API exposes that SHAKE256 object, since
 *    the caller then performs the hashing externally.
 */

/* ========================================================================= */
/*
 * Error codes.
 *
 * Most functions in this API that may fail for some reason return an
 * 'int' value which will be 0 on success, or a negative error code.
 * The macros below define the error codes. In the interest of forward
 * compatibility, callers should be prepared to receive additional error
 * codes not included in the list below.
 */

/*
 * HAWK_ERR_RANDOM is returned when the library tries to use an
 * OS-provided RNG, but either none is supported, or that RNG fails.
 */
#define HAWK_ERR_RANDOM     -1

/*
 * HAWK_ERR_SIZE is returned when a buffer has been provided to
 * the library but is too small to receive the intended value.
 */
#define HAWK_ERR_SIZE       -2

/*
 * HAWK_ERR_FORMAT is returned when decoding of an external object
 * (public key, secret key, signature) fails.
 */
#define HAWK_ERR_FORMAT     -3

/*
 * HAWK_ERR_BADSIG is returned when verifying a signature, the signature
 * is validly encoded, but its value does not match the provided message
 * and public key.
 */
#define HAWK_ERR_BADSIG     -4

/*
 * HAWK_ERR_BADARG is returned when a provided parameter is not in
 * a valid range.
 */
#define HAWK_ERR_BADARG     -5

/*
 * HAWK_ERR_INTERNAL is returned when some internal computation failed.
 */
#define HAWK_ERR_INTERNAL   -6

/* ========================================================================= */
/*
 * Signature formats.
 */

/*
 * Variable-size signature. This format produces the most compact
 * signatures on average, but the signature size may vary depending
 * on secret key, signed data, and random seed.
 */
#define HAWK_SIG_COMPACT   1

/*
 * Fixed-size signature. This format produces is equivalent to the
 * "compact" format, but includes padding to a known fixed size
 * (specified by HAWK_SIG_PADDED_SIZE). With this format, the
 * signature generation loops until an appropriate signature size is
 * achieved (such looping is uncommon) and adds the padding bytes;
 * the verification functions check the presence and contents of the
 * padding bytes.
 */
#define HAWK_SIG_PADDED       2

/* ========================================================================= */
/*
 * Sizes.
 *
 * The sizes are expressed in bytes. Each size depends on the Hawk degree,
 * which is provided logarithmically: use logn=9 for Hawk-512 and logn=10 for
 * Hawk-1024. Valid values for logn range from 1 to 10 (values 1 to 8
 * correspond to reduced variants of Hawk that do not provided adequate
 * security and are meant for research purposes only).
 *
 * The sizes are provided as macros that evaluate to constant expressions, as
 * long as the 'logn' parameter is itself a constant expression. Moreover, all
 * sizes are monotonic (for each size category, increasing logn cannot result
 * in a shorter length).
 *
 * Note: each macro may evaluate its argument 'logn' several times.
 */

/*
 * Secret key size (in bytes). The size is exact.
 */
#define HAWK_SECKEY_SIZE(logn) \
	((logn) <= 1 ? 6u : (1u + (((logn) == 10 ? 10u : 9u) << ((logn) - 2))))

/*
 * Public key size (in bytes). The size is an upper bound on the allowed size,
 * the average is lower than this.
 */
extern const size_t HAWK_PUBKEY_SIZE[11];

/*
 * Currently, signature sizes are only determined for logn = 9 and logn = 10,
 * as these are determined by simulations in each case, giving an average size
 * and standard deviation. Currently the compact maxsize is avg + 12 stdddev,
 * while the padded size is avg + 6 stddev.
 *
 * The padded signature is therefore expected to fail with probability at most
 * 2 10^{-9} times and if this fails, signing is restarted (while reusing the
 * old hash).
 *
 * The compact signature may fail with probability at most 10^{-32} assuming
 * the sizes are normally distributed with the found average and standard
 * deviation.
 */

/*
 * Maximum practical signature size (in bytes) when using the COMPACT
 * format.
 */
#define HAWK_SIG_COMPACT_MAXSIZE(logn) \
	((logn) == 10 ? 1264u : 595u)

/*
 * Signature size (in bytes) when using the PADDED format. The size is exact.
 */
#define HAWK_SIG_PADDED_SIZE(logn) \
	((logn) == 10 ? 1229u : 568u)

/*
 * The number of bytes that are required for the salt that gets prepended to a
 * message before hashing.
 *
 * From GPV08, having \lambda security bits and a transcript size of Q_s,
 * having k = \lambda + \log_2(Q_s) salt bits gives a hash collision of message
 * to happen with probability at most 2^-k.
 * So for n = 512, take k = 192 and for n = 1024, take k = 320.
 */
#define HAWK_SALT_SIZE(logn) ((logn) <= 9 ? 24u : 40u)

/*
 * The number of bytes that are required for the hash when signing a message or
 * verifying a signature using the hash-then-sign principle on which HAWK is
 * based.
 */
#define HAWK_HASH_SIZE(logn) ((logn) <= 2 ? 2u : (1u << ((logn) - 2)))

/* ========================================================================= */
/*
 * Temporary buffer sizes.
 */

/*
 * Temporary buffer size for key pair generation.
 */
#define HAWK_TMPSIZE_KEYGEN(logn) \
	((48u << (logn)) + 7)

/*
 * Temporary buffer size for computing the public key from the secret key.
 */
#define HAWK_TMPSIZE_MAKEPUB(logn) \
	((24u << (logn)) + 3)

/*
 * Temporary buffer size for dynamically generating a signature.
 */
#define HAWK_TMPSIZE_SIGNDYN(logn) \
	(HAWK_HASH_SIZE(logn) + (53u << (logn)) + 7) // 5 in hawk.c, 6*8 in sign
#define HAWK_TMPSIZE_SIGNDYN_NTT(logn) \
	(HAWK_HASH_SIZE(logn) + (29u << (logn)) + 3) // 5 in hawk.c, 6*4 in sign

/*
 * Temporary buffer size for generating a signature with an expanded key.
 */
#define HAWK_TMPSIZE_SIGN(logn) \
	(HAWK_HASH_SIZE(logn) + (26u << (logn)) + 7)

/*
 * Size of an expanded secret key.
 */
#if HAWK_RECOVER_CHECK
/* For the recover-check, we also need to store 1/q00. */
#define HAWK_EXPANDEDKEY_SIZE(logn) \
	((36u << (logn)) + 8)
#else
/* Store f, g, F, G in FFT-representation */
#define HAWK_EXPANDEDKEY_SIZE(logn) \
	((32u << (logn)) + 8)
#endif

/*
 * Temporary buffer size for verifying a signature.
 */

#define HAWK_TMPSIZE_VERIFY(logn) \
	(HAWK_HASH_SIZE(logn) + (46u << (logn)) + 7)
#define HAWK_TMPSIZE_VERIFY_NTT(logn) \
	(HAWK_HASH_SIZE(logn) + (34u << (logn)) + 7)

/* ========================================================================= */
/*
 * SHAKE256.
 */

/*
 * Context for a SHAKE256 computation. Contents are opaque.
 * Contents are pure data with no pointer; they need not be released
 * explicitly and don't reference any other allocated resource. The
 * caller is responsible for allocating the context structure itself,
 * typically on the stack.
 */
typedef struct {
	uint64_t opaque_contents[26];
} shake256_context;

/*
 * Initialize a SHAKE256 context to its initial state. The state is
 * then ready to receive data (with shake256_inject()).
 */
void shake256_init(shake256_context *sc);

/*
 * Inject some data bytes into the SHAKE256 context ("absorb" operation).
 * This function can be called several times, to inject several chunks
 * of data of arbitrary length.
 */
void shake256_inject(shake256_context *sc, const void *data, size_t len);

/*
 * Flip the SHAKE256 state to output mode. After this call, shake256_inject()
 * can no longer be called on the context, but shake256_extract() can be
 * called.
 *
 * Flipping is one-way; a given context can be converted back to input
 * mode only by initializing it again, which forgets all previously
 * injected data.
 */
void shake256_flip(shake256_context *sc);

/*
 * Extract bytes from the SHAKE256 context ("squeeze" operation). The
 * context must have been flipped to output mode (with shake256_flip()).
 * Arbitrary amounts of data can be extracted, in one or several calls
 * to this function.
 */
void shake256_extract(shake256_context *sc, void *out, size_t len);

/*
 * Initialize a SHAKE256 context as a PRNG from the provided seed.
 * This initializes the context, injects the seed, then flips the context
 * to output mode to make it ready to produce bytes.
 */
void shake256_init_prng_from_seed(shake256_context *sc, const void *seed,
	size_t seed_len);

/*
 * Initialize a SHAKE256 context as a PRNG, using an initial seed from
 * the OS-provided RNG. If there is no known/supported OS-provided RNG,
 * or if that RNG fails, then the context is not properly initialized
 * and HAWK_ERR_RANDOM is returned.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int shake256_init_prng_from_system(shake256_context *sc);

/* ========================================================================= */
/*
 * Key pair generation.
 */

/*
 * Generate a new keypair.
 *
 * The logarithm of the Hawk degree (logn) must be in the 1 to 10
 * range; values 1 to 8 correspond to reduced versions of Hawk that do
 * not provide adequate security and are meant for research purposes
 * only.
 *
 * The source of randomness is the provided SHAKE256 context *rng, which
 * must have been already initialized, seeded, and set to output mode (see
 * shake256_init_prng_from_seed() and shake256_init_prng_from_system()).
 *
 * The new secret key is written in the buffer pointed to by seckey.
 * The size of that buffer must be specified in seckey_len; if that
 * size is too low, then this function fails with HAWK_ERR_SIZE. The
 * actual secret key length can be obtained from the HAWK_SECKEY_SIZE()
 * macro.
 *
 * If pubkey is not NULL, then the new public key is written in the buffer
 * pointed to by pubkey. The size of that buffer must be specified in
 * pubkey_len; if that size is too low, then this function fails with
 * HAWK_ERR_SIZE. The actual public key length can be obtained from the
 * HAWK_PUBKEY_SIZE() macro.
 *
 * If pubkey is NULL then pubkey_len is ignored; the secret key will
 * still be generated and written to seckey[], but the public key
 * won't be written anywhere. The public key can be recomputed later on
 * from the secret key with hawk_make_public().
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len
 * MUST be at least HAWK_TMPSIZE_KEYGEN(logn) bytes.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_keygen_make(shake256_context *rng, unsigned logn, void *seckey,
	size_t seckey_len, void *pubkey, size_t pubkey_len, void *tmp,
	size_t tmp_len);

/*
 * Recompute the public key from the secret key.
 *
 * The secret key is provided encoded. This function decodes the
 * secret key and verifies that its length (in bytes) is exactly
 * the provided value seckey_len (trailing extra bytes are not
 * tolerated).
 *
 * The public key is written in the buffer pointed to by pubkey. The
 * size (in bytes) of the pubkey buffer must be provided in pubkey_len;
 * if it is too short for the public key, then HAWK_ERR_SIZE is
 * returned. The actual public key size can be obtained from the
 * HAWK_PUBKEY_SIZE() macro.
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len
 * MUST be at least HAWK_TMPSIZE_MAKEPUB(logn) bytes.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_make_public(void *pubkey, size_t pubkey_len, const void *seckey,
	size_t seckey_len, void *tmp, size_t tmp_len);

/*
 * Get the Hawk degree from an encoded secret key, public key or
 * signature. Returned value is the logarithm of the degree (1 to 10),
 * or a negative error code.
 */
int hawk_get_logn(const void *obj, size_t len);

/* ========================================================================= */
/*
 * Signature generation.
 */

/*
 * Expand a secret key. The provided Hawk secret key (seckey, of size
 * seckey_len bytes) is decoded and expanded into expanded_key[].
 *
 * The expanded_key[] buffer has size expanded_key_len, which MUST be at least
 * HAWK_EXPANDEDKEY_SIZE(logn) bytes (where 'logn' qualifies the Hawk degree
 * encoded in the secret key and can be obtained with hawk_get_logn()).
 * Expanded key contents have an internal, implementation-specific format.
 * Expanded keys may be moved in RAM only if their 8-byte alignment remains
 * unchanged.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_expand_seckey(void *expanded_key, size_t expanded_key_len,
	const void *seckey, size_t seckey_len);

/* Helper for hawk_sign_dyn_finish */
int hawk_sign_dyn(shake256_context *rng, void *sig, size_t *sig_len,
	int sig_type, const void *seckey, size_t seckey_len, const void *data,
	size_t data_len, void *tmp, size_t tmp_len);

/* Helper for hawk_sign_finish */
int hawk_sign(shake256_context *rng, void *sig, size_t *sig_len,
	int sig_type, const void *expanded_key, const void *data, size_t data_len,
	void *tmp, size_t tmp_len);

/* ========================================================================= */
/*
 * Signature generation, streamed API.
 *
 * In the streamed API, the caller performs the data hashing externally. An
 * initialization function (hawk_sign_start()) is first called; it initializes
 * the SHAKE256 context. The caller must then inject the data to sign in the
 * SHAKE256 context, and finally call one of the hawk_sign_XXX_finish()
 * functions to finalize the signature generation.
 */

/*
 * Start a signature generation context: the *hash_data context is initialized.
 */
void hawk_sign_start(shake256_context *hash_data);

/*
 * Finish a signature generation operation, using the secret key held in
 * seckey[] of length seckey_len bytes. The hashed message is provided as the
 * SHAKE256 context *hash_data, which must still be in input mode (i.e. not yet
 * flipped to output mode). During signing, a salt is generated, written in
 * salt[] and added to (a copy of) *hash_data, after which output is taken from
 * it. The salt length is 24 or 40 bytes, see the HAWK_SALT_SIZE macro.
 *
 * The source of randomness is the provided SHAKE256 context *rng, which
 * must have been already initialized, seeded, and set to output mode (see
 * shake256_init_prng_from_seed() and shake256_init_prng_from_system()).
 *
 * The signature is written in sig[]. The caller must set *sig_len to the
 * maximum size of sig[]; if the signature computation is successful, then
 * *sig_len will be set to the actual length of the signature. The signature
 * length depends on the signature type, which is specified with the sig_type
 * parameter to one of the three defined values HAWK_SIG_COMPACT or
 * HAWK_SIG_PADDED; for the latter, the signature length is fixed (for a given
 * Hawk degree).
 *
 * Regardless of the signature type, the process is constant-time with regard
 * to the secret key.
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len MUST be
 * at least HAWK_TMPSIZE_SIGNDYN(logn) bytes when HAWK_AVX is defined, else at
 * least HAWK_TMPSIZE_SIGNDYN_NTT(logn).
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_sign_dyn_finish(shake256_context *rng, void *sig, size_t *sig_len,
	int sig_type, const void *seckey, size_t seckey_len,
	shake256_context *hash_data, void *salt, void *tmp, size_t tmp_len);

/*
 * Finish a signature generation operation, using the expanded secret key held
 * in expanded_key[] (as obtained from hawk_expand_seckey()). The hashed
 * message is provided as the SHAKE256 context *hash_data, which must still be
 * in input mode (i.e. not yet flipped to output mode). During signing, a salt
 * is generated, written in salt[] and added to (a copy of) *hash_data, after
 * which output is taken from it. The salt length is 24 or 40 bytes, see the
 * HAWK_SALT_SIZE macro.
 *
 * The source of randomness is the provided SHAKE256 context *rng, which must
 * have been already initialized, seeded, and set to output mode (see
 * shake256_init_prng_from_seed() and shake256_init_prng_from_system()).
 *
 * The signature is written in sig[]. The caller must set *sig_len to the
 * maximum size of sig[]; if the signature computation is successful, then
 * *sig_len will be set to the actual length of the signature. The signature
 * length depends on the signature type, which is specified with the sig_type
 * parameter to one of the three defined values HAWK_SIG_COMPACT or
 * HAWK_SIG_PADDED; for the latter, the signature length is fixed (for a given
 * Hawk degree).
 *
 * Regardless of the signature type, the process is constant-time with regard
 * to the secret key.
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len MUST be
 * at least HAWK_TMPSIZE_SIGN(logn) bytes.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_sign_finish(shake256_context *rng, void *sig, size_t *sig_len,
	int sig_type, const void *expanded_key, shake256_context *hash_data,
	void *salt, void *tmp, size_t tmp_len);

/* ========================================================================= */
/*
 * Signature verification.
 */

/* Helper for hawk_verify_finish */
int hawk_verify(const void *sig, size_t sig_len, int sig_type,
	const void *pubkey, size_t pubkey_len, const void *data, size_t data_len,
	void *tmp, size_t tmp_len);


/*
 * Start a streamed signature verification. The provided SHAKE256 context
 * *hash_data is initialized. The caller shall then inject the message data
 * into the SHAKE256 context, and finally call hawk_verify_finish().
 */
void hawk_verify_start(shake256_context *hash_data);

/*
 * Finish a streamed signature verification. The signature sig[] (of length
 * sig_len bytes) is verified against the provided public key pubkey[] (of
 * length pubkey_len bytes) and the hashed message. The hashed message is
 * provided as a SHAKE256 context *hash_data; that context must have received
 * the message itself, and still be in input mode (not yet flipped to output
 * mode). *hash_data is modified by the verification process, as salt is added
 * from the signature.
 *
 * The sig_type parameter must be zero, or one of HAWK_SIG_COMPACT or
 * HAWK_SIG_PADDED. The function verifies that the provided signature has the
 * correct format. If sig_type is zero, then the signature format is inferred
 * from the signature header byte; note that in that case, the signature is
 * malleable (since a signature value can be transcoded to other formats).
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len MUST be
 * at least HAWK_TMPSIZE_VERIFY(logn) bytes.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_verify_finish(const void *sig, size_t sig_len, int sig_type,
	const void *pubkey, size_t pubkey_len, shake256_context *hash_data,
	void *tmp, size_t tmp_len);

/* ========================================================================= */

#ifdef __cplusplus
}
#endif

#endif // HAWK_H__
