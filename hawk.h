/*
 * External Hawk API.
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

/* ==================================================================== */
/*
 * Hawk API Notes
 * ----------------
 *
 *
 * HAWK DEGREE
 *
 * Hawk is parameterized by a degree, which is a power of two. Formally,
 * one value is possible: 512 for Hawk-512. This implementation also supports
 * lower degrees, from 2 to 256; these reduced variants do not provide adequate
 * security and should be used for research purposes only.
 *
 * In all functions and macros defined below, the degree is provided
 * logarithmically as the 'logn' parameter: logn ranges from 1 to 9,
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
 * possibly distinct from each other. For all degrees from 8 to 512
 * (logn = 3 to 9), the following sizes are in ascending order:
 *
 *    HAWK_TMPSIZE_MAKEPUB
 *    HAWK_TMPSIZE_VERIFY
 *    HAWK_TMPSIZE_KEYGEN
 *    HAWK_TMPSIZE_SIGN
 *    HAWK_TMPSIZE_EXPANDPRIV
 *    HAWK_TMPSIZE_SIGNDYN
 *
 * i.e. a temporary buffer large enough for computing signatures with
 * an expanded key ("SIGN") will also be large enough for a
 * key pair generation ("KEYGEN"). For logn = 1 or 2, the same order
 * holds, except that the KEYGEN buffer is larger.
 *
 * Here are the actual values for the temporary buffer sizes (in bytes):
 *
 * TODO: update table
 * degree  mkpub  verify  keygen      sign  expkey  signdyn
 *     2      13      17     285       107     111      163
 *     4      25      33     291       207     215      319
 *     8      49      65     303       407     423      631
 *    16      97     129     503       807     839     1255
 *    32     193     257     999      1607    1671     2503
 *    64     385     513    1991      3207    3335     4999
 *   128     769    1025    3975      6407    6663     9991
 *   256    1537    2049    7943     12807   13319    19975
 *   512    3073    4097   15879     25607   26631    39943
 *
 * Take care that the "expkey" column here qualifies the temporary buffer
 * for the key expansion process, but NOT the expanded key itself (which
 * has size HAWK_EXPANDEDKEY_SIZE(logn) and is larger than that).
 *
 *
 * FORMATS
 *
 * Public and private keys are exchanged using a compressed Gaussian encoding,
 * part of compress.c. Their respective sizes (for a given degree) follow a
 * distribution that is close to a normal distribution.
 *
 * The HAWK_MIN_SECKEY_SIZE and HAWK_MAX_SECKEY_SIZE arrays give the lower and
 * upper bound on the encoded secret key size in bytes.
 *
 * The HAWK_MIN_PUBKEY_SIZE and HAWK_MAX_PUBKEY_SIZE arrays give the lower and
 * upper bound on the encoded public key size in bytes.
 *
 * Below are the average sizes in number of bytes of secret keys and public
 * keys, taken over 10000 samples. To arrive at the minimal and maximal sizes,
 * add -,+5 std.dev. to the average and round to the average.
 *
 * TODO: run test_sizes with fixed compress.c, or where flushing last byte only happens at the end...
 * TODO: also account for the 1 header byte containing type of key & logn.
 * logn | Average +/- stddev | min  | max
 * -----+- Secret key -------+------+-----
 *    1 |    3.45 +/-   0.64 |    3 |    7
 *    2 |    7.49 +/-   0.90 |    6 |   11
 *    3 |   15.49 +/-   1.33 |   12 |   21
 *    4 |   29.63 +/-   1.50 |   25 |   36
 *    5 |   57.83 +/-   2.49 |   50 |   69
 *    6 |  106.48 +/-   2.88 |   96 |  118
 *    7 |  203.70 +/-   4.94 |  187 |  223
 *    8 |  367.20 +/-   5.64 |  346 |  391
 *    9 |  678.76 +/-   6.66 |  653 |  703
 * -----+- Public key -------+------+-----
 *    1 |    5.27 +/-   0.52 |    0 |    8
 *    2 |    9.59 +/-   0.71 |    9 |   14
 *    3 |   18.00 +/-   1.03 |   16 |   26
 *    4 |   34.46 +/-   1.36 |   31 |   42
 *    5 |   66.70 +/-   2.06 |   60 |   77
 *    6 |  130.40 +/-   2.69 |  122 |  141
 *    7 |  255.17 +/-   4.17 |  242 |  274
 *    8 |  501.60 +/-   5.32 |  485 |  528
 *    9 |  984.52 +/-   8.32 |  955 | 1019
 *
 *
 * There are two formats for signatures:
 *
 *   - COMPRESSED: this is the default format, which yields the shortest
 *     signatures on average. However, the size is variable (see below)
 *     though within a limited range.
 *
 *   - PADDED: this is the compressed format, but with extra padding bytes
 *     to obtain a fixed size known at compile-time. The size depends only
 *     on the degree; the HAWK_SIG_PADDED_SIZE macro computes it. The
 *     signature process enforces that size by restarting the process
 *     until an appropriate size is obtained (such restarts are uncommon
 *     enough that the computational overhead is negligible).
 *
 * The signature format is selected by the 'sig_type' parameter to
 * the signature generation and verification functions.
 *
 * Actual signature size has been measured over 10000 signatures for each
 * degree (100 random keys, 100 signatures per key):
 *
 * TODO: update table.
 * degree     ct   padded  compressed (with std. dev)  comp_max
 *     2      44      44       44.00 (+/- 0.00)            44
 *     4      47      47       46.03 (+/- 0.17)            47
 *     8      52      52       50.97 (+/- 0.26)            52
 *    16      65      63       60.45 (+/- 0.52)            64
 *    32      89      82       79.53 (+/- 0.68)            86
 *    64     137     122      117.69 (+/- 0.94)           130
 *   128     233     200      193.96 (+/- 1.30)           219
 *   256     425     356      346.53 (+/- 1.84)           397
 *   512     809     666      651.59 (+/- 2.55)           752
 *
 * with:
 *   degree = Hawk degree = 2^logn
 *   padded = HAWK_SIG_PADDED_SIZE(logn)  (size of a PADDED signature)
 *   compressed = measured average length of a COMPRESSED signature
 *   v_max = HAWK_SIG_COMPRESSED_MAXSIZE(logn)  (maximum theoretical
 *           size of a COMPRESSED signature)
 * All lengths are in bytes.
 *
 * A private key, in its encoded format, can be used as parameter to
 * hawk_sign_dyn(). An "expanded private key" is computed with
 * hawk_expand_privkey(), to be used with hawk_sign(). The
 * expanded private key is much larger than the encoded private key, and
 * its format is not portable. Its size (in bytes) is provided by
 * HAWK_EXPANDEDKEY_SIZE. There are no specific alignment requirements
 * on expanded keys, except that the alignment of a given expanded key
 * must not change (i.e. if an expanded key is moved from address addr1
 * to address addr2, then it must hold that addr1 = addr2 mod 8).
 * Expanded private keys are meant to be used when several signatures are
 * to be computed with the same private key: amortized cost per signature
 * is about halved when using expanded private keys (for short messages,
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

/* ==================================================================== */
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
 * (public key, private key, signature) fails.
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

/* ==================================================================== */
/*
 * Signature formats.
 */

/*
 * Variable-size signature. This format produces the most compact
 * signatures on average, but the signature size may vary depending
 * on private key, signed data, and random seed.
 */
#define HAWK_SIG_COMPRESSED   1

/*
 * Fixed-size signature. This format produces is equivalent to the
 * "compressed" format, but includes padding to a known fixed size
 * (specified by HAWK_SIG_PADDED_SIZE). With this format, the
 * signature generation loops until an appropriate signature size is
 * achieved (such looping is uncommon) and adds the padding bytes;
 * the verification functions check the presence and contents of the
 * padding bytes.
 */
#define HAWK_SIG_PADDED       2

/* ==================================================================== */
/*
 * Sizes.
 *
 * The sizes are expressed in bytes. Each size depends on the Hawk
 * degree, which is provided logarithmically: use logn=9 for Hawk-512. Valid
 * values for logn range from 1 to 9 (values 1 to 8 correspond to reduced
 * variants of Hawk that do not provided adequate security and are meant for
 * research purposes only).
 *
 * The sizes are provided as macros that evaluate to constant
 * expressions, as long as the 'logn' parameter is itself a constant
 * expression. Moreover, all sizes are monotonic (for each size category,
 * increasing logn cannot result in a shorter length).
 *
 * Note: each macro may evaluate its argument 'logn' several times.
 */

extern const size_t HAWK_MIN_SECKEY_SIZE[10];
extern const size_t HAWK_MAX_SECKEY_SIZE[10];
extern const size_t HAWK_MIN_PUBKEY_SIZE[10];
extern const size_t HAWK_MAX_PUBKEY_SIZE[10];


/*
 * Maximum signature size (in bytes) when using the COMPRESSED format.
 * In practice, the signature will be shorter.
 */
#define HAWK_SIG_COMPRESSED_MAXSIZE(logn) \
	(((((11u << (logn)) + (101u >> (10 - (logn)))) \
	+ 7) >> 3) + 41)

/*
 * Signature size (in bytes) when using the PADDED format. The size
 * is exact.
 */
#define HAWK_SIG_PADDED_SIZE(logn) \
	(44u + 3 * (256u >> (10 - (logn))) + 2 * (128u >> (10 - (logn))) \
	+ 3 * (64u >> (10 - (logn))) + 2 * (16u >> (10 - (logn))) \
	- 2 * (2u >> (10 - (logn))) - 8 * (1u >> (10 - (logn))))

/*
 * Temporary buffer size for key pair generation.
 */
#define HAWK_TMPSIZE_KEYGEN(logn) \
	(((logn) <= 3 ? 272u : (48u << (logn))) + (28u << (logn)) + 7)

/*
 * Temporary buffer size for computing the public key from the private key.
 */
/* #define HAWK_TMPSIZE_MAKEPUB(logn) \
	((6u << (logn)) + 1) */

/*
 * Temporary buffer size for generating a signature with an expanded key.
 */
#define HAWK_TMPSIZE_SIGN(logn) \
	((logn) <= 3 ? 1U : (1U << ((logn)-3))) + (2u << (logn)) + (26u << (logn)) + 7
	// ((50u << (logn)) + 7)

/*
 * Temporary buffer size for expanding a private key.
 */
#define HAWK_TMPSIZE_EXPANDPRIV(logn) \
	(3u << (logn))

/*
 * Size of an expanded private key.
 */
#define HAWK_EXPANDEDKEY_SIZE(logn) \
	((36u << (logn)) + 8)

/*
 * Temporary buffer size for verifying a signature.
 */
#define HAWK_TMPSIZE_VERIFY(logn) \
	(logn <= 3 ? 1U : (1U << (logn-3))) + (26u << (logn)) + (16u << (logn))

/* ==================================================================== */
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

/* ==================================================================== */
/*
 * Key pair generation.
 */

/*
 * Generate a new keypair.
 *
 * The logarithm of the Hawk degree (logn) must be in the 1 to 9
 * range; values 1 to 8 correspond to reduced versions of Hawk that do
 * not provide adequate security and are meant for research purposes
 * only.
 *
 * The source of randomness is the provided SHAKE256 context *rng, which
 * must have been already initialized, seeded, and set to output mode (see
 * shake256_init_prng_from_seed() and shake256_init_prng_from_system()).
 *
 * The new private key is written in the buffer pointed to by privkey.
 * The size of that buffer must be specified in privkey_len; if that
 * size is too low, then this function fails with HAWK_ERR_SIZE. The
 * actual private key length can be obtained from the HAWK_SECKEY_SIZE()
 * macro.
 *
 * If pubkey is not NULL, then the new public key is written in the buffer
 * pointed to by pubkey. The size of that buffer must be specified in
 * pubkey_len; if that size is too low, then this function fails with
 * HAWK_ERR_SIZE. The actual public key length can be obtained from the
 * HAWK_PUBKEY_SIZE() macro.
 *
 * If pubkey is NULL then pubkey_len is ignored; the private key will
 * still be generated and written to privkey[], but the public key
 * won't be written anywhere. The public key can be recomputed later on
 * from the private key with hawk_make_public().
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len
 * MUST be at least HAWK_TMPSIZE_KEYGEN(logn) bytes.
 *
 * Returned value: 0 on success, or a negative error code.
 *
 *
 * Note/TODO: currently one MUST provide a public key of sufficient large size.
 * TODO: update documentation: *privkey and *pubkey get changed to the actual
 * length of secret/private key respectively.
 */
int hawk_keygen_make(shake256_context *rng, unsigned logn, void *privkey,
	size_t *privkey_len, void *pubkey, size_t *pubkey_len, void *tmp,
	size_t tmp_len);

/*
 * Recompute the public key from the private key.
 *
 * The private key is provided encoded. This function decodes the
 * private key and verifies that its length (in bytes) is exactly
 * the provided value privkey_len (trailing extra bytes are not
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
// int hawk_make_public(void *pubkey, size_t pubkey_len, const void *privkey,
//	size_t privkey_len, void *tmp, size_t tmp_len);

/*
 * Get the Hawk degree from an encoded private key, public key or
 * signature. Returned value is the logarithm of the degree (1 to 9),
 * or a negative error code.
 */
int hawk_get_logn(void *obj, size_t len);

/* ==================================================================== */
/*
 * Signature generation.
 */

/*
 * Sign the data provided in buffer data[] (of length data_len bytes),
 * using the private key held in privkey[] (of length privkey_len bytes).
 *
 * The source of randomness is the provided SHAKE256 context *rng, which
 * must have been already initialized, seeded, and set to output mode (see
 * shake256_init_prng_from_seed() and shake256_init_prng_from_system()).
 *
 * The signature is written in sig[]. The caller must set *sig_len to
 * the maximum size of sig[]; if the signature computation is
 * successful, then *sig_len will be set to the actual length of the
 * signature. The signature length depends on the signature type,
 * which is specified with the sig_type parameter to one of the three
 * defined values HAWK_SIG_COMPRESSED or HAWK_SIG_PADDED; for the latter, the
 * signature length is fixed (for a given Hawk degree).
 *
 * Regardless of the signature type, the process is constant-time with regard
 * to the private key.
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len
 * MUST be at least HAWK_TMPSIZE_SIGNDYN(logn) bytes.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_sign_dyn(shake256_context *rng, void *sig, size_t *sig_len,
	int sig_type, const void *privkey, size_t privkey_len, const void *data,
	size_t data_len, void *tmp, size_t tmp_len);

/*
 * Expand a private key. The provided Hawk private key (privkey, of
 * size privkey_len bytes) is decoded and expanded into expanded_key[].
 *
 * The expanded_key[] buffer has size expanded_key_len, which MUST be at
 * least HAWK_EXPANDEDKEY_SIZE(logn) bytes (where 'logn' qualifies the
 * Hawk degree encoded in the private key and can be obtained with
 * hawk_get_logn()). Expanded key contents have an internal,
 * implementation-specific format. Expanded keys may be moved in RAM
 * only if their 8-byte alignment remains unchanged.
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len
 * MUST be at least HAWK_TMPSIZE_EXPANDPRIV(logn) bytes.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_expand_privkey(void *expanded_key, size_t expanded_key_len,
	const void *privkey, size_t privkey_len, void *tmp, size_t tmp_len);

/*
 * Sign the data provided in buffer data[] (of length data_len bytes),
 * using the expanded private key held in expanded_key[], as generated
 * by hawk_expand_privkey().
 *
 * The source of randomness is the provided SHAKE256 context *rng, which
 * must have been already initialized, seeded, and set to output mode (see
 * shake256_init_prng_from_seed() and shake256_init_prng_from_system()).
 *
 * The signature is written in sig[]. The caller must set *sig_len to
 * the maximum size of sig[]; if the signature computation is
 * successful, then *sig_len will be set to the actual length of the
 * signature. The signature length depends on the signature type,
 * which is specified with the sig_type parameter to one of the three
 * defined values HAWK_SIG_COMPRESSED or HAWK_SIG_PADDED; for the latter, the
 * signature length is fixed (for a given Hawk degree).
 *
 * Regardless of the signature type, the process is constant-time with
 * regard to the private key.
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len
 * MUST be at least HAWK_TMPSIZE_SIGN(logn) bytes.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_sign(shake256_context *rng, void *sig, size_t *sig_len,
	int sig_type, const void *expanded_key, const void *data, size_t data_len,
	void *tmp, size_t tmp_len);

/* ==================================================================== */
/*
 * Signature generation, streamed API.
 *
 * In the streamed API, the caller performs the data hashing externally.
 * An initialization function (hawk_sign_start()) is first called; it
 * generates and returns a random 40-byte salt value; it also initializes
 * a SHAKE256 context and injects the salt value in that context. The
 * caller must then inject the data to sign in the SHAKE256 context, and
 * finally call hawk_sign_dyn_finish() or hawk_sign_finish() to
 * finalize the signature generation.
 */

/*
 * Start a signature generation context.
 *
 * A 40-byte salt is generated and written in salt[]. The *hash_data
 * context is also initialized, and the salt is injected in that context.
 *
 * The source of randomness is the provided SHAKE256 context *rng, which
 * must have been already initialized, seeded, and set to output mode (see
 * shake256_init_prng_from_seed() and shake256_init_prng_from_system()).
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_sign_start(shake256_context *rng, void *salt,
	shake256_context *hash_data);

/*
 * Finish a signature generation operation, using the private key held
 * in privkey[] (of length privkey_len bytes). The hashed salt + message
 * is provided as the SHAKE256 context *hash_data, which must still be
 * in input mode (i.e. not yet flipped to output mode). That context is
 * modified in the process.
 *
 * The salt value (which was used at the start of the hashing process,
 * usually as part of a hawk_sign_start() call) must be provided again,
 * because it is encoded into the signature. The salt length is 40 bytes.
 *
 * The source of randomness is the provided SHAKE256 context *rng, which
 * must have been already initialized, seeded, and set to output mode (see
 * shake256_init_prng_from_seed() and shake256_init_prng_from_system()).
 *
 * The signature is written in sig[]. The caller must set *sig_len to
 * the maximum size of sig[]; if the signature computation is
 * successful, then *sig_len will be set to the actual length of the
 * signature. The signature length depends on the signature type,
 * which is specified with the sig_type parameter to one of the three
 * defined values HAWK_SIG_COMPRESSED or HAWK_SIG_PADDED; for the latter, the
 * signature length is fixed (for a given Hawk degree).
 *
 * Regardless of the signature type, the process is constant-time with
 * regard to the private key.
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len
 * MUST be at least HAWK_TMPSIZE_SIGNDYN(logn) bytes.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_sign_dyn_finish(shake256_context *rng, void *sig, size_t *sig_len,
	int sig_type, const void *privkey, size_t privkey_len,
	shake256_context *hash_data, const void *salt, void *tmp, size_t tmp_len);

/*
 * Finish a signature generation operation, using the expanded private
 * key held in expanded_key[] (as obtained from
 * hawk_expand_privkey()). The hashed salt + message is provided as
 * the SHAKE256 context *hash_data, which must still be in input mode
 * (i.e. not yet flipped to output mode). That context is modified in
 * the process.
 *
 * The salt value (which was used at the start of the hashing process,
 * usually as part of a hawk_sign_start() call) must be provided again,
 * because it is encoded into the signature. The salt length is 40 bytes.
 *
 * The source of randomness is the provided SHAKE256 context *rng, which
 * must have been already initialized, seeded, and set to output mode (see
 * shake256_init_prng_from_seed() and shake256_init_prng_from_system()).
 *
 * The signature is written in sig[]. The caller must set *sig_len to
 * the maximum size of sig[]; if the signature computation is
 * successful, then *sig_len will be set to the actual length of the
 * signature. The signature length depends on the signature type,
 * which is specified with the sig_type parameter to one of the three
 * defined values HAWK_SIG_COMPRESSED or HAWK_SIG_PADDED; for the latter, the
 * signature length is fixed (for a given Hawk degree).
 *
 * Regardless of the signature type, the process is constant-time with
 * regard to the private key.
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len
 * MUST be at least HAWK_TMPSIZE_SIGN(logn) bytes.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_sign_finish(shake256_context *rng, void *sig, size_t *sig_len,
	int sig_type, const void *expanded_key, shake256_context *hash_data,
	const void *salt, void *tmp, size_t tmp_len);

/* ==================================================================== */
/*
 * Signature verification.
 */

/*
 * Verify the signature sig[] (of length sig_len bytes) with regards to
 * the provided public key pubkey[] (of length pubkey_len bytes) and the
 * message data[] (of length data_len bytes).
 *
 * The sig_type parameter must be zero, or one of HAWK_SIG_COMPRESSED or
 * HAWK_SIG_PADDED. The function verifies that the provided signature has the
 * correct format. If sig_type is zero, then the signature format is inferred
 * from the signature header byte; note that in that case, the signature is
 * malleable (since a signature value can be transcoded to other formats).
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len
 * MUST be at least HAWK_TMPSIZE_VERIFY(logn) bytes.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_verify(const void *sig, size_t sig_len, int sig_type,
	const void *pubkey, size_t pubkey_len, const void *data, size_t data_len,
	void *tmp, size_t tmp_len);

/*
 * Start a streamed signature verification. The provided SHAKE256 context
 * *hash_data is initialized, and the salt value (extracted from the
 * signature) is injected into it. The caller shall then inject the
 * message data into the SHAKE256 context, and finally call
 * hawk_verify_finish().
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_verify_start(shake256_context *hash_data, const void *sig,
	size_t sig_len);

/*
 * Finish a streamed signature verification. The signature sig[] (of
 * length sig_len bytes) is verified against the provided public key
 * pubkey[] (of length pubkey_len bytes) and the hashed message. The
 * hashed message is provided as a SHAKE256 context *hash_data;
 * that context must have received the salt and the message itself
 * (usually, the context is initialized and the salt injected as
 * part of a hawk_verify_start() call), and still be in input
 * mode (not yet flipped to output mode). *hash_data is modified by
 * the verification process.
 *
 * The sig_type parameter must be zero, or one of HAWK_SIG_COMPRESSED or
 * HAWK_SIG_PADDED. The function verifies that the provided signature has the
 * correct format. If sig_type is zero, then the signature format is inferred
 * from the signature header byte; note that in that case, the signature is
 * malleable (since a signature value can be transcoded to other formats).
 *
 * The tmp[] buffer is used to hold temporary values. Its size tmp_len
 * MUST be at least HAWK_TMPSIZE_VERIFY(logn) bytes.
 *
 * Returned value: 0 on success, or a negative error code.
 */
int hawk_verify_finish(const void *sig, size_t sig_len, int sig_type,
	const void *pubkey, size_t pubkey_len, shake256_context *hash_data,
	void *tmp, size_t tmp_len);

/* ==================================================================== */

#ifdef __cplusplus
}
#endif

#endif // HAWK_H__
