/*
 * Implementation of the external Hawk API.
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

#include "hawk.h"
#include "inner.h"

#ifdef __AVX2__
/*
 * Use floating points when the AVX2 instruction set is available, since a CPU
 * supporting AVX2 probably has a FPU.
 */
#define HAWK_AVX
#endif

/* see hawk.h */
const size_t HAWK_PUBKEY_SIZE[11] = {
	0 /* unused */, 8, 10, 15, 27, 54, 108, 231, 485, 1043, 2396
};

/* ========================================================================= */

/* see hawk.h */
void
shake256_init(shake256_context *sc)
{
	inner_shake256_init((inner_shake256_context *)sc);
}

/* see hawk.h */
void
shake256_inject(shake256_context *sc, const void *data, size_t len)
{
	inner_shake256_inject((inner_shake256_context *)sc, data, len);
}

/* see hawk.h */
void
shake256_flip(shake256_context *sc)
{
	inner_shake256_flip((inner_shake256_context *)sc);
}

/* see hawk.h */
void
shake256_extract(shake256_context *sc, void *out, size_t len)
{
	inner_shake256_extract((inner_shake256_context *)sc, out, len);
}

/* see hawk.h */
void
shake256_init_prng_from_seed(shake256_context *sc, const void *seed,
	size_t seed_len)
{
	shake256_init(sc);
	shake256_inject(sc, seed, seed_len);
}

/* see hawk.h */
int
shake256_init_prng_from_system(shake256_context *sc)
{
	uint8_t seed[48];

	if (!Zf(get_seed)(seed, sizeof seed)) {
		return HAWK_ERR_RANDOM;
	}
	shake256_init(sc);
	shake256_inject(sc, seed, sizeof seed);
	return 0;
}

/* ========================================================================= */

static inline int16_t *
align_i16(void *tmp)
{
	uint8_t *atmp;

	atmp = tmp;
	if (((uintptr_t)atmp & 1u) != 0) {
		atmp ++;
	}
	return (int16_t *)atmp;
}

static inline int32_t *
align_i32(void *tmp)
{
	uint8_t *atmp;
	unsigned off;

	atmp = tmp;
	off = (uintptr_t)atmp & 3u;
	if (off != 0) {
		atmp += 4u - off;
	}
	return (int32_t *)atmp;
}

static inline fpr *
align_fpr(void *tmp)
{
	uint8_t *atmp;
	unsigned off;

	atmp = tmp;
	off = (uintptr_t)atmp & 7u;
	if (off != 0) {
		atmp += 8u - off;
	}
	return (fpr *)atmp;
}

/* ========================================================================= */

/* see hawk.h */
int
hawk_keygen_make(shake256_context *rng, unsigned logn, void *seckey,
	size_t seckey_len, void *pubkey, size_t pubkey_len, void *tmp,
	size_t tmp_len)
{
	int8_t *f, *g, *F, *G;
	int16_t *iq00, *iq01;
	uint8_t *sk, *pk, *atmp;
	size_t u, n, sk_len, pk_len;
	unsigned oldcw;

	/*
	 * Check parameters.
	 */
	if (logn < 1 || logn > 10) {
		return HAWK_ERR_BADARG;
	}

	/*
	 * Check that the seckey and pubkey buffers are at least as large as the
	 * allowed encoded sizes, and check the temporary buffer size.
	 */
	if (seckey_len != HAWK_SECKEY_SIZE(logn)
		|| (pubkey != NULL && pubkey_len < HAWK_PUBKEY_SIZE[logn])
		|| tmp_len < HAWK_TMPSIZE_KEYGEN(logn)) {
		return HAWK_ERR_SIZE;
	}

	/*
	 * Prepare buffers and generate secret key.
	 * The buffers for iq00 and iq01 overlap with those of f, g, F, G but this
	 * is fine as we encode the public key after the secret key is already
	 * encoded.
	 */
	n = MKN(logn);

	f = (int8_t *)tmp;
	g = f + n;
	F = g + n;
	G = F + n;
	iq00 = align_i16(G + n);
	iq01 = iq00 + n;
	atmp = (uint8_t *)align_fpr(iq01 + n);

	/*
	 * Fix the first byte of secret key and secret key.
	 */
	sk = seckey;
	sk[0] = 0x50 + logn;

	if (pubkey == NULL) {
		/* Should not be possible */
		if (tmp_len < (8u << logn) + 7 + HAWK_PUBKEY_SIZE[logn]) {
			return HAWK_ERR_SIZE;
		}
		pk = atmp;
	} else {
		pk = pubkey;
	}
	pk[0] = 0x00 + logn;

	do {
		oldcw = set_fpu_cw(2);
		Zf(keygen)((inner_shake256_context *)rng,
			f, g, F, G, iq00, iq01, logn, atmp);
		set_fpu_cw(oldcw);

		/*
		 * For constant-time code, the secret key has a simple encoding of a
		 * fixed size.
		 */
		sk_len = Zf(encode_seckey)(sk + 1, HAWK_SECKEY_SIZE(logn) - 1,
			f, g, F, logn);
		if (sk_len != HAWK_SECKEY_SIZE(logn) - 1) {
			return HAWK_ERR_INTERNAL;
		}

		pk_len = Zf(encode_pubkey)(pk + 1, HAWK_PUBKEY_SIZE[logn] - 1,
			iq00, iq01, logn);

		/*
		 * Retry key-generation as the secret key or public key cannot be
		 * encoded, is shorter or is larger than the allowed size. This only
		 * happens with negligible probability.
		 */
	} while (pk_len == 0);

	if (pubkey != NULL) {
		for (u = 1; u < pubkey_len; u ++) {
			/*
			 * Pad the public key with zeros in constant time, by doing:
			 * if (u >= 1 + pk_len) pk[u] = 0;
			 */
			pk[u] &= -(((uint32_t)(u - 1 - pk_len)) >> 31);
		}
	}

	return 0;
}

/* see hawk.h */
int
hawk_make_public(void *pubkey, size_t pubkey_len, const void *seckey,
	size_t seckey_len, void *tmp, size_t tmp_len)
{
	unsigned logn;
	size_t u, n, pk_len;
	uint8_t *pk, *atmp;
	const uint8_t *sk;
	int8_t *f, *g, *F;
	int16_t *iq00, *iq01;

	/*
	 * Get degree from secret key header byte, and check parameters.
	 */
	if (seckey_len == 0) {
		return HAWK_ERR_FORMAT;
	}
	sk = seckey;
	if ((sk[0] & 0xF0) != 0x50) {
		return HAWK_ERR_FORMAT;
	}
	logn = sk[0] & 0x0F;
	if (logn < 1 || logn > 10 || seckey_len != HAWK_SECKEY_SIZE(logn)) {
		return HAWK_ERR_FORMAT;
	}
	if (pubkey_len < HAWK_PUBKEY_SIZE[logn]
		|| tmp_len < HAWK_TMPSIZE_MAKEPUB(logn)) {
		return HAWK_ERR_SIZE;
	}

	/*
	 * Decode secret key (f and g).
	 */
	n = MKN(logn);
	f = (int8_t *)tmp;
	g = f + n;
	F = g + n;

	if (Zf(decode_seckey)(f, g, F, seckey + 1, seckey_len - 1, logn) == 0) {
		return HAWK_ERR_FORMAT;
	}

	/*
	 * Compute public key.
	 */
	iq00 = align_i16(tmp);
	iq01 = iq00 + n;
	atmp = (uint8_t *)align_i32(iq01 + n);

	Zf(make_public)(f, g, F, NULL, iq00, iq01, NULL, logn, atmp);

	/*
	 * Encode public key.
	 */
	pk = pubkey;
	pk[0] = 0x00 + logn;

	pk_len = Zf(encode_pubkey)(pk + 1, HAWK_PUBKEY_SIZE[logn] - 1,
		iq00, iq01, logn);
	if (pk_len == 0) {
		return HAWK_ERR_FORMAT;
	}

	for (u = 1; u < pubkey_len; u ++) {
		/*
		 * Pad the public key with zeros in constant time, by doing:
		 * if (u >= 1 + pk_len) pk[u] = 0;
		 */
		pk[u] &= -(((uint8_t)(u - 1 - pk_len)) >> 7);
	}
	return 0;

}

/* see hawk.h */
int
hawk_get_logn(const void *obj, size_t len)
{
	unsigned logn;

	if (len == 0) {
		return HAWK_ERR_FORMAT;
	}
	logn = *(const uint8_t *)obj & 0x0F;
	if (logn < 1 || logn > 10) {
		return HAWK_ERR_FORMAT;
	}
	return logn;
}

/* see hawk.h */
int
hawk_expand_seckey(void *expanded_key, size_t expanded_key_len,
	const void *seckey, size_t seckey_len)
{
	unsigned logn, oldcw;
	const uint8_t *sk;
	int8_t *f, *g, *F;
	size_t n;
	fpr *expkey;

	/*
	 * Get degree from secret key header byte, and check parameters.
	 */
	if (seckey_len == 0) {
		return HAWK_ERR_FORMAT;
	}
	sk = seckey;
	if ((sk[0] & 0xF0) != 0x50) {
		return HAWK_ERR_FORMAT;
	}
	logn = sk[0] & 0x0F;
	if (logn < 1 || logn > 10 || seckey_len != HAWK_SECKEY_SIZE(logn)) {
		return HAWK_ERR_FORMAT;
	}
	if (expanded_key_len < HAWK_EXPANDEDKEY_SIZE(logn))
	{
		return HAWK_ERR_SIZE;
	}

	/*
	 * Expand secret key.
	 */
	*(uint8_t *)expanded_key = logn;
	expkey = align_fpr((uint8_t *)expanded_key + 1);

	/*
	 * Decode secret key elements, and complete secret key.
	 */
	n = MKN(logn);
	f = (int8_t *)(expkey + 3*n);
	g = f + n;
	F = g + n;
	if (Zf(decode_seckey)(f, g, F, sk + 1, seckey_len - 1, logn) == 0) {
		return HAWK_ERR_FORMAT;
	}

	oldcw = set_fpu_cw(2);
	Zf(expand_seckey)(expkey, f, g, F, logn);
	set_fpu_cw(oldcw);

	return 0;
}

/* ========================================================================= */

/* see hawk.h */
int
hawk_sign_dyn(shake256_context *rng, void *sig, size_t *sig_len, int sig_type,
	const void *seckey, size_t seckey_len, const void *data, size_t data_len,
	void *tmp, size_t tmp_len)
{
	shake256_context hd;
	uint8_t salt[40];

	hawk_sign_start(&hd);
	shake256_inject(&hd, data, data_len);
	return hawk_sign_dyn_finish(rng, sig, sig_len, sig_type,
		seckey, seckey_len, &hd, salt, tmp, tmp_len);
}

/* see hawk.h */
int
hawk_sign(shake256_context *rng, void *sig, size_t *sig_len, int sig_type,
	const void *expanded_key, const void *data, size_t data_len, void *tmp,
	size_t tmp_len)
{
	shake256_context hd;
	uint8_t salt[40];

	hawk_sign_start(&hd);
	shake256_inject(&hd, data, data_len);
	return hawk_sign_finish(rng, sig, sig_len, sig_type,
		expanded_key, &hd, salt, tmp, tmp_len);
}

/* ========================================================================= */

/* see hawk.h */
void
hawk_sign_start(shake256_context *hash_data)
{
	shake256_init(hash_data);
}

/* see hawk.h */
int
hawk_sign_dyn_finish(shake256_context *rng, void *sig, size_t *sig_len,
	int sig_type, const void *seckey, size_t seckey_len,
	shake256_context *hash_data, void *salt, void *tmp, size_t tmp_len)
{
#ifdef HAWK_AVX
	unsigned oldcw;
#endif
	unsigned logn;
	size_t n, u, v, es_len;
	int8_t *f, *g, *F;
	uint8_t header_byte, *es, *hm, *atmp;
	int16_t *sv;
	inner_shake256_context hash_state;

	/*
	 * Get degree from secret key header byte, and check parameters.
	 */
	if (seckey_len == 0) {
		return HAWK_ERR_FORMAT;
	}
	header_byte = ((uint8_t *)seckey)[0];
	if ((header_byte & 0xF0) != 0x50) {
		return HAWK_ERR_FORMAT;
	}
	logn = header_byte & 0x0F;
	if (logn < 1 || logn > 10 || seckey_len != HAWK_SECKEY_SIZE(logn)) {
		return HAWK_ERR_FORMAT;
	}
#ifdef HAWK_AVX
	if (tmp_len < HAWK_TMPSIZE_SIGNDYN(logn)) {
		return HAWK_ERR_SIZE;
	}
#else
	if (tmp_len < HAWK_TMPSIZE_SIGNDYN_NTT(logn)) {
		return HAWK_ERR_SIZE;
	}
#endif
	es_len = *sig_len;
	if (es_len < 1 + HAWK_SALT_SIZE(logn)) {
		return HAWK_ERR_SIZE;
	}

	switch (sig_type) {
	case HAWK_SIG_COMPACT:
		break;
	case HAWK_SIG_PADDED:
		if (es_len < HAWK_SIG_PADDED_SIZE(logn)) {
			return HAWK_ERR_SIZE;
		}
		es_len = HAWK_SIG_PADDED_SIZE(logn);
		break;
	default:
		return HAWK_ERR_BADARG;
	}

	/*
	 * Decode secret key elements, and complete secret key.
	 */
	n = MKN(logn);
	f = (int8_t *)tmp;
	g = f + n;
	F = g + n;
	hm = (uint8_t *)(F + n);
	sv = align_i16(hm + HAWK_HASH_SIZE(logn));
#ifdef HAWK_AVX
	atmp = (uint8_t *)align_fpr(sv + n);
#else
	atmp = (uint8_t *)align_i32(sv + n);
#endif

	if (Zf(decode_seckey)(f, g, F, seckey + 1, seckey_len - 1, logn) == 0) {
		return HAWK_ERR_FORMAT;
	}

	/*
	 * Fix the first byte (containing logn).
	 */
	es = sig;
	es[0] = 0x30 + logn;
	u = 1 + HAWK_SALT_SIZE(logn);

	for (;;) {
		/*
		 * Make a copy of the current state of hash_data, as we add a salt
		 * perhaps multiple times in case signing fails. Add the salt, flip the
		 * state and then hash salt + message to a point.
		 */
		memcpy(&hash_state, (inner_shake256_context *)hash_data,
			sizeof *hash_data);
		shake256_extract(rng, salt, HAWK_SALT_SIZE(logn));
		inner_shake256_inject(&hash_state, salt, HAWK_SALT_SIZE(logn));
		inner_shake256_flip(&hash_state);
		inner_shake256_extract(&hash_state, hm, HAWK_HASH_SIZE(logn));

		/*
		 * If signature generation fails, restart. This does not break
		 * constant-time discipline, since the norm-check fails with a small
		 * almost-constant probability and the recovery check works almost
		 * always.
		 */
#ifdef HAWK_AVX
		oldcw = set_fpu_cw(2);
		if (!Zf(sign_dyn)((inner_shake256_context *)rng, sv, f, g, F, NULL,
				hm, logn, atmp)) {
			set_fpu_cw(oldcw);
			continue;
		}
		set_fpu_cw(oldcw);
#else
		if (!Zf(sign_NTT)((inner_shake256_context *)rng, sv, f, g, F, NULL,
				hm, logn, atmp)) {
			continue;
		}
#endif

		/*
		 * Fix the bytes for the salt.
		 */
		memcpy(es + 1, salt, HAWK_SALT_SIZE(logn));

		v = Zf(encode_sig)(es + u, es_len - u, sv, logn);
		if (sig_type == HAWK_SIG_COMPACT) {
			if (v == 0) {
				return HAWK_ERR_SIZE;
			} else {
				*sig_len = u + v;
				return 0;
			}
		} else {
			/*
			 * sig_type == HAWK_SIG_PADDED
			 */
			if (v == 0) {
				continue;
			}
			if (u + v < es_len) {
				/*
				 * Pad with zeros
				 */
				memset(es + u + v, 0, es_len - (u + v));
			}
			*sig_len = es_len;
			return 0;
		}
	}
}

/* see hawk.h */
int
hawk_sign_finish(shake256_context *rng, void *sig, size_t *sig_len,
	int sig_type, const void *expanded_key, shake256_context *hash_data,
	void *salt, void *tmp, size_t tmp_len)
{
	unsigned logn, oldcw;
	uint8_t *es, *hm, *atmp;
	const fpr *expkey;
	int16_t *sv;
	size_t u, v, n, es_len;
	inner_shake256_context hash_state;

	/*
	 * Get degree from secret key header byte, and check parameters.
	 */
	logn = *(const uint8_t *)expanded_key;
	if (logn < 1 || logn > 10) {
		return HAWK_ERR_FORMAT;
	}
	if (tmp_len < HAWK_TMPSIZE_SIGN(logn)) {
		return HAWK_ERR_SIZE;
	}
	es_len = *sig_len;
	if (es_len < 1 + HAWK_SALT_SIZE(logn)) {
		return HAWK_ERR_SIZE;
	}

	switch (sig_type) {
	case HAWK_SIG_COMPACT:
		break;
	case HAWK_SIG_PADDED:
		if (es_len < HAWK_SIG_PADDED_SIZE(logn)) {
			return HAWK_ERR_SIZE;
		}
		es_len = HAWK_SIG_PADDED_SIZE(logn);
		break;
	default:
		return HAWK_ERR_BADARG;
	}

	n = MKN(logn);
	hm = (uint8_t *)tmp;
	sv = align_i16(hm + HAWK_HASH_SIZE(logn));
	atmp = (uint8_t *)align_fpr(sv + n);

	expkey = (const fpr *)align_fpr((uint8_t *)expanded_key + 1);

	/*
	 * Fix the first byte (containing logn).
	 */
	es = sig;
	es[0] = 0x30 + logn;
	u = 1 + HAWK_SALT_SIZE(logn);

	for (;;) {
		/*
		 * Make a copy of the current state of hash_data, as we add a salt
		 * perhaps multiple times in case signing fails. Add the salt, flip the
		 * state and then hash salt + message to a point.
		 */
		memcpy(&hash_state, (inner_shake256_context *)hash_data,
			sizeof *hash_data);
		shake256_extract(rng, salt, HAWK_SALT_SIZE(logn));
		inner_shake256_inject(&hash_state, salt, HAWK_SALT_SIZE(logn));
		inner_shake256_flip(&hash_state);
		inner_shake256_extract(&hash_state, hm, HAWK_HASH_SIZE(logn));

		/*
		 * If signature generation fails, restart. This does not break
		 * constant-time discipline, since the norm-check fails with a small
		 * almost-constant probability and the recovery check works almost
		 * always.
		 */
		oldcw = set_fpu_cw(2);
		if (!Zf(sign)((inner_shake256_context *)rng, sv, expkey, hm, logn,
				atmp)) {
			set_fpu_cw(oldcw);
			continue;
		}
		set_fpu_cw(oldcw);

		/*
		 * Fix the bytes for the salt.
		 */
		memcpy(es + 1, salt, HAWK_SALT_SIZE(logn));

		v = Zf(encode_sig)(es + u, es_len - u, sv, logn);
		if (sig_type == HAWK_SIG_COMPACT) {
			if (v == 0) {
				return HAWK_ERR_SIZE;
			} else {
				*sig_len = u + v;
				return 0;
			}
		} else {
			/*
			 * sig_type == HAWK_SIG_PADDED
			 */
			if (v == 0) {
				continue;
			}
			if (u + v < es_len) {
				/*
				 * Pad with zeros
				 */
				memset(es + u + v, 0, es_len - (u + v));
			}
			*sig_len = es_len;
			return 0;
		}
	}
}

/* ==========================================================================*/

/* see hawk.h */
int
hawk_verify(const void *sig, size_t sig_len, int sig_type, const void *pubkey,
	size_t pubkey_len, const void *data, size_t data_len, void *tmp,
	size_t tmp_len)
{
	shake256_context hd;
	hawk_verify_start(&hd);
	shake256_inject(&hd, data, data_len);
	return hawk_verify_finish(sig, sig_len, sig_type,
		pubkey, pubkey_len, &hd, tmp, tmp_len);
}

/* ========================================================================= */

/* see hawk.h */
void
hawk_verify_start(shake256_context *hash_data)
{
	shake256_init(hash_data);
}

/* see hawk.h */
int
hawk_verify_finish(const void *sig, size_t sig_len, int sig_type,
	const void *pubkey, size_t pubkey_len, shake256_context *hash_data,
	void *tmp, size_t tmp_len)
{
	unsigned logn, salt_len;
	uint8_t *hm, *atmp;
	int16_t *iq00, *iq01;
#ifdef HAWK_AVX
	fpr *q00, *q01, *q11;
#endif
	const uint8_t *pk, *es;
	size_t u, v, n;
	int16_t *sv;

	/*
	 * Get Hawk degree from public key; verify consistency with
	 * signature value, and check parameters.
	 */
	if (sig_len == 0 || pubkey_len == 0) {
		return HAWK_ERR_FORMAT;
	}
	es = sig;
	pk = pubkey;
	if ((pk[0] & 0xF0) != 0x00) {
		return HAWK_ERR_FORMAT;
	}
	logn = pk[0] & 0x0F;
	if (logn < 1 || logn > 10) {
		return HAWK_ERR_FORMAT;
	}
	if ((es[0] & 0x0F) != logn || (es[0] & 0xF0) != 0x30) {
		return HAWK_ERR_BADSIG;
	}

	salt_len = HAWK_SALT_SIZE(logn);
	if (sig_len < 1 + salt_len) {
		return HAWK_ERR_FORMAT;
	}

	switch (sig_type) {
	case 0:
	case HAWK_SIG_COMPACT:
		break;
	case HAWK_SIG_PADDED:
		if (sig_len != HAWK_SIG_PADDED_SIZE(logn)) {
			return HAWK_ERR_FORMAT;
		}
		break;
	default:
		return HAWK_ERR_BADARG;
	}

	if (pubkey_len < HAWK_PUBKEY_SIZE[logn]) {
		return HAWK_ERR_FORMAT;
	}
#ifdef HAWK_AVX
	if (tmp_len < HAWK_TMPSIZE_VERIFY(logn)) {
		return HAWK_ERR_SIZE;
	}
#else
	if (tmp_len < HAWK_TMPSIZE_VERIFY_NTT(logn)) {
		return HAWK_ERR_SIZE;
	}
#endif

	n = MKN(logn);
	hm = (uint8_t *)tmp;
	sv = align_i16(hm + HAWK_HASH_SIZE(logn));
	iq00 = sv + n;
	iq01 = iq00 + n;

	/*
	 * Decode public key.
	 */
	if (Zf(decode_pubkey)(iq00, iq01, pk + 1, pubkey_len - 1, logn) == 0)
	{
		return HAWK_ERR_FORMAT;
	}

	/*
	 * Decode signature value.
	 */
	u = 1 + salt_len;
	v = Zf(decode_sig)(sv, es + u, sig_len - u, logn);
	if (v == 0) {
		return HAWK_ERR_FORMAT;
	}
	u += v;
	if (u != sig_len) {
		/*
		 * Extra bytes of value 0 are tolerated only for the "padded" format.
		 */
		if (sig_type == HAWK_SIG_PADDED)
		{
			while (u < sig_len) {
				if (es[u] != 0) {
					return HAWK_ERR_FORMAT;
				}
				u ++;
			}
		} else {
			return HAWK_ERR_FORMAT;
		}
	}

	/*
	 * Add the salt, flip the state and then hash salt + message to a point.
	 */
	shake256_inject(hash_data, es + 1, salt_len);
	shake256_flip(hash_data);
	inner_shake256_extract((inner_shake256_context *)hash_data, hm,
		HAWK_HASH_SIZE(logn));

#ifdef HAWK_AVX
	q00 = (fpr *)align_fpr(sv + n);
	q01 = q00 + n;
	q11 = q01 + n;
	atmp = (uint8_t *)(q11 + n);

	/*
	 * Construct full public key and verify signature.
	 */
	Zf(complete_pubkey)(iq00, iq01, q00, q01, q11, logn);
	if (!Zf(verify)(hm, sv, q00, q01, q11, logn, atmp)) {
		return HAWK_ERR_BADSIG;
	}
#else
	atmp = (uint8_t *)align_fpr(iq01 + n);

	/*
	 * Verify signature.
	 */
	if (!Zf(verify_NTT)(hm, sv, iq00, iq01, logn, atmp)) {
		return HAWK_ERR_BADSIG;
	}
#endif
	return 0;
}
