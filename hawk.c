/*
 * Implementation of the external Hawk API.
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

/*
 * These are the sizes that are allowed for the secret key and private key for
 * various values of logn.
 */
const size_t HAWK_MIN_SECKEY_SIZE[10] = {
	0 /* unused */, 1, 3, 9, 23, 46, 93, 179, 339, 646
};

const size_t HAWK_MAX_SECKEY_SIZE[10] = {
	0 /* unused */, 6, 11, 22, 37, 70, 120, 228, 395, 712
};

const size_t HAWK_MIN_PUBKEY_SIZE[10] = {
	0 /* unused */, 3, 7, 13, 28, 57, 117, 235, 475, 943
};

const size_t HAWK_MAX_PUBKEY_SIZE[10] = {
	0 /* unused */, 7, 13, 23, 41, 77, 143, 276, 528, 1026
};

// TODO: set 5 in the compress.c code (depending on logn!).
#define SIG_LOBITS(logn) (5)

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

/* see hawk.h */
int
hawk_keygen_make(shake256_context *rng, unsigned logn, void *privkey,
	size_t *privkey_len, void *pubkey, size_t *pubkey_len, void *tmp,
	size_t tmp_len)
{
	int8_t *f, *g, *F, *G;
	int16_t *q00_num, *q10_num;
	fpr *q00, *q10, *q11;
	uint8_t *atmp;
	size_t n, u, sk_len, pk_len;
	uint8_t *sk, *pk;
	unsigned oldcw;

	/*
	 * Check parameters.
	 */
	if (logn < 1 || logn > 9) {
		return HAWK_ERR_BADARG;
	}

	/*
	 * Check that the seckey and pubkey buffers are at least as large as the
	 * allowed encoded sizes, and check the temporary buffer size.
	 */
	if (*privkey_len < HAWK_MAX_SECKEY_SIZE[logn]
		|| (pubkey != NULL && *pubkey_len < HAWK_MAX_PUBKEY_SIZE[logn])
		|| tmp_len < HAWK_TMPSIZE_KEYGEN(logn)) {
		return HAWK_ERR_SIZE;
	}

	/*
	 * Prepare buffers and generate private key.
	 */
	n = MKN(logn);
	f = (int8_t *)tmp;
	g = f + n;
	F = g + n;
	G = F + n;
	q00 = align_fpr(G + n);
	q10 = q00 + n;
	q11 = q10 + n;

	/*
	 * These buffers overlap with f, g, F, G but this is fine as we encode the
	 * public key after the private key is already encoded. Also, the first
	 * byte of q00 may overlap with q10_num because of alignment, but the data
	 * of q00[0] is first extracted and only afterwards is the last value of
	 * q10_num set.
	 */
	q00_num = align_i16(tmp);
	q10_num = q00_num + n;

	atmp = (uint8_t *)(q11 + n);

	/*
	 * Fix the first byte of secret key and private key.
	 */
	sk = privkey;
	sk[0] = 0x50 + logn;

	pk = pubkey;
	pk_len = 0;
	if (pubkey != NULL) {
		pk[0] = 0x00 + logn;
	}

	do {
		oldcw = set_fpu_cw(2);
		Zf(keygen)((inner_shake256_context *)rng, f, g, F, G, q00, q10, q11,
			logn, atmp);
		set_fpu_cw(oldcw);

		sk_len = Zf(encode_seckey)(sk + 1, *privkey_len - 1, f, g, F, logn);

		/*
		 * Destroy the private key basis [[f,g], [F,G]] to store q00, q10.
		 */
		Zf(iFFT)(q00, logn);
		Zf(iFFT)(q10, logn);

		for (u = 0; u < n; u++) {
			q00_num[u] = fpr_rint(q00[u]);
		}
		for (u = 0; u < n; u++) {
			q10_num[u] = fpr_rint(q10[u]);
		}

		if (pubkey != NULL) {
			pk_len = Zf(encode_pubkey)(pk + 1, *pubkey_len - 1, q00_num, q10_num, logn);
		}
		/*
		 * Retry key-generation as the secret key or public key cannot be
		 * encoded. This only happens with negligible probability.
		 */
	} while (sk_len == 0 || pk_len == 0);

	/*
	 * Do not forgot that there is one header byte in sk and pk.
	 */
	*privkey_len = 1 + sk_len;
	if (pubkey != NULL) {
		*pubkey_len = 1 + pk_len;
	}

	return 0;
}

/* see hawk.h */
int
hawk_make_public(void *pubkey, size_t *pubkey_len, const void *privkey,
	size_t privkey_len, void *tmp, size_t tmp_len)
{
	uint8_t *pk, *atmp;
	const uint8_t *sk;
	unsigned logn;
	size_t v, n;
	int8_t *f, *g, *F;
	int16_t *q00, *q10;
	fpr *bq00, *bq10;

	/*
	 * Get degree from private key header byte, and check
	 * parameters.
	 */
	if (privkey_len == 0) {
		return HAWK_ERR_FORMAT;
	}
	sk = privkey;
	if ((sk[0] & 0xF0) != 0x50) {
		return HAWK_ERR_FORMAT;
	}
	logn = sk[0] & 0x0F;
	if (logn < 1 || logn > 9) {
		return HAWK_ERR_FORMAT;
	}
	if (privkey_len < HAWK_MIN_SECKEY_SIZE[logn]
		|| privkey_len > HAWK_MAX_SECKEY_SIZE[logn]) {
		return HAWK_ERR_FORMAT;
	}
	if (*pubkey_len < HAWK_MAX_PUBKEY_SIZE[logn]
		|| tmp_len < HAWK_TMPSIZE_MAKEPUB(logn))
	{
		return HAWK_ERR_SIZE;
	}

	/*
	 * Decode private key (f and g).
	 */
	n = MKN(logn);
	f = (int8_t *)tmp;
	g = f + n;
	F = g + n;

	if (Zf(decode_seckey)(f, g, F, privkey + 1, privkey_len - 1, logn) == 0) {
		return HAWK_ERR_FORMAT;
	}

	/*
	 * Compute public key.
	 */
	q00 = align_i16(tmp);
	q10 = q00 + n;
	bq00 = align_fpr(q10 + n);
	bq10 = bq00 + n;
	atmp = (uint8_t *)(bq10 + n);

	Zf(make_public)(f, g, F, NULL, bq00, bq10, NULL, logn, atmp);

	Zf(fft_to_int16)(q10, bq10, logn);
	Zf(fft_to_int16)(q00, bq00, logn);

	/*
	 * Encode public key.
	 */
	pk = pubkey;
	pk[0] = 0x00 + logn;
	v = Zf(encode_pubkey)(pk + 1, *pubkey_len - 1, q00, q10, logn);
	if (v == 0) {
		return HAWK_ERR_FORMAT;
	}
	*pubkey_len = 1 + v;
	return 0;

}

/* see hawk.h */
int
hawk_get_logn(void *obj, size_t len)
{
	int logn;

	if (len == 0) {
		return HAWK_ERR_FORMAT;
	}
	logn = *(uint8_t *)obj & 0x0F;
	if (logn < 1 || logn > 9) {
		return HAWK_ERR_FORMAT;
	}
	return logn;
}

/* see hawk.h */
int
hawk_sign_start(shake256_context *rng, void *salt, shake256_context *hash_data)
{
	shake256_extract(rng, salt, 40);
	shake256_init(hash_data);
	shake256_inject(hash_data, salt, 40);
	return 0;
}

/* see hawk.h */
int
hawk_expand_privkey(void *expanded_key, size_t expanded_key_len,
	const void *privkey, size_t privkey_len, void *tmp, size_t tmp_len)
{
	unsigned logn;
	const uint8_t *sk;
	int8_t *f, *g, *F;
	size_t n, v;
	fpr *expkey;
	unsigned oldcw;

	/*
	 * Get degree from private key header byte, and check
	 * parameters.
	 */
	if (privkey_len == 0) {
		return HAWK_ERR_FORMAT;
	}
	sk = privkey;
	if ((sk[0] & 0xF0) != 0x50) {
		return HAWK_ERR_FORMAT;
	}
	logn = sk[0] & 0x0F;
	if (logn < 1 || logn > 9) {
		return HAWK_ERR_FORMAT;
	}
	if (privkey_len < HAWK_MIN_SECKEY_SIZE[logn]
		|| privkey_len > HAWK_MAX_SECKEY_SIZE[logn]) {
		return HAWK_ERR_FORMAT;
	}
	if (expanded_key_len < HAWK_EXPANDEDKEY_SIZE(logn)
		|| tmp_len < HAWK_TMPSIZE_EXPANDPRIV(logn))
	{
		return HAWK_ERR_SIZE;
	}

	/*
	 * Decode private key elements, and complete private key.
	 */
	n = MKN(logn);
	f = (int8_t *)tmp;
	g = f + n;
	F = g + n;

	v = Zf(decode_seckey)(f, g, F, sk + 1, privkey_len - 1, logn);
	if (v == 0) {
		return HAWK_ERR_FORMAT;
	}

	/*
	 * Expand private key.
	 */
	*(uint8_t *)expanded_key = logn;
	expkey = align_fpr((uint8_t *)expanded_key + 1);

	oldcw = set_fpu_cw(2);
	Zf(expand_seckey)(expkey, f, g, F, logn);
	set_fpu_cw(oldcw);

	return 0;
}

/* see hawk.h */
int
hawk_sign_finish(shake256_context *rng, void *sig, size_t *sig_len,
	int sig_type, const void *expanded_key, shake256_context *hash_data,
	const void *salt, void *tmp, size_t tmp_len)
{
	unsigned logn;
	uint8_t *es, *hm, *atmp;
	const fpr *expkey;
	int16_t *sv;
	size_t u, v, n, es_len;
	unsigned oldcw;

	/*
	 * Get degree from private key header byte, and check
	 * parameters.
	 */
	logn = *(const uint8_t *)expanded_key;
	if (logn < 1 || logn > 9) {
		return HAWK_ERR_FORMAT;
	}
	if (tmp_len < HAWK_TMPSIZE_SIGN(logn)) {
		return HAWK_ERR_SIZE;
	}
	es_len = *sig_len;
	if (es_len < 41) {
		return HAWK_ERR_SIZE;
	}
	expkey = (const fpr *)align_fpr((uint8_t *)expanded_key + 1);
	switch (sig_type) {
	case HAWK_SIG_COMPRESSED:
		break;
	case HAWK_SIG_PADDED:
		if (es_len < HAWK_SIG_PADDED_SIZE(logn)) {
			return HAWK_ERR_SIZE;
		}
		break;
	default:
		return HAWK_ERR_BADARG;
	}

	n = MKN(logn);
	hm = (uint8_t *)tmp;
	sv = align_i16(hm + HAWK_HASH_SIZE(logn));
	atmp = (uint8_t *)align_fpr(sv + n);

	/*
	 * Hash message to a point.
	 */
	shake256_flip(hash_data);
	inner_shake256_extract((inner_shake256_context *)hash_data, hm,
		HAWK_HASH_SIZE(logn));

	/*
	 * Fix the first byte (containing logn) and 40 bytes for the salt first.
	 */
	es = sig;
	memcpy(es + 1, salt, 40);
	u = 41;
	es[0] = 0x30 + logn;

	if (sig_type == HAWK_SIG_COMPRESSED) {
		oldcw = set_fpu_cw(2);
		Zf(sign)((inner_shake256_context *)rng, sv, expkey, hm, logn, atmp);
		set_fpu_cw(oldcw);

		v = Zf(encode_sig)(es + u, es_len - u, sv, logn, SIG_LOBITS(logn));
		if (v == 0) {
			return HAWK_ERR_SIZE;
		} else {
			*sig_len = u + v;
			return 0;
		}
	}

	/*
	 * Now, sig_type is HAWK_SIG_PADDED.
	 * Compute the signature until one is found that is encodable.
	 */
	es_len = HAWK_SIG_PADDED_SIZE(logn);

	do {
		oldcw = set_fpu_cw(2);
		Zf(sign)((inner_shake256_context *)rng, sv, expkey, hm, logn, atmp);
		set_fpu_cw(oldcw);

		v = Zf(encode_sig)(es + u, es_len - u, sv, logn, SIG_LOBITS(logn));
		/*
		 * If v = 0, the signature does not fit and loop.
		 */
	} while (v == 0);

	if (u + v < es_len) {
		/*
		 * Pad with zeros
		 */
		memset(es + u + v, 0, es_len - (u + v));
	}
	*sig_len = es_len;
	return 0;
}

/* see hawk.h */
int
hawk_sign(shake256_context *rng, void *sig, size_t *sig_len, int sig_type,
	const void *expanded_key, const void *data, size_t data_len, void *tmp,
	size_t tmp_len)
{
	shake256_context hd;
	uint8_t salt[40];
	int r;

	r = hawk_sign_start(rng, salt, &hd);
	if (r != 0) {
		return r;
	}
	shake256_inject(&hd, data, data_len);
	return hawk_sign_finish(rng, sig, sig_len, sig_type, expanded_key, &hd,
		salt, tmp, tmp_len);
}

/* see hawk.h */
int
hawk_sign_dyn_finish(shake256_context *rng, void *sig, size_t *sig_len,
	int sig_type, const void *privkey, size_t privkey_len,
	shake256_context *hash_data, const void *salt, void *tmp, size_t tmp_len)
{
	unsigned logn, oldcw;
	size_t n, u, v, es_len;
	int8_t *f, *g, *F;
	uint8_t header_byte, *es, *hm, *atmp;
	int16_t *sv;

	/*
	 * Get degree from private key header byte, and check
	 * parameters.
	 */
	if (privkey_len == 0) {
		return HAWK_ERR_FORMAT;
	}
	header_byte = ((uint8_t *)privkey)[0];
	if ((header_byte & 0xF0) != 0x50) {
		return HAWK_ERR_FORMAT;
	}
	logn = header_byte & 0x0F;
	if (logn < 1 || logn > 9) {
		return HAWK_ERR_FORMAT;
	}
	if (privkey_len < HAWK_MIN_SECKEY_SIZE[logn]
		|| privkey_len > HAWK_MAX_SECKEY_SIZE[logn]) {
		return HAWK_ERR_FORMAT;
	}
	if (tmp_len < HAWK_TMPSIZE_SIGNDYN(logn)) {
		return HAWK_ERR_SIZE;
	}
	es_len = *sig_len;
	if (es_len < 41) {
		return HAWK_ERR_SIZE;
	}

	switch (sig_type) {
	case HAWK_SIG_COMPRESSED:
		break;
	case HAWK_SIG_PADDED:
		if (es_len < HAWK_SIG_PADDED_SIZE(logn)) {
			return HAWK_ERR_SIZE;
		}
		break;
	default:
		return HAWK_ERR_BADARG;
	}

	/*
	 * Decode private key elements, and complete private key.
	 */
	n = MKN(logn);
	f = (int8_t *)tmp;
	g = f + n;
	F = g + n;
	hm = (uint8_t *)(F + n);
	sv = align_i16(hm + HAWK_HASH_SIZE(logn));
	atmp = (uint8_t *)align_fpr(sv + n);

	v = Zf(decode_seckey)(f, g, F, privkey + 1, privkey_len - 1, logn);
	if (v == 0) {
		return HAWK_ERR_FORMAT;
	}

	/*
	 * Hash message to a point.
	 */
	shake256_flip(hash_data);
	inner_shake256_extract((inner_shake256_context *)hash_data, hm,
		HAWK_HASH_SIZE(logn));

	/*
	 * Fix the first byte (containing logn) and 40 bytes for the salt first.
	 */
	es = sig;
	memcpy(es + 1, salt, 40);
	u = 41;
	es[0] = 0x30 + logn;

	if (sig_type == HAWK_SIG_COMPRESSED) {
		oldcw = set_fpu_cw(2);
		Zf(sign_dyn)((inner_shake256_context *)rng, sv, f, g, F, NULL, hm,
			logn, atmp);
		set_fpu_cw(oldcw);

		v = Zf(encode_sig)(es + u, es_len - u, sv, logn, SIG_LOBITS(logn));
		if (v == 0) {
			return HAWK_ERR_SIZE;
		} else {
			*sig_len = u + v;
			return 0;
		}
	}

	/*
	 * Now, sig_type is HAWK_SIG_PADDED.
	 * Compute the signature until one is found that is encodable.
	 */
	es_len = HAWK_SIG_PADDED_SIZE(logn);

	do {
		oldcw = set_fpu_cw(2);
		Zf(sign_dyn)((inner_shake256_context *)rng, sv, f, g, F, NULL, hm,
			logn, atmp);
		set_fpu_cw(oldcw);

		v = Zf(encode_sig)(es + u, es_len - u, sv, logn, SIG_LOBITS(logn));
		/*
		 * If v = 0, the signature does not fit and loop.
		 */
	} while (v == 0);

	if (u + v < es_len) {
		/*
		 * Pad with zeros
		 */
		memset(es + u + v, 0, es_len - (u + v));
	}
	*sig_len = es_len;
	return 0;

}

/* see hawk.h */
int
hawk_sign_dyn(shake256_context *rng, void *sig, size_t *sig_len, int sig_type,
	const void *privkey, size_t privkey_len, const void *data, size_t data_len,
	void *tmp, size_t tmp_len)
{
	shake256_context hd;
	uint8_t salt[40];
	int r;

	r = hawk_sign_start(rng, salt, &hd);
	if (r != 0) {
		return r;
	}
	shake256_inject(&hd, data, data_len);
	return hawk_sign_dyn_finish(rng, sig, sig_len, sig_type,
		privkey, privkey_len, &hd, salt, tmp, tmp_len);
}


/* see hawk.h */
int
hawk_verify_start(shake256_context *hash_data, const void *sig, size_t sig_len)
{
	if (sig_len < 41) {
		return HAWK_ERR_FORMAT;
	}
	shake256_init(hash_data);
	/*
	 * First inject the salt in the SHAKE context.
	 */
	shake256_inject(hash_data, (const uint8_t *)sig + 1, 40);
	return 0;
}

/* see hawk.h */
int
hawk_verify_finish(const void *sig, size_t sig_len, int sig_type,
	const void *pubkey, size_t pubkey_len, shake256_context *hash_data,
	void *tmp, size_t tmp_len)
{
	unsigned logn;
	uint8_t *hm, *atmp;
	int16_t *q00_num, *q10_num;
	fpr *q00, *q10, *q11;
	const uint8_t *pk, *es;
	size_t u, v, n;
	int16_t *sv;

	/*
	 * Get Hawk degree from public key; verify consistency with
	 * signature value, and check parameters.
	 */
	if (sig_len < 41 || pubkey_len == 0) {
		return HAWK_ERR_FORMAT;
	}
	es = sig;
	pk = pubkey;
	if ((pk[0] & 0xF0) != 0x00) {
		return HAWK_ERR_FORMAT;
	}
	logn = pk[0] & 0x0F;
	if (logn < 1 || logn > 9) {
		return HAWK_ERR_FORMAT;
	}
	if ((es[0] & 0x0F) != logn) {
		return HAWK_ERR_BADSIG;
	}

	if ((es[0] & 0xF0) != 0x30) {
		return HAWK_ERR_BADSIG;
	}

	switch (sig_type) {
	case 0:
	case HAWK_SIG_COMPRESSED:
		break;
	case HAWK_SIG_PADDED:
		if (sig_len != HAWK_SIG_PADDED_SIZE(logn)) {
			return HAWK_ERR_FORMAT;
		}
		break;
	default:
		return HAWK_ERR_BADARG;
	}

	if (pubkey_len < HAWK_MIN_PUBKEY_SIZE[logn]
		|| pubkey_len > HAWK_MAX_PUBKEY_SIZE[logn]) {
		return HAWK_ERR_FORMAT;
	}
	if (tmp_len < HAWK_TMPSIZE_VERIFY(logn)) {
		return HAWK_ERR_SIZE;
	}

	n = MKN(logn);
	hm = (uint8_t *)tmp;
	sv = align_i16(hm + HAWK_HASH_SIZE(logn));

	q00_num = sv + n;
	q10_num = q00_num + n;

	q00 = (fpr *)align_fpr(sv + n);
	q10 = q00 + n;
	q11 = q10 + n;
	atmp = (uint8_t *)(q11 + n);

	/*
	 * Decode public key.
	 */
	if (Zf(decode_pubkey)(q00_num, q10_num, pk + 1, pubkey_len - 1, logn) == 0)
	{
		return HAWK_ERR_FORMAT;
	}

	/*
	 * Decode signature value.
	 */
	u = 41;
	v = Zf(decode_sig)(sv, es + u, sig_len - u, logn, SIG_LOBITS(logn));
	if (v == 0) {
		return HAWK_ERR_FORMAT;
	}
	u += v;
	if (u != sig_len) {
		/*
		 * Extra bytes of value 0 are tolerated only for the
		 * "padded" format.
		 */
		if ((sig_type == 0 && sig_len == HAWK_SIG_PADDED_SIZE(logn))
			|| sig_type == HAWK_SIG_PADDED)
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
	 * Hash message to point.
	 */
	shake256_flip(hash_data);
	inner_shake256_extract((inner_shake256_context *)hash_data, hm,
		HAWK_HASH_SIZE(logn));

	/*
	 * Construct full public key.
	 */
	Zf(complete_pubkey)(q00_num, q10_num, q00, q10, q11, logn);

	/*
	 * Verify signature.
	 */
	if (!Zf(verify_simple_rounding_fft)(hm, sv, q00, q10, q11, logn, atmp)) {
		return HAWK_ERR_BADSIG;
	}
	return 0;
}

/* see hawk.h */
int
hawk_verify(const void *sig, size_t sig_len, int sig_type, const void *pubkey,
	size_t pubkey_len, const void *data, size_t data_len, void *tmp,
	size_t tmp_len)
{
	shake256_context hd;
	int r;

	r = hawk_verify_start(&hd, sig, sig_len);
	if (r < 0) {
		return r;
	}
	shake256_inject(&hd, data, data_len);
	return hawk_verify_finish(sig, sig_len, sig_type,
		pubkey, pubkey_len, &hd, tmp, tmp_len);
}
