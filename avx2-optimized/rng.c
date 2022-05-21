/*
 * PRNG and interface to the system RNG.
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
 * Include relevant system header files. For Win32, this will also need
 * linking with advapi32.dll, which we trigger with an appropriate #pragma.
 */
#if HAWK_RAND_GETENTROPY
#include <unistd.h>
#endif
#if HAWK_RAND_URANDOM
#include <sys/types.h>
#if !HAWK_RAND_GETENTROPY
#include <unistd.h>
#endif
#include <fcntl.h>
#include <errno.h>
#endif
#if HAWK_RAND_WIN32
#include <windows.h>
#include <wincrypt.h>
#pragma comment(lib, "advapi32")
#endif

/* see inner.h */
int
Zf(get_seed)(void *seed, size_t len)
{
	(void)seed;
	if (len == 0) {
		return 1;
	}
#if HAWK_RAND_GETENTROPY
	if (getentropy(seed, len) == 0) {
		return 1;
	}
#endif
#if HAWK_RAND_URANDOM
	{
		int f;

		f = open("/dev/urandom", O_RDONLY);
		if (f >= 0) {
			while (len > 0) {
				ssize_t rlen;

				rlen = read(f, seed, len);
				if (rlen < 0) {
					if (errno == EINTR) {
						continue;
					}
					break;
				}
				seed = (uint8_t *)seed + rlen;
				len -= (size_t)rlen;
			}
			close(f);
			if (len == 0) {
				return 1;
			}
		}
	}
#endif
#if HAWK_RAND_WIN32
	{
		HCRYPTPROV hp;

		if (CryptAcquireContext(&hp, 0, 0, PROV_RSA_FULL,
			CRYPT_VERIFYCONTEXT | CRYPT_SILENT))
		{
			BOOL r;

			r = CryptGenRandom(hp, (DWORD)len, seed);
			CryptReleaseContext(hp, 0);
			if (r) {
				return 1;
			}
		}
	}
#endif
	return 0;
}


/* see inner.h */
void
Zf(prng_init)(prng *p, inner_shake256_context *src)
{
	inner_shake256_extract(src, p->state.d, 56);
	Zf(prng_refill)(p);
}

/*
 * PRNG based on ChaCha20.
 *
 * State consists in key (32 bytes) then IV (16 bytes) and block counter
 * (8 bytes). Normally, we should not care about local endianness (this
 * is for a PRNG), but for the NIST competition we need reproducible KAT
 * vectors that work across architectures, so we enforce little-endian
 * interpretation where applicable. Moreover, output words are "spread
 * out" over the output buffer with the interleaving pattern that is
 * naturally obtained from the AVX2 implementation that runs eight
 * ChaCha20 instances in parallel.
 *
 * The block counter is XORed into the first 8 bytes of the IV.
 */
TARGET_AVX2
void
Zf(prng_refill)(prng *p)
{

	static const uint32_t CW[] = {
		0x61707865, 0x3320646e, 0x79622d32, 0x6b206574
	};

	uint64_t cc;
	size_t u;
	int i;
	uint32_t *sw;
	union {
		uint32_t w[16];
		__m256i y[2];  /* for alignment */
	} t;
	__m256i state[16], init[16];

	sw = (uint32_t *)p->state.d;

	/*
	 * XOR next counter values into state.
	 */
	cc = *(uint64_t *)(p->state.d + 48);
	for (u = 0; u < 8; u ++) {
		t.w[u] = (uint32_t)(cc + u);
		t.w[u + 8] = (uint32_t)((cc + u) >> 32);
	}
	*(uint64_t *)(p->state.d + 48) = cc + 8;

	/*
	 * Load state.
	 */
	for (u = 0; u < 4; u ++) {
		state[u] = init[u] =
			_mm256_broadcastd_epi32(_mm_cvtsi32_si128(CW[u]));
	}
	for (u = 0; u < 10; u ++) {
		state[u + 4] = init[u + 4] =
			_mm256_broadcastd_epi32(_mm_cvtsi32_si128(sw[u]));
	}
	state[14] = init[14] = _mm256_xor_si256(
		_mm256_broadcastd_epi32(_mm_cvtsi32_si128(sw[10])),
		_mm256_loadu_si256((__m256i *)&t.w[0]));
	state[15] = init[15] = _mm256_xor_si256(
		_mm256_broadcastd_epi32(_mm_cvtsi32_si128(sw[11])),
		_mm256_loadu_si256((__m256i *)&t.w[8]));

	/*
	 * Do all rounds.
	 */
	for (i = 0; i < 10; i ++) {

#define QROUND(a, b, c, d)   do { \
		state[a] = _mm256_add_epi32(state[a], state[b]); \
		state[d] = _mm256_xor_si256(state[d], state[a]); \
		state[d] = _mm256_or_si256( \
			_mm256_slli_epi32(state[d], 16), \
			_mm256_srli_epi32(state[d], 16)); \
		state[c] = _mm256_add_epi32(state[c], state[d]); \
		state[b] = _mm256_xor_si256(state[b], state[c]); \
		state[b] = _mm256_or_si256( \
			_mm256_slli_epi32(state[b], 12), \
			_mm256_srli_epi32(state[b], 20)); \
		state[a] = _mm256_add_epi32(state[a], state[b]); \
		state[d] = _mm256_xor_si256(state[d], state[a]); \
		state[d] = _mm256_or_si256( \
			_mm256_slli_epi32(state[d],  8), \
			_mm256_srli_epi32(state[d], 24)); \
		state[c] = _mm256_add_epi32(state[c], state[d]); \
		state[b] = _mm256_xor_si256(state[b], state[c]); \
		state[b] = _mm256_or_si256( \
			_mm256_slli_epi32(state[b], 7), \
			_mm256_srli_epi32(state[b], 25)); \
	} while (0)

		QROUND( 0,  4,  8, 12);
		QROUND( 1,  5,  9, 13);
		QROUND( 2,  6, 10, 14);
		QROUND( 3,  7, 11, 15);
		QROUND( 0,  5, 10, 15);
		QROUND( 1,  6, 11, 12);
		QROUND( 2,  7,  8, 13);
		QROUND( 3,  4,  9, 14);

#undef QROUND

	}

	/*
	 * Add initial state back and encode the result in the destination
	 * buffer. We can dump the AVX2 values "as is" because the non-AVX2
	 * code uses a compatible order of values.
	 */
	for (u = 0; u < 16; u ++) {
		_mm256_storeu_si256((__m256i *)&p->buf.d[u << 5],
			_mm256_add_epi32(state[u], init[u]));
	}


	p->ptr = 0;
}

/* see inner.h */
void
Zf(prng_get_bytes)(prng *p, void *dst, size_t len)
{
	uint8_t *buf;

	buf = (uint8_t *)dst;
	while (len > 0) {
		size_t clen;

		clen = (sizeof p->buf.d) - p->ptr;
		if (clen > len) {
			clen = len;
		}
		memcpy(buf, p->buf.d, clen);
		buf += clen;
		len -= clen;
		p->ptr += clen;
		if (p->ptr == sizeof p->buf.d) {
			Zf(prng_refill)(p);
		}
	}
}

/* see inner.h */
void
Zf(prng_get_80_bits)(prng *p, uint16_t *bit16, uint64_t *bit64)
{
	size_t u;

	/*
	 * This function requires 10 bytes to be present in the buffer.
	 */
	u = p->ptr;
	if (u + 10 >= sizeof p->buf.d) {
		/*
		 * If there are less than 10 bytes in the buffer, first refill it and
		 * then read the bytes from the freshly filled buffer.
		 */
		Zf(prng_refill)(p);
		u = 0;
	}

	/*
	 * On systems that use little-endian encoding and allow unaligned accesses,
	 * simply read the data where it is.
	 */
	*bit16 = (uint16_t)p->buf.d[u+8] | ((uint16_t)p->buf.d[u+9] << 8);
	*bit64 = *(uint64_t *)&p->buf.d[u];
	p->ptr = u + 10;
}
