/*
 * Discrete gaussian sampling
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
 * Sample an integer value along a half-gaussian distribution centered
 * on zero and standard deviation 1.8205, with a precision of 72 bits.
 */
TARGET_AVX2
int
Zf(gaussian0_sampler)(prng *p)
{

	/*
	 * High words.
	 */
	static const union {
		uint16_t u16[16];
		__m256i ymm[1];
	} rhi15 = {
		{
			0x51FB, 0x2A69, 0x113E, 0x0568,
			0x014A, 0x003B, 0x0008, 0x0000,
			0x0000, 0x0000, 0x0000, 0x0000,
			0x0000, 0x0000, 0x0000, 0x0000
		}
	};

	static const union {
		uint64_t u64[20];
		__m256i ymm[5];
	} rlo57 = {
		{
			0x1F42ED3AC391802, 0x12B181F3F7DDB82,
			0x1CDD0934829C1FF, 0x1754377C7994AE4,
			0x1846CAEF33F1F6F, 0x14AC754ED74BD5F,
			0x024DD542B776AE4, 0x1A1FFDC65AD63DA,
			0x01F80D88A7B6428, 0x001C3FDB2040C69,
			0x00012CF24D031FB, 0x00000949F8B091F,
			0x0000003665DA998, 0x00000000EBF6EBB,
			0x0000000002F5D7E, 0x000000000007098,
			0x0000000000000C6, 0x000000000000001,
			0x000000000000000, 0x000000000000000
		}
	};

	uint64_t lo;
	unsigned hi;
	__m256i xhi, rhi, gthi, eqhi, eqm;
	__m256i xlo, gtlo0, gtlo1, gtlo2, gtlo3, gtlo4;
	__m128i t, zt;
	int r;

	/*
	 * Get a 72-bit random value and split it into a low part
	 * (57 bits) and a high part (15 bits)
	 */
	lo = prng_get_u64(p);
	hi = prng_get_u8(p);
	hi = (hi << 7) | (unsigned)(lo >> 57);
	lo &= 0x1FFFFFFFFFFFFFF;

	/*
	 * Broadcast the high part and compare it with the relevant
	 * values. We need both a "greater than" and an "equal"
	 * comparisons.
	 */
	xhi = _mm256_broadcastw_epi16(_mm_cvtsi32_si128(hi));
	rhi = _mm256_loadu_si256(&rhi15.ymm[0]);
	gthi = _mm256_cmpgt_epi16(rhi, xhi);
	eqhi = _mm256_cmpeq_epi16(rhi, xhi);

	/*
	 * The result is the number of 72-bit values (among the list of 19)
	 * which are greater than the 72-bit random value. We first count
	 * all non-zero 16-bit elements in the first eight of gthi. Such
	 * elements have value -1 or 0, so we first negate them.
	 */
	t = _mm_srli_epi16(_mm256_castsi256_si128(gthi), 15);
	zt = _mm_setzero_si128();
	t = _mm_hadd_epi16(t, zt);
	t = _mm_hadd_epi16(t, zt);
	t = _mm_hadd_epi16(t, zt);
	r = _mm_cvtsi128_si32(t);

	/*
	 * We must look at the low bits for all values for which the
	 * high bits are an "equal" match; values 8-18 all have the
	 * same high bits (0).
	 * On 32-bit systems, 'lo' really is two registers, requiring
	 * some extra code.
	 */
#if defined(__x86_64__) || defined(_M_X64)
	xlo = _mm256_broadcastq_epi64(_mm_cvtsi64_si128(*(int64_t *)&lo));
#else
	{
		uint32_t e0, e1;
		int32_t f0, f1;

		e0 = (uint32_t)lo;
		e1 = (uint32_t)(lo >> 32);
		f0 = *(int32_t *)&e0;
		f1 = *(int32_t *)&e1;
		xlo = _mm256_set_epi32(f1, f0, f1, f0, f1, f0, f1, f0);
	}
#endif
	gtlo0 = _mm256_cmpgt_epi64(_mm256_loadu_si256(&rlo57.ymm[0]), xlo); 
	gtlo1 = _mm256_cmpgt_epi64(_mm256_loadu_si256(&rlo57.ymm[1]), xlo); 
	gtlo2 = _mm256_cmpgt_epi64(_mm256_loadu_si256(&rlo57.ymm[2]), xlo); 
	gtlo3 = _mm256_cmpgt_epi64(_mm256_loadu_si256(&rlo57.ymm[3]), xlo); 
	gtlo4 = _mm256_cmpgt_epi64(_mm256_loadu_si256(&rlo57.ymm[4]), xlo); 

	/*
	 * Keep only comparison results that correspond to the non-zero
	 * elements in eqhi.
	 */
	gtlo0 = _mm256_and_si256(gtlo0, _mm256_cvtepi16_epi64(
		_mm256_castsi256_si128(eqhi)));
	gtlo1 = _mm256_and_si256(gtlo1, _mm256_cvtepi16_epi64(
		_mm256_castsi256_si128(_mm256_bsrli_epi128(eqhi, 8))));
	eqm = _mm256_permute4x64_epi64(eqhi, 0xFF);
	gtlo2 = _mm256_and_si256(gtlo2, eqm);
	gtlo3 = _mm256_and_si256(gtlo3, eqm);
	gtlo4 = _mm256_and_si256(gtlo4, eqm);

	/*
	 * Add all values to count the total number of "-1" elements.
	 * Since the first eight "high" words are all different, only
	 * one element (at most) in gtlo0:gtlo1 can be non-zero; however,
	 * if the high word of the random value is zero, then many
	 * elements of gtlo2:gtlo3:gtlo4 can be non-zero.
	 */
	gtlo0 = _mm256_or_si256(gtlo0, gtlo1);
	gtlo0 = _mm256_add_epi64(
		_mm256_add_epi64(gtlo0, gtlo2),
		_mm256_add_epi64(gtlo3, gtlo4));
	t = _mm_add_epi64(
		_mm256_castsi256_si128(gtlo0),
		_mm256_extracti128_si256(gtlo0, 1));
	t = _mm_add_epi64(t, _mm_srli_si128(t, 8));
	r -= _mm_cvtsi128_si32(t);

	return r;

}

/*
 * Sample a bit with probability exp(-x) for some x >= 0.
 */
TARGET_AVX2
static int
BerExp(prng *p, fpr x, fpr ccs)
{
	int s, i;
	fpr r;
	uint32_t sw, w;
	uint64_t z;

	/*
	 * Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
	 * and 0 <= r < log(2). Since x >= 0, we can use fpr_trunc().
	 */
	s = (int)fpr_trunc(fpr_mul(x, fpr_inv_log2));
	r = fpr_sub(x, fpr_mul(fpr_of(s), fpr_log2));

	/*
	 * It may happen (quite rarely) that s >= 64; if sigma = 1.2
	 * (the minimum value for sigma), r = 0 and b = 1, then we get
	 * s >= 64 if the half-Gaussian produced a z >= 13, which happens
	 * with probability about 0.000000000230383991, which is
	 * approximatively equal to 2^(-32). In any case, if s >= 64,
	 * then BerExp will be non-zero with probability less than
	 * 2^(-64), so we can simply saturate s at 63.
	 */
	sw = (uint32_t)s;
	sw ^= (sw ^ 63) & -((63 - sw) >> 31);
	s = (int)sw;

	/*
	 * Compute exp(-r); we know that 0 <= r < log(2) at this point, so
	 * we can use fpr_expm_p63(), which yields a result scaled to 2^63.
	 * We scale it up to 2^64, then right-shift it by s bits because
	 * we really want exp(-x) = 2^(-s)*exp(-r).
	 *
	 * The "-1" operation makes sure that the value fits on 64 bits
	 * (i.e. if r = 0, we may get 2^64, and we prefer 2^64-1 in that
	 * case). The bias is negligible since fpr_expm_p63() only computes
	 * with 51 bits of precision or so.
	 */
	z = ((fpr_expm_p63(r, ccs) << 1) - 1) >> s;

	/*
	 * Sample a bit with probability exp(-x). Since x = s*log(2) + r,
	 * exp(-x) = 2^-s * exp(-r), we compare lazily exp(-x) with the
	 * PRNG output to limit its consumption, the sign of the difference
	 * yields the expected result.
	 */
	i = 64;
	do {
		i -= 8;
		w = prng_get_u8(p) - ((uint32_t)(z >> i) & 0xFF);
	} while (!w && i > 0);
	return (int)(w >> 31);
}

/*
 * The sampler produces a random integer that follows a discrete Gaussian
 * distribution, centered on mu, and with standard deviation sigma. The
 * provided parameter isigma is equal to 1/sigma.
 *
 * The value of sigma MUST lie between 1 and 2 (i.e. isigma lies between
 * 0.5 and 1); in Falcon, sigma should always be between 1.2 and 1.9.
 */
TARGET_AVX2
int
Zf(sampler)(void *ctx, fpr mu, fpr isigma)
{
	sampler_context *spc;
	int s;
	fpr r, dss, ccs;

	spc = (sampler_context *)ctx;

	/*
	 * Center is mu. We compute mu = s + r where s is an integer
	 * and 0 <= r < 1.
	 */
	s = (int)fpr_floor(mu);
	r = fpr_sub(mu, fpr_of(s));

	/*
	 * dss = 1/(2*sigma^2) = 0.5*(isigma^2).
	 */
	dss = fpr_half(fpr_sqr(isigma));

	/*
	 * ccs = sigma_min / sigma = sigma_min * isigma.
	 */
	ccs = fpr_mul(isigma, spc->sigma_min);

	/*
	 * We now need to sample on center r.
	 */
	for (;;) {
		int z0, z, b;
		fpr x;

		/*
		 * Sample z for a Gaussian distribution. Then get a
		 * random bit b to turn the sampling into a bimodal
		 * distribution: if b = 1, we use z+1, otherwise we
		 * use -z. We thus have two situations:
		 *
		 *  - b = 1: z >= 1 and sampled against a Gaussian
		 *    centered on 1.
		 *  - b = 0: z <= 0 and sampled against a Gaussian
		 *    centered on 0.
		 */
		z0 = Zf(gaussian0_sampler)(&spc->p);
		b = (int)prng_get_u8(&spc->p) & 1;
		z = b + ((b << 1) - 1) * z0;

		/*
		 * Rejection sampling. We want a Gaussian centered on r;
		 * but we sampled against a Gaussian centered on b (0 or
		 * 1). But we know that z is always in the range where
		 * our sampling distribution is greater than the Gaussian
		 * distribution, so rejection works.
		 *
		 * We got z with distribution:
		 *    G(z) = exp(-((z-b)^2)/(2*sigma0^2))
		 * We target distribution:
		 *    S(z) = exp(-((z-r)^2)/(2*sigma^2))
		 * Rejection sampling works by keeping the value z with
		 * probability S(z)/G(z), and starting again otherwise.
		 * This requires S(z) <= G(z), which is the case here.
		 * Thus, we simply need to keep our z with probability:
		 *    P = exp(-x)
		 * where:
		 *    x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
		 *
		 * Here, we scale up the Bernouilli distribution, which
		 * makes rejection more probable, but makes rejection
		 * rate sufficiently decorrelated from the Gaussian
		 * center and standard deviation that the whole sampler
		 * can be said to be constant-time.
		 */
		x = fpr_mul(fpr_sqr(fpr_sub(fpr_of(z), r)), dss);
		x = fpr_sub(x, fpr_mul(fpr_of(z0 * z0), fpr_inv_2sqrsigma0));
		if (BerExp(&spc->p, x, ccs)) {
			/*
			 * Rejection sampling was centered on r, but the
			 * actual center is mu = s + r.
			 */
			return s + z;
		}
	}
}
