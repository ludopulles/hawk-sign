#include <assert.h>
#include <stdio.h>
// x86_64 specific:
#include <sys/time.h>

#include "keygen.c"

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

void fpr_to_int16(int16_t *buf, fpr *p, size_t logn) {
	unsigned n, u;

	n = MKN(logn);
	for (u = 0; u < n; u ++) {
		int val = fpr_rint(p[u]);
		assert(-(1 << 15) <= val && val < (1 << 15));
		buf[u] = (int16_t) val;
	}
}

void poly_output(fpr *p, size_t logn) {
	/* fpr sqnorm = fpr_zero;
	for (size_t u = 0; u < MKN(logn); u++)
		sqnorm = fpr_add(sqnorm, fpr_sqr(p[u]));
	
	printf("norm %ld ", fpr_rint(fpr_sqrt(sqnorm))); */

	for (size_t u = 0; u < MKN(logn); u++) {
		if (u) printf(" ");
		printf("%ld", fpr_rint(p[u]));
	}
	printf("\n");
}

/* see inner.h */
size_t
Zf(comp_encode_q00)(
	void *out, size_t max_out_len,
	const int16_t *x, unsigned logn)
{
	uint8_t *buf;
	size_t n, u, v;
	uint32_t acc;
	unsigned acc_len;

	n = (size_t)1 << logn;
	buf = (uint8_t *)out;

	/*
	 * Make sure that all values are within the -2047..+2047 range.
	 * Moreover, q00 = q00^* so q00 + q00^* must be in fact an integer.
	 */
	for (u = 1; u < n/2; u ++) {
		if (x[u] < -2047 || x[u] > +2047 || x[u] + x[n-u] != 0) {
			return 0;
		}
	}

	// acc = 0;
	// acc_len = 0;
	v = 0;

	/*
	 * The first value of q00 is of the order sigma_kg^2 2d,
	 * but definitely << 8d when sigma_kg = 1.425, by the standard
	 * Laurent-Massart bound (see https://en.wikipedia.org/wiki/Chi-squared_distribution#Concentration).
	 */
	acc_len = 3 + logn;
	acc = (uint32_t) x[0];

	int32_t max_val = 1 << acc_len;
	if (x[0] <= -max_val || x[0] >= max_val) {
		return 0;
	}

	/*
	 * Produce all full bytes.
	 */
	while (acc_len >= 8) {
		acc_len -= 8;
		if (buf != NULL) {
			if (v >= max_out_len) {
				return 0;
			}
			buf[v] = (uint8_t)(acc >> acc_len);
		}
		v ++;
	}

	for (u = 1; u < n/2; u ++) {
		int t;
		unsigned w;

		/*
		 * Get sign and absolute value of next integer; push the
		 * sign bit.
		 */
		acc <<= 1;
		t = x[u];
		if (t < 0) {
			t = -t;
			acc |= 1;
		}
		w = (unsigned)t;

		/*
		 * Push the low `lo_bits` bits of the absolute value.
		 */

		const int lo_bits = 5;
		acc <<= lo_bits;
		acc |= w & ((1U << lo_bits) - 1);
		w >>= lo_bits;

		/*
		 * We pushed exactly `lo_bits + 1` bits.
		 */
		acc_len += (lo_bits + 1);

		/*
		 * Push as many zeros as necessary, then a one. Since the
		 * absolute value is at most 2047, w can only range up to
		 * 15 at this point, thus we will add at most 16 bits
		 * here. With the 8 bits above and possibly up to 7 bits
		 * from previous iterations, we may go up to 31 bits, which
		 * will fit in the accumulator, which is an uint32_t.
		 */
		acc <<= (w + 1);
		acc |= 1;
		acc_len += w + 1;

		/*
		 * Produce all full bytes.
		 */
		while (acc_len >= 8) {
			acc_len -= 8;
			if (buf != NULL) {
				if (v >= max_out_len) {
					return 0;
				}
				buf[v] = (uint8_t)(acc >> acc_len);
			}
			v ++;
		}
	}

	/*
	 * Flush remaining bits (if any).
	 */
	if (acc_len > 0) {
		if (buf != NULL) {
			if (v >= max_out_len) {
				return 0;
			}
			buf[v] = (uint8_t)(acc << (8 - acc_len));
		}
		v ++;
	}

	return v;
}

/* see inner.h */
size_t
Zf(comp_encode_q10)(
	void *out, size_t max_out_len,
	const int16_t *x, unsigned logn)
{
	uint8_t *buf;
	size_t n, u, v;
	uint64_t acc;
	unsigned acc_len;

	n = (size_t)1 << logn;
	buf = (uint8_t *)out;

	/*
	 * Make sure that all values are within the -4095..+4095 range.
	 */
	for (u = 0; u < n; u ++) {
		if (x[u] < -4095 || x[u] > +4095) {
			return 0;
		}
	}

	acc = 0;
	acc_len = 0;
	v = 0;
	for (u = 0; u < n; u ++) {
		int t;
		unsigned w;

		/*
		 * Get sign and absolute value of next integer; push the
		 * sign bit.
		 */
		acc <<= 1;
		t = x[u];
		if (t < 0) {
			t = -t;
			acc |= 1;
		}
		w = (unsigned)t;

		/*
		 * Push the low `lo_bits` bits of the absolute value.
		 */

		const int lo_bits = 9;

		acc <<= lo_bits;
		acc |= w & ((1U << lo_bits) - 1);
		w >>= lo_bits;

		/*
		 * We pushed exactly `lo_bits + 1` bits.
		 */
		acc_len += (lo_bits + 1);

		// TODO: perhaps this still works, but perhaps we need uint64_t...

		/*
		 * Push as many zeros as necessary, then a one. Since the
		 * absolute value is at most 4095, w can only range up to
		 * 7 at this point, thus we will add at most 8 bits
		 * here. With the 10 bits above and possibly up to 7 bits
		 * from previous iterations, we may go up to 25 bits, which
		 * will fit in the accumulator, which is an uint32_t.
		 */
		acc <<= (w + 1);
		acc |= 1;
		acc_len += w + 1;

		/*
		 * Produce all full bytes.
		 */
		while (acc_len >= 8) {
			acc_len -= 8;
			if (buf != NULL) {
				if (v >= max_out_len) {
					return 0;
				}
				buf[v] = (uint8_t)(acc >> acc_len);
			}
			v ++;
		}
	}

	/*
	 * Flush remaining bits (if any).
	 */
	if (acc_len > 0) {
		if (buf != NULL) {
			if (v >= max_out_len) {
				return 0;
			}
			buf[v] = (uint8_t)(acc << (8 - acc_len));
		}
		v ++;
	}

	return v;
}

/* see inner.h for keygen */
void
keygen_count_fails(inner_shake256_context *rng,
	int8_t *restrict f, int8_t *restrict g, // secret key
	int8_t *restrict F, int8_t *restrict G, // secret key
	fpr *restrict q00, fpr *restrict q10, fpr *restrict q11, // public key
	unsigned logn, uint8_t *restrict tmp, fpr isigma_kg, int *num_fails)
{
	size_t n = MKN(logn);
	for (;;) {
		fpr *rt1, *rt2;
		sampler_context spc;
		void *samp_ctx;
		spc.sigma_min = fpr_sigma_min[logn];
		Zf(prng_init)(&spc.p, rng);
		samp_ctx = &spc;
		poly_small_mkgauss(samp_ctx, f, logn, isigma_kg, 128);
		poly_small_mkgauss(samp_ctx, g, logn, isigma_kg, 128);
		if (!solve_NTRU(logn, F, G, f, g, 128, (uint32_t *)tmp)) {
			(*num_fails)++;
			continue;
		}

		rt1 = (fpr *)tmp;
		rt2 = rt1 + n;
		poly_small_to_fp(q00, f, logn);
		poly_small_to_fp(rt1, g, logn);
		poly_small_to_fp(q11, F, logn);
		poly_small_to_fp(rt2, G, logn);
		Zf(FFT)(q00, logn); // f
		Zf(FFT)(rt1, logn); // g
		Zf(FFT)(q11, logn); // F
		Zf(FFT)(rt2, logn); // G
		Zf(poly_add_muladj_fft)(q10, q11, rt2, q00, rt1, logn);
		Zf(poly_mulselfadj_fft)(q00, logn); // f*adj(f)
		Zf(poly_mulselfadj_fft)(rt1, logn); // g*adj(g)
		Zf(poly_add)(q00, rt1, logn);
		Zf(poly_mulselfadj_fft)(q11, logn); // F*adj(F)
		Zf(poly_mulselfadj_fft)(rt2, logn); // G*adj(G)
		Zf(poly_add)(q11, rt2, logn);
		break;
	}
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
const size_t logn = 9, n = MKN(logn);

void measure_keygen(fpr isigma_kg) {
	uint8_t b[28 << logn]; // 14 kB temporary memory, 17.5 kB total
	int8_t f[n], g[n], F[n], G[n];
	fpr q00[n], q10[n], q11[n];
	fpr q00i[n], q10i[n], q11i[n];
	int16_t q00n[n], q10n[n];
	unsigned char seed[48];
	inner_shake256_context sc;

	struct timeval t0, t1;
	const int n_repetitions = 1000;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	gettimeofday(&t0, NULL);

	/* long long sumval[16], sumsqval[16];
	memset(sumval, 0, sizeof sumval);
	memset(sumsqval, 0, sizeof sumsqval); */

	int fails = 0;
	for (int i = 0; i < n_repetitions; i++) {
		// Generate key pair.
		keygen_count_fails(&sc, f, g, F, G, q00, q10, q11, logn, b, isigma_kg, &fails);

		memcpy(q00i, q00, sizeof(q00));
		memcpy(q10i, q10, sizeof(q10));
		memcpy(q11i, q11, sizeof(q11));

		Zf(iFFT)(q00i, logn);
		Zf(iFFT)(q10i, logn);
		Zf(iFFT)(q11i, logn);

		fpr_to_int16(q00n, q00i, logn);
		fpr_to_int16(q10n, q10i, logn);

		// poly_output(q00i, logn);
		// poly_output(q10i, logn);
		// poly_output(q11i, logn);

		assert(q00n[n / 2] == 0);
		for (size_t u = 1; u < n; u++) {
			assert(q00n[n - u] == -q00n[u]);
		}

		size_t pubkey_sz =
			Zf(comp_encode_q00)(NULL, 0, q00n, logn) +
			Zf(comp_encode_q10)(NULL, 0, q10n, logn);
		if (pubkey_sz >= 1024)
			printf("Public key size (bits): \t%zu\n", pubkey_sz);

		/* for (size_t lobits = 1; lobits < 16; lobits++) {
			size_t s = Zf(comp_encode_q10)(NULL, 0, q10n, logn, lobits);
			if (s == 0) {
				printf("Failed for lobits=%zu\n", lobits);
			}
			sumval[lobits] += s;
			sumsqval[lobits] += s*s;
		} */
	}

	gettimeofday(&t1, NULL);
	double kg_duration = (double)time_diff(&t0, &t1) / n_repetitions; // (in us)
	printf("Average time per keygen: %.3f ms\n", kg_duration / 1000.0);
	// This requires catching failed attempts
	printf("Probability failure: %.2f%%\n", 100.0 * fails / n_repetitions);

/*
	// Gather statistics
	for (size_t lobits = 1; lobits < 16; lobits++) {
		double avg = ((double) sumval[lobits]) / n_repetitions;
		double var = ((double) sumsqval[lobits]) / n_repetitions - avg*avg;
		printf("Average (lobits=%zu): %.2f\t stddev %.2f\n", lobits, avg, sqrt(var));
	}
*/
}

int8_t valid_sigma(fpr sigma_sig) {
	return !fpr_lt(sigma_sig, fpr_sigma_min[logn])
		&& fpr_lt(sigma_sig, fpr_div(fpr_of(18205), fpr_of(10000)));
}

int main() {
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	const fpr sigma_kg  = fpr_div(fpr_of(1425), fpr_of(1000));
	assert(valid_sigma(sigma_kg));
	fpr isigma_kg = fpr_inv(sigma_kg);
	measure_keygen(isigma_kg);
	return 0;
}
