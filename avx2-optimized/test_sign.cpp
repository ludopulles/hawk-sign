#include <cassert>
#include <cstdio>
#include <vector>
#include <atomic>
#include <thread>

// x86_64 specific:
#include<sys/time.h>

extern "C" {
	#ifndef restrict
		#define restrict
	#endif

	#include "inner.h"
}

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

long long time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

/* see inner.h */
size_t
Zf(comp_encode)(
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
	 */
	for (u = 0; u < n; u ++) {
		if (x[u] < -2047 || x[u] > +2047) {
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

		constexpr int lo_bits = 5;

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


// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
constexpr size_t logn = 9, n = MKN(logn);

struct WorkerResult {
	int num_signed;
	int num_invalid;
	int num_babai_fail;
	int num_sig_differ;

	WorkerResult() =default;
};

WorkerResult measure_signatures(fpr isigma_kg, fpr isigma_sig, fpr verif_bound) {
	uint8_t b[42 << logn];
	int8_t f[n], g[n], F[n], G[n], h[n];
	int16_t s0[n], reconstructed_s0[n], s1[n];
	fpr q00[n], q10[n], q11[n];
	unsigned char seed[48];
	inner_shake256_context sc;
	const int n_repetitions = 1000;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	WorkerResult result;
	result.num_signed = n_repetitions;

	for (int rep = 0; rep < n_repetitions; rep++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, logn, b, isigma_kg);

		// make a signature of a random message...
		randombytes((unsigned char *)h, sizeof h);

		// Compute the signature.
		Zf(complete_sign)(&sc, s0, s1, f, g, F, G, h, logn, isigma_sig, b);

		printf("%zu ", Zf(comp_encode)(NULL, 0, s1, logn));
		fflush(stdout);

		if (!Zf(complete_verify)(h, s0, s1, q00, q10, q11, logn, verif_bound, b))
			result.num_invalid++;

		if (!Zf(verify)(h, reconstructed_s0, s1, q00, q10, q11, logn, verif_bound, b))
			result.num_babai_fail++;
		else
			assert(Zf(complete_verify)(h, reconstructed_s0, s1, q00, q10, q11, logn, verif_bound, b));

		for (size_t u = 0; u < n; u++) {
			if (s0[u] != reconstructed_s0[u]) {
				result.num_sig_differ++;
				break;
			}
		}
	}
	return result;
}

std::atomic<int> tot_signed(0), tot_invalid(0), tot_babai_fail(0), tot_sig_differ(0);

constexpr fpr sigma_kg  = { v: 1.425 };
constexpr fpr sigma_sig = { v: 1.292 };
constexpr fpr verif_margin = { v: sigma_kg.v / sigma_sig.v };

fpr getverif_bound() {
	return fpr_mul(fpr_sqr(fpr_mul(verif_margin, fpr_double(sigma_sig))), fpr_double(fpr_sqr(fpr_of(n))));
}

void work() {
	WorkerResult result = measure_signatures(fpr_inv(sigma_kg), fpr_inv(sigma_sig), getverif_bound());

	tot_signed += result.num_signed;
	tot_invalid += result.num_invalid;
	tot_babai_fail += result.num_babai_fail;
	tot_sig_differ += result.num_sig_differ;;
}

int main() {
	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	unsigned seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %u\n", seed);
	srand(seed);

	int nthreads = 4;
	if (nthreads == 1) {
		work();
	} else {
		std::vector<std::thread*> pool(nthreads);
		for (int i = 0; i < nthreads; i++) {
			pool[i] = new std::thread(work);
		}

		for (int i = 0; i < nthreads; i++) {
			pool[i]->join();
		}
	}

	printf("# Signatures signed:      %d\n", static_cast<int>(tot_signed));
	printf("# Signatures invalid:     %d\n", static_cast<int>(tot_invalid));
	printf("# Babai roundings failed: %d\n", static_cast<int>(tot_babai_fail));
	printf("# Babai != original s0:   %d\n", static_cast<int>(tot_sig_differ));

	return 0;
}
