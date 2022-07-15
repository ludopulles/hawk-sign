/*
 * Generates a list of signatures
 */
#include <cassert>
#include <climits>
#include <cstdio>
#include <vector>
#include <mutex>
#include <thread>

// x86_64 specific:
#include<sys/time.h>

extern "C" {
	#define restrict
	#include "../inner.h"
}

#define MAXLOGN (10)
#define MAXN MKN(MAXLOGN)

constexpr size_t num_sigs = 1000;
constexpr size_t num_threads = 4;
constexpr size_t signs_per_kg = 10;

int16_t sig[num_sigs][MAXN];

void generate_sigs(unsigned logn, const unsigned char *seed, int seed_len, size_t from, size_t to)
{
	union {
		uint8_t b[44 * MAXN];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;

	int8_t f[MAXN], g[MAXN], F[MAXN], G[MAXN];
	uint8_t h[MAXN/4];
	int16_t iq00[MAXN], iq01[MAXN];
	fpr exp_sk[EXPANDED_SECKEY_SIZE(MAXLOGN)];
	inner_shake256_context sc;

	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, seed_len);
	inner_shake256_flip(&sc);

	int nth = signs_per_kg - 1;
	// Generate key pair.
	for (size_t i = from; i < to; i++) {
		if (++nth == signs_per_kg) {
			Zf(keygen)(&sc, f, g, F, G, iq00, iq01, logn, tmp.b);
			Zf(expand_seckey)(exp_sk, f, g, F, logn);
			nth = 0;
		}

		// make a signature of a random message...
		inner_shake256_extract(&sc, h, sizeof h);

		// Compute the signature.
		while (!Zf(sign)(&sc, sig[i], exp_sk, h, logn, tmp.b)) {}
	}
}

int main(int argc, char **argv) {
	unsigned char seed[48] = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	inner_shake256_context sc;

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);

	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	std::thread* pool[num_threads - 1];
	unsigned char seeds[num_threads * 48];

	inner_shake256_extract(&sc, seeds, num_threads * 48);

	size_t job_parts[num_threads + 1];
	job_parts[0] = 0;
	for (size_t i = 1; i <= num_threads; i++) {
		job_parts[i] = i * num_sigs / num_threads;
	}

	if (argc != 2) {
		printf("Usage: ./file logn\n");
		return 1;
	}
	unsigned logn = atoi(argv[1]);
	assert(1 <= logn && logn <= MAXLOGN);

	for (size_t i = 1; i < num_threads; i++) {
		pool[i - 1] = new std::thread(generate_sigs, logn, seeds + i * 48, 48, job_parts[i], job_parts[i + 1]);
	}
	generate_sigs(logn, seeds, 48, job_parts[0], job_parts[1]);

	for (size_t i = 0; i < num_threads - 1; i++) {
		pool[i]->join();
		delete pool[i];
	}

	for (size_t i = 0; i < num_sigs; i++) {
		for (size_t u = 0; u < MKN(logn); u++)
			printf("%d ", sig[i][u]);
		printf("\n");
	}

	return 0;
}
