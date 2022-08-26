/*
 * Tests:
 * - if signing gives a valid signature,
 * - if (0, 0, .... 0) is not a valid signature.
 */
#include <cassert>
#include <cstdint>
#include <cstdio>

#include <mutex>
#include <thread>

extern "C" {
	#define restrict
	#include "../inner.h"
}

#define MAXLOGN (10)
#define MAXN MKN(MAXLOGN)

struct WorkerResult {
	long long iterations, fails;

	WorkerResult() : iterations(0), fails(0) {}

	void combine(const WorkerResult &res) {
		iterations += res.iterations;
		fails += res.fails;
	}
};

const long long num_kg = 50, signs_per_kg = 10;

void print_poly(const int8_t *p, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	printf("{ ");
	for (u = 0; u < n; u++) {
		if (u) printf(",");
		printf("%d", p[u]);
	}
	printf(" },\n");
}

std::mutex mx;

WorkerResult run(unsigned logn, const unsigned char *seed, size_t seed_len)
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
	int16_t sig[MAXN];
	inner_shake256_context sc;
	WorkerResult result;

	// Initialize a RNG.
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, seed_len);
	inner_shake256_flip(&sc);

	// Generate key pair.
	for (size_t u = 0; u < num_kg; u++) {
		Zf(keygen)(&sc, f, g, F, G, iq00, iq01, logn, tmp.b);
		Zf(expand_seckey)(exp_sk, f, g, F, logn);

		long long last_fail = result.fails;
		for (size_t rep = 0; rep < signs_per_kg; rep++) {
			// Make a signature of a random message.
			inner_shake256_extract(&sc, h, sizeof h);

			// Compute the signature.
			while (!Zf(sign)(&sc, sig, exp_sk, h, logn, tmp.b)) {}

			result.fails += !Zf(verify_NTT)(h, sig, iq00, iq01, logn, tmp.b);
		}

		if (last_fail < result.fails) {
			// There were one or more failures, maybe this depends on the public key?

			/* acquire mutex lock */
			std::lock_guard<std::mutex> guard(mx);

			printf("\n");
			print_poly(f, logn);
			print_poly(g, logn);
			print_poly(F, logn);
			print_poly(G, logn);
			printf("\n");
			fflush(stdout);

			/* release mutex lock */
		}
	}

	result.iterations = num_kg * signs_per_kg;
	return result;
}

WorkerResult tot;

void work(unsigned logn, const unsigned char *seed, int seed_len)
{
	WorkerResult result = run(logn, seed, seed_len);

	/* acquire mutex lock */ {
		std::lock_guard<std::mutex> guard(mx);
		tot.combine(result);
	}
}

int main() {
	unsigned char seed[48] = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	inner_shake256_context sc;

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);

	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	const int nthreads = 2;
	std::thread* pool[nthreads-1];
	unsigned char seeds[nthreads * 48];

	inner_shake256_extract(&sc, seeds, nthreads * 48);

	time_t my_time = time(NULL);
	printf("%s", ctime(&my_time));
	fflush(stdout);

	for (unsigned logn = 1; logn <= MAXLOGN; logn++) {
		tot = WorkerResult();

		for (int i = 0; i < nthreads-1; i++) {
			pool[i] = new std::thread(work, logn, seeds + i * 48, 48);
		}
		work(logn, seeds + (nthreads - 1) * 48, 48);
		for (int i = 0; i < nthreads-1; i++) pool[i]->join(), delete pool[i];

		printf("%u: %lld/%lld ~ %.6f%% fails\n", logn, tot.fails, tot.iterations, 100.0 * tot.fails / tot.iterations);

		my_time = time(NULL);
		printf("%s", ctime(&my_time));
		fflush(stdout);
	}

	return 0;
}
