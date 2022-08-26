// c99 -W -Wall -O3 -march=native test_sizes.c -o test_sizes build/* -lm

#include <cassert>
#include <cstdint>
#include <cstdio>

#include <mutex>
#include <thread>

extern "C" {
	#define restrict
	#include "../inner.h"
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
#define MAXLOGN (10)
#define MAXN MKN(MAXLOGN)

size_t sksize[MAXLOGN+1];

void find_seckeysizes(inner_shake256_context *sc) {
	uint8_t b[48 << MAXLOGN];
	int8_t f[MAXN], g[MAXN], F[MAXN], G[MAXN];
	int16_t iq00[MAXN], iq01[MAXN];

	for (unsigned logn = 1; logn <= MAXLOGN; logn++) {
		// Generate key pair.
		Zf(keygen)(sc, f, g, F, G, iq00, iq01, logn, b);

		// Take the header byte into account:
		sksize[logn] = 1 + Zf(encode_seckey)(NULL, 0, f, g, F, logn);

		for (int test = 10; test--; ) {
			// Check if it is the same a second time
			Zf(keygen)(sc, f, g, F, G, iq00, iq01, logn, b);
			assert(sksize[logn] == 1 + Zf(encode_seckey)(NULL, 0, f, g, F, logn));
		}
	}
}

constexpr int nrepetitions = 125;
constexpr int nthreads = 4;

struct WorkerResult {
	long long iterations;
	long long pksum, pksq, pkmin, pkmax;
	long long fgsum, fgsq, FGsum, FGsq;
	long long q00sum, q00sq;
	long long q01sum, q01sq;
	long long q11sum, q11sq;

	WorkerResult() : iterations(0),
		pksum(0), pksq(0), pkmin(1000000), pkmax(0),
		fgsum(0), fgsq(0), FGsum(0), FGsq(0),
		q00sum(0), q00sq(0), q01sum(0), q01sq(0), q11sum(0), q11sq(0) {}

	void combine(const WorkerResult &res) {
		iterations += res.iterations;

		pksum += res.pksum;
		pksq += res.pksq;
		if (res.pkmin < pkmin) pkmin = res.pkmin;
		if (res.pkmax > pkmax) pkmax = res.pkmax;
		fgsum += res.fgsum;
		fgsq += res.fgsq;
		FGsum += res.FGsum;
		FGsq += res.FGsq;

		q00sum += res.q00sum;
		q00sq += res.q00sq;
		q01sum += res.q01sum;
		q01sq += res.q01sq;
		q11sum += res.q11sum;
		q11sq += res.q11sq;
	}
};

WorkerResult run(unsigned logn, const unsigned char *seed, size_t seed_len)
{
	union {
		uint8_t b[44 * MAXN];
		uint64_t dummy_u64;
		fpr dummy_fpr;
	} tmp;

	int8_t f[MAXN], g[MAXN], F[MAXN], G[MAXN];
	fpr q00[MAXN], q01[MAXN], q11[MAXN];
	int16_t iq00[MAXN], iq01[MAXN];
	inner_shake256_context sc;
	WorkerResult res;

	// Initialize a RNG.
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, seed_len);
	inner_shake256_flip(&sc);

	for (size_t i = 0; i < nrepetitions; i++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, iq00, iq01, logn, tmp.b);

		for (size_t u = 0; u < MKN(logn); u++) {
			res.fgsq += f[u]*f[u] + g[u]*g[u];
			res.FGsq += F[u]*F[u] + G[u]*G[u];
		}

		assert(iq00[0] >= 0);
		// Pay attention: q11 does NOT fit into a int16_t.

		// Take the header byte into account:
		int pk_sz = 1 + Zf(encode_pubkey)(NULL, 0, iq00, iq01, logn);

		res.pksum += pk_sz;
		res.pksq += pk_sz * pk_sz;
		if (pk_sz < res.pkmin) res.pkmin = pk_sz;
		if (pk_sz > res.pkmax) res.pkmax = pk_sz;

		Zf(complete_pubkey)(iq00, iq01, q00, q01, q11, logn);
		Zf(iFFT)(q00, logn);
		Zf(iFFT)(q01, logn);
		Zf(iFFT)(q11, logn);

		for  (size_t u = 1; u < MKN(logn - 1); u++) {
			long long x = fpr_rint(q00[u]);
			long long z = fpr_rint(q11[u]);
			res.q00sum += x;
			res.q00sq += x*x;
			res.q11sum += z;
			res.q11sq += z*z;
		}

		for  (size_t u = 0; u < MKN(logn); u++) {
			long long y = fpr_rint(q01[u]);
			res.q01sum += y;
			res.q01sq += y*y;
		}
	}

	res.iterations = nrepetitions;
	return res;
}

std::mutex mx;
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

	// Check secret key sizes quickly
	find_seckeysizes(&sc);

	std::thread* pool[nthreads-1];
	unsigned char seeds[nthreads * 48];
	inner_shake256_extract(&sc, seeds, nthreads * 48);

#define VAR0(sq, N) (N == 0 ? 0.0 : sqrt((double)(sq) / (double)(N)))

	WorkerResult totals[MAXLOGN + 1];
	printf("logn | sigma f,g | sigma F,G | sigma q00 | sigma q01 | sigma q11 | |sk|\n");
	for (unsigned logn = 1; logn <= MAXLOGN; logn++) {
		tot = WorkerResult();

		// Do the work
		for (int i = 0; i < nthreads-1; i++)
			pool[i] = new std::thread(work, logn, seeds + i * 48, 48);
		work(logn, seeds + (nthreads - 1) * 48, 48);
		for (int i = 0; i < nthreads-1; i++)
			pool[i]->join(), delete pool[i];

		// Store the result
		totals[logn] = tot;

		// Print current values
		int n = MKN(logn);
		int ncc = n/2 - 1, reps = totals[logn].iterations;
		// printf("logn = %d: sigma_{f,g} ~ %.8f, sigma_{F,G} ~ %.8f, sig_{q00} ~ %.8f, sig_{q01} ~ %.8f, sig_{q11} ~ %.8f\n",
		printf("%4d | %9.7f | %9.6f | %9.4f | %9.4f | %9.3f | %4d\n",
			(int) logn,
			VAR0(totals[logn].fgsq, 2 * reps * n),
			VAR0(totals[logn].FGsq, 2 * reps * n),
			VAR0(totals[logn].q00sq, reps * ncc),
			VAR0(totals[logn].q01sq, reps * n),
			VAR0(totals[logn].q11sq, reps * ncc),
			(int) sksize[logn]
		);
		fflush(stdout);
	}

	printf("\nIterations performed: %u\n", (unsigned) totals[MAXLOGN].iterations);
	printf("-----+- Public key -------+------+------+---------\n");
	printf("logn | Average +/- stddev | min  | max  | 6sigma  \n");
	double avg, std;

	// Public key
	for (size_t logn = 1; logn <= MAXLOGN; logn++) {
		int reps = totals[logn].iterations;
		avg = (double) totals[logn].pksum / reps;
		std = sqrt( (double) totals[logn].pksq / reps - avg*avg );
		printf("%4d | %7.2f +/- %6.2f | %4lld | %4lld | %7.2f\n",
			(int) logn, avg, std, totals[logn].pkmin, totals[logn].pkmax, avg + 6 * std);
	}
	return 0;
}
