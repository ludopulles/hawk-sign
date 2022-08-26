#include <cassert>
#include <climits>
#include <cstdio>
#include <fstream>
#include <vector>
#include <string>

#include <map>
#include <unordered_map>

// x86_64 specific:
#include<sys/time.h>

extern "C" {
	#define restrict
	#include "../inner.h"
}

typedef long long ll;
using namespace std;

namespace std {
	template<class T1, class T2>
	struct hash<pair<T1,T2>> {
	public:
		size_t operator()(const pair<T1,T2> &p) const {
			size_t x = hash<T1>()(p.first), y = hash<T2>()(p.second);
			return x ^ (y + 0x9e3779b9 + (x<<6) + (x>>2));
		}
	};

	template<class T>
	struct hash<vector<T>> {
	public:
		size_t operator()(const vector<T> &vec) const {
			size_t seed = vec.size();
			for (const T &x : vec) {
				seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};
}

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

ll time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
struct HawkKeyPair {
	size_t logn, n;
	int num_occuring = 0; // # of times, we find this in the text file.

	vector<int8_t> f, g, F, G;
	vector<fpr> q00, q01, q11;

	bool read(ifstream &in, size_t logn);
	void complete(uint8_t *tmp);
};

unordered_map< pair<vector<int8_t>, vector<int8_t>>, int> lookup_fg;
vector<HawkKeyPair> bases;

void HawkKeyPair::complete(uint8_t *tmp)
{
	F.resize(n);
	G.resize(n);
	q00.resize(n);
	q01.resize(n);
	q11.resize(n);

	assert(Zf(complete_private)(&f[0], &g[0], &F[0], &G[0], &q00[0], &q01[0], &q11[0], logn, tmp) != 0);
}

bool HawkKeyPair::read(ifstream &in, size_t _logn)
{
	logn = _logn;
	n = MKN(logn);

	f.resize(n);
	g.resize(n);

	// read f,g and complete basis
	for (int8_t &x : f) {
		int _x; in >> _x; x = (int8_t) _x;
	}
	for (int8_t &x : g) {
		int _x; in >> _x; x = (int8_t) _x;
	}
	return !!in;
}

void output_poly(int8_t *f, size_t logn) {
	for (size_t u = 0; u < MKN(logn); u++)
		printf("%d ", f[u]);
	printf("\n");
}

#define SQR(x) ((x) * (x))

int approximate_fail_prob(inner_shake256_context *rng,
	int8_t *f, int8_t *g, int8_t *F, int8_t *G,
	fpr *q00, fpr *q01, fpr *q11,
	unsigned logn, uint8_t *tmp)
{
	int16_t sig[512];
	uint8_t h[512/4];

	size_t n, u;
	n = MKN(logn);
	uint32_t bound = (uint32_t)(SQR(1.1 * 2 * 1.278) * (2*n));

	fpr expkey[EXPANDED_SECKEY_SIZE(9)];
	Zf(expand_seckey)(expkey, f, g, F, logn);

	int fails = 0;
	for (u = 0; u < 1024; u++) {
		inner_shake256_extract(rng, h, n / 4);
		// Make sure that sign.c may fail on generating a signature that does not decompress correctly.
		Zf(sign)(rng, sig, expkey, h, logn, tmp);
		if (!Zf(verify)(h, sig, q00, q01, q11, logn, tmp))
			fails++;
	}
	return fails;
}


int main(int argc, char **argv) {
	// Initialize a RNG.
	inner_shake256_context sc;
	inner_shake256_init(&sc);
	struct timeval tv;
	gettimeofday(&tv, NULL);
	inner_shake256_inject(&sc, (unsigned char *)&tv, sizeof tv);
	inner_shake256_flip(&sc);

	if (argc <= 1) {
		printf("Usage: ./a.out file\n");
		return 1;
	}

	ifstream in(argv[1]);

	string lines[3];
	getline(in, lines[0]); // Seed: blabla
	getline(in, lines[1]); // Using X threads.
	getline(in, lines[2]); // logn = Y, #sims = Z

	size_t logn = lines[2][7] - '0'; // only works up to 9, so we'll manage up to this NIST-1.
	size_t n = MKN(logn);
	vector<uint8_t> tmp(50 * n);

	HawkKeyPair newkey;
	size_t counter = 0;
	while (newkey.read(in, logn)) {
		pair<vector<int8_t>, vector<int8_t>> fg(newkey.f, newkey.g);
		auto it = lookup_fg.find(fg);
		if (it != lookup_fg.end()) {
			int idx = it->second;
			bases[idx].num_occuring++;
		} else {
			newkey.complete(&tmp[0]);
			newkey.num_occuring = 1;
			lookup_fg[fg] = bases.size();
			bases.push_back(newkey);
		}
		counter++;
	}

	// Make a selection of the bad bases that we want to analyse:
	if (argc > 3 && strcmp(argv[2], "--select") == 0) {
		int threshold = atoi(argv[3]);
		for (int i = 0; i < 3; i++) printf("%s\n", lines[i].c_str());

		for (HawkKeyPair &key : bases) {
			if (key.num_occuring >= threshold) {
				for (size_t u = 0; u < n; u++) printf("%d ", key.f[u]);
				printf("\n");
				for (size_t u = 0; u < n; u++) printf("%d ", key.g[u]);
				printf("\n\n");
			}
		}
		return 0;
	}

	// printf("Input is read: %d unique bases out of %d fails.\n", (int) bases.size(), (int) counter);

	fpr *inv_q00 = (fpr *)&tmp[0];

	// printf("The following are all bad bases:\n");
	for (HawkKeyPair &key : bases) {
		for (size_t u = 0; u < n/2; u++) {
			inv_q00[u] = fpr_inv(key.q00[u]);
			inv_q00[u + n/2] = fpr_zero;
		}

		Zf(iFFT)(inv_q00, logn);
		// TODO: compare (1/q00)[0] to 1/(q00[0])....

		double cst_term = inv_q00[0].v, norm2 = 0.0;
		for (size_t u = 0; u < n; u++)
			norm2 += fpr_sqr(inv_q00[u]).v;
		norm2 = sqrt(norm2);
		printf("(f*,g*)/(ff*+gg*) coeff^2: %.8f vs || 1/(ff* + gg*) ||: %.8f\n", cst_term, norm2);
		// printf("%.8f %d\n", cst_term, key.num_occuring);
	}
	// exit(0);

	// Now compare this to the average private key
	int8_t f[512], g[512], F[512], G[512];
	fpr q00[512], q01[512], q11[512];

	// printf("On average: ");
	const int n_repetitions = 1000;

	fpr avg = fpr_zero, var = fpr_zero;
	double avgN = 0.0, varN = 0.0;
	for (int i = 0; i < n_repetitions; i++) {
		// Generate key pair.
		Zf(keygen)(&sc, f, g, F, G, q00, q01, q11, logn, &tmp[0]);

		for (size_t u = 0; u < n/2; u++) {
			q00[u] = fpr_inv(q00[u]);
			q00[u + n/2] = fpr_zero;
		}
		Zf(iFFT)(q00, logn);

		fpr cstQ00 = q00[0];

		// Zf(FFT)(q00, logn); // restore value.
		// for (size_t u = 0; u < n/2; u++) q00[u] = fpr_inv(q00[u]);
		// printf("%.8f %d\n", cstQ00.v, approximate_fail_prob(&sc, f, g, F, G, q00, q01, q11, logn, &tmp[0]));

		double norm2 = 0.0;
		for (size_t u = 0; u < n; u++) norm2 += fpr_sqr(q00[u]).v;
		norm2 = sqrt(norm2);

		avg = fpr_add(avg, cstQ00);
		var = fpr_add(var, fpr_sqr(cstQ00));
		avgN += norm2;
		varN += norm2 * norm2;
	}
	// printf("\n");

	// Expectation: avg ~ sqrt(2n) sigma_pk / (2n sigma_pk^2)
	//	= 1 / (sqrt(2n) sigma_pk).
	avg = fpr_div(avg, fpr_of(n_repetitions));
	var = fpr_sub(fpr_div(var, fpr_of(n_repetitions)), fpr_sqr(avg));

	avgN /= n_repetitions;
	varN = varN / n_repetitions - avgN * avgN;

	printf("sum of coefficients^2 of (f*, g*) / q00 = %.8f +/- %.8f\n",
		avg.v, sqrt(var.v));
	printf("||1/q00|| = %.8f +/- %.8f\n", avgN, sqrt(varN));
	return 0;
}
