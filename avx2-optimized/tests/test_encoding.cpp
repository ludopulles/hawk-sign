/*
 * Report on the size of signatures
 */
#include <cassert>
#include <climits>
#include <cstdio>
#include <iostream>
#include <vector>
#include <mutex>
#include <thread>
#include <vector>

// x86_64 specific:
#include<sys/time.h>
using namespace std;

extern "C" {
	#define restrict
	#include "../inner.h"
}

#define MAXLOGN (10)
#define MAXN MKN(MAXLOGN)



#define S1_LOBITS(logn) ((logn) == 10 ? 6u : 5u)
#define BOUND_S1(logn) ((logn) == 10 ? 1024 : 512)

size_t
encode_sig(void *out, size_t max_out_len, const int16_t *s1, unsigned logn,
	size_t lo_bits)
{
	uint8_t *buf;
	size_t n, u, v;
	uint64_t acc;
	unsigned acc_len;

	n = MKN(logn);
	buf = (uint8_t *)out;

	/*
	 * Make sure no coefficient is too large.
	 */
	for (u = 0; u < n; u ++) {
		if (s1[u] < -BOUND_S1(logn) || s1[u] >= BOUND_S1(logn)) return 0;
	}

	acc = 0;
	acc_len = 0;
	v = 0;
	for (u = 0; u < n; u ++) {
		uint16_t w;

		/*
		 * Push sign bit
		 */
		w = (uint16_t) s1[u];
		acc = (acc << 1) | (w >> 15);
		w ^= -(w >> 15);

		/*
		 * Push the lowest `lo_bits` bits of w which is equal to |x| - [x < 0].
		 */
		acc <<= lo_bits;
		acc |= w & ((1U << lo_bits) - 1);
		w >>= lo_bits;
		acc_len += lo_bits + 1;

		acc <<= (w + 1);
		acc |= 1;
		acc_len += w + 1;

		assert(acc_len < 64);

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


int head_fs[256] = {};

int main(int argc, char **argv) {
	if (argc != 2) {
		printf("Usage: ./file logn\n");
		return 1;
	}
	unsigned logn = atoi(argv[1]);
	assert(1 <= logn && logn <= MAXLOGN);

	size_t num_sig;
	cin >> num_sig;

	vector<vector<int16_t>> sigs(num_sig, vector<int16_t>(1 << logn));
	for (vector<int16_t> &v : sigs)
		for (int16_t &x : v) cin >> x;

	int sumsz = 0, squsz = 0;
	for (int i = 0; i < num_sig; i++) {
		int x = encode_sig(NULL, 0, &*sigs[i].begin(), logn, S1_LOBITS(logn));
		sumsz += x;
		squsz += x*x;
	}

	double avg = double(sumsz) / num_sig;
	double var = double(squsz) / num_sig - avg * avg;

	printf("%.3f +/- %.3f\n", avg, sqrt(var));

	for (int i = 0; i < num_sig; i++) {
		for (int16_t x : sigs[i]) {
			uint16_t w = (uint16_t) x;
			w ^= -(w >> 15);
			assert((w >> S1_LOBITS(logn)) < 256);
			head_fs[w >> S1_LOBITS(logn)]++;
		}
	}

	for (size_t u = 0; u < 256; u++) {
		if (head_fs[u]) printf("%d: %d\n", u, head_fs[u]);
	}
	printf("\n");

	return 0;
}
