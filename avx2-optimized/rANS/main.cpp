/*
 * Report on the size of signatures
 */
#include<bits/stdc++.h>

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

static const size_t low_bits_q00[11] = {
	0 /* unused */, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6
};
static const size_t low_bits_q01[11] = {
	0 /* unused */, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
};

size_t encode_q00(const int16_t *q00, unsigned logn)
{
	size_t n, u, v;
	uint16_t w;
	int16_t bound;
	unsigned acc_len;

	n = MKN(logn);
	bound = (int16_t)(1U << Zf(bits_q00)[logn]);

	/*
	 * Make sure no coefficient is too large.
	 */
	for (u = 1; u < n/2; u ++) {
		if (q00[u] < -bound || q00[u] >= bound) {
			return 0;
		}
	}

	acc_len = 0;
	v = 2;
	for (u = 1; u < n/2; u ++) {
		w = (uint16_t) q00[u];
		w ^= -(w >> 15);

		w >>= low_bits_q00[logn];
		acc_len += low_bits_q00[logn] + 1;
		acc_len += w + 1;

		while (acc_len >= 8) {
			acc_len -= 8;
			v ++;
		}
	}
	if (acc_len > 0) v ++;
	return v;
}

size_t encode_q01(const int16_t *q01, unsigned logn)
{
	size_t n, u, v;
	uint16_t w;
	int16_t bound;
	unsigned acc_len;

	n = MKN(logn);
	bound = (int16_t)(1U << Zf(bits_q01)[logn]);

	for (u = 1; u < n/2; u ++) {
		if (q01[u] < -bound || q01[u] >= bound) {
			return 0;
		}
	}

	acc_len = 0;
	v = 0;
	for (u = 0; u < n; u ++) {
		w = (uint16_t) q01[u];
		w ^= -(w >> 15);

		w >>= low_bits_q01[logn];
		acc_len += low_bits_q01[logn] + 1;
		acc_len += w + 1;

		while (acc_len >= 8) {
			acc_len -= 8;
			v ++;
		}
	}

	if (acc_len > 0) v ++;
	return v;
}


#include "platform.h"
#include "rans_byte.h"

// This is just the sample program. All the meat is in rans_byte.h.

// ---- Stats
#define NUM_FREQ (16 * 1024)

struct SymbolStats
{
    uint32_t freqs[NUM_FREQ];
    uint32_t cum_freqs[NUM_FREQ + 1];

    void init_freqs();
    void normalize_freqs(uint32_t target_total);
};

void SymbolStats::init_freqs()
{
    for (int i=0; i < NUM_FREQ; i++)
        freqs[i] = 0;
}

void SymbolStats::normalize_freqs(uint32_t target_total)
{
    assert(target_total >= NUM_FREQ);

    cum_freqs[0] = 0;
    for (int i=0; i < NUM_FREQ; i++)
        cum_freqs[i+1] = cum_freqs[i] + freqs[i];
    uint32_t cur_total = cum_freqs[NUM_FREQ];

    // resample distribution based on cumulative freqs
    for (int i = 1; i <= NUM_FREQ; i++)
        cum_freqs[i] = ((uint64_t)target_total * cum_freqs[i])/cur_total;

    // if we nuked any non-0 frequency symbol to 0, we need to steal
    // the range to make the frequency nonzero from elsewhere.
    //
    // this is not at all optimal, i'm just doing the first thing that comes to mind.
    for (int i=0; i < NUM_FREQ; i++) {
        if (freqs[i] && cum_freqs[i+1] == cum_freqs[i]) {
            // symbol i was set to zero freq

            // find best symbol to steal frequency from (try to steal from low-freq ones)
            uint32_t best_freq = ~0u;
            int best_steal = -1;
            for (int j=0; j < NUM_FREQ; j++) {
                uint32_t freq = cum_freqs[j+1] - cum_freqs[j];
                if (freq > 1 && freq < best_freq) {
                    best_freq = freq;
                    best_steal = j;
                }
            }
            assert(best_steal != -1);

            // and steal from it!
            if (best_steal < i) {
                for (int j = best_steal + 1; j <= i; j++)
                    cum_freqs[j]--;
            } else {
                assert(best_steal > i);
                for (int j = i + 1; j <= best_steal; j++)
                    cum_freqs[j]++;
            }
        }
    }

    // calculate updated freqs and make sure we didn't screw anything up
    assert(cum_freqs[0] == 0 && cum_freqs[NUM_FREQ] == target_total);
    for (int i=0; i < NUM_FREQ; i++) {
        if (freqs[i] == 0)
            assert(cum_freqs[i+1] == cum_freqs[i]);
        else
            assert(cum_freqs[i+1] > cum_freqs[i]);

        // calc updated freq
        freqs[i] = cum_freqs[i+1] - cum_freqs[i];
    }
}

void report_ANS_size(const vector<vector<int>> &data, int othersize) {
    static const uint32_t prob_bits = 16;
    static const uint32_t prob_scale = 1 << prob_bits;

    SymbolStats stats;
    stats.init_freqs();
	for (const vector<int> &v : data) {
		for (const int &x : v) {
			assert(0 <= x && x < NUM_FREQ);
			stats.freqs[x]++;
		}
	}
    stats.normalize_freqs(prob_scale);

    static size_t out_max_size = 32<<20; // 32MB
    uint8_t* out_buf = new uint8_t[out_max_size];

    // try rANS encode
    uint8_t *rans_begin;
    RansEncSymbol esyms[NUM_FREQ];
    RansDecSymbol dsyms[NUM_FREQ];

    for (int i=0; i < NUM_FREQ; i++) {
        RansEncSymbolInit(&esyms[i], stats.cum_freqs[i], stats.freqs[i], prob_bits);
        RansDecSymbolInit(&dsyms[i], stats.cum_freqs[i], stats.freqs[i]);
    }

    // ---- regular rANS encode/decode. Typical usage.
	long long sumsz = 0, squsz = 0;

    printf("rANS encode:\n");

    for (const vector<int> &v : data) {
        RansState rans;
        RansEncInit(&rans);

        uint8_t* ptr = out_buf + out_max_size; // *end* of output buffer
        for (size_t i = v.size(); i > 0; ) { // NB: working in reverse!
            int s = v[--i];
            RansEncPutSymbol(&rans, &ptr, &esyms[s]);
        }
        RansEncFlush(&rans, &ptr);
        rans_begin = ptr;

		int x = othersize + (int) (out_buf + out_max_size - rans_begin);
		// printf("rANS: %d bytes (%d total)\n", head_bytes, x);

		sumsz += x;
		squsz += x*x;
	}

    delete[] out_buf;

	int num_sig = data.size();
	double avg = double(sumsz) / num_sig;
	double var = double(squsz) / num_sig - avg * avg;
	printf("%.3f +/- %.3f\n", avg, sqrt(var));

}

int main(int argc, char **argv) {
	if (argc != 3) {
		printf("Usage: ./file logn nsamples\n");
		return 1;
	}

	unsigned logn = atoi(argv[1]);
	assert(1 <= logn && logn <= MAXLOGN);
	int nsamples = atoi(argv[2]);
	assert(1 <= nsamples && nsamples <= 1000 * 1000);

	int n = MKN(logn);
	unsigned char seed[48];
	inner_shake256_context sc;

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	vector<vector<int16_t>> iq00(nsamples, vector<int16_t>(n)), iq01(nsamples, vector<int16_t>(n));
	for (int s = 0; s < nsamples; s++) {
		int8_t f[MAXN], g[MAXN], F[MAXN], G[MAXN];
		uint8_t tmp[44 * MAXN];

		Zf(keygen)(&sc, f, g, F, G, &*iq00[s].begin(), &*iq01[s].begin(), logn, tmp);
	}

	long long sumsz = 0, squsz = 0;
	for (int s = 0; s < nsamples; s++) {
		int x = encode_q00(&*iq00[s].begin(), logn);
		sumsz += x;
		squsz += x * x;
	}

	double avg = double(sumsz) / nsamples;
	double var = double(squsz) / nsamples - avg * avg;
	printf("Golomb-Rice(q00) ~ %.3f +/- %.3f\n", avg, sqrt(var));

	sumsz = 0; squsz = 0;
	for (int s = 0; s < nsamples; s++) {
		int x = encode_q01(&*iq01[s].begin(), logn);
		sumsz += x;
		squsz += x * x;
	}

	avg = double(sumsz) / nsamples;
	var = double(squsz) / nsamples - avg * avg;
	printf("Golomb-Rice(q01) ~ %.3f +/- %.3f\n", avg, sqrt(var));

	for (int saving = 0; saving <= (int) low_bits_q00[logn] && saving <= (int) low_bits_q01[logn]; saving++) {
		vector<vector<int>> dq00(nsamples, vector<int>(n/2 - 1)), dq01(nsamples, vector<int>(n));
		for (int s = 0; s < nsamples; s++) {
			for (int i = 1; i < n/2; i++) {
				uint16_t w = (uint16_t) iq00[s][i];
				w ^= -(w >> 15);
				dq00[s][i - 1] = w >> (low_bits_q00[logn] - saving);
			}
			for (int i = 0; i < n; i++) {
				uint16_t w = (uint16_t) iq01[s][i];
				w ^= -(w >> 15);
				dq01[s][i] = w >> (low_bits_q01[logn] - saving);
			}
		}

		report_ANS_size(dq00, 2 + ((low_bits_q00[logn] - saving + 1) * (n/2 - 1) + 7) / 8);
		report_ANS_size(dq01, ((low_bits_q01[logn] - saving + 1) * n + 7) / 8);
	}
	return 0;
}
