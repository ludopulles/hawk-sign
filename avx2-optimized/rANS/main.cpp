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

#include "platform.h"
#include "rans_byte.h"

// This is just the sample program. All the meat is in rans_byte.h.

// ---- Stats

struct SymbolStats
{
    uint32_t freqs[256];
    uint32_t cum_freqs[257];

    void init_freqs();
    void normalize_freqs(uint32_t target_total);
};

void SymbolStats::init_freqs()
{
    for (int i=0; i < 256; i++)
        freqs[i] = 0;
}

void SymbolStats::normalize_freqs(uint32_t target_total)
{
    assert(target_total >= 256);

    cum_freqs[0] = 0;
    for (int i=0; i < 256; i++)
        cum_freqs[i+1] = cum_freqs[i] + freqs[i];
    uint32_t cur_total = cum_freqs[256];

    // resample distribution based on cumulative freqs
    for (int i = 1; i <= 256; i++)
        cum_freqs[i] = ((uint64_t)target_total * cum_freqs[i])/cur_total;

    // if we nuked any non-0 frequency symbol to 0, we need to steal
    // the range to make the frequency nonzero from elsewhere.
    //
    // this is not at all optimal, i'm just doing the first thing that comes to mind.
    for (int i=0; i < 256; i++) {
        if (freqs[i] && cum_freqs[i+1] == cum_freqs[i]) {
            // symbol i was set to zero freq

            // find best symbol to steal frequency from (try to steal from low-freq ones)
            uint32_t best_freq = ~0u;
            int best_steal = -1;
            for (int j=0; j < 256; j++) {
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
    assert(cum_freqs[0] == 0 && cum_freqs[256] == target_total);
    for (int i=0; i < 256; i++) {
        if (freqs[i] == 0)
            assert(cum_freqs[i+1] == cum_freqs[i]);
        else
            assert(cum_freqs[i+1] > cum_freqs[i]);

        // calc updated freq
        freqs[i] = cum_freqs[i+1] - cum_freqs[i];
    }
}

void report_ANS_size(unsigned logn, const vector<vector<int>> &data) {

    static const uint32_t prob_bits = 20;
    static const uint32_t prob_scale = 1 << prob_bits;

    SymbolStats stats;
    stats.init_freqs();
	for (const vector<int> &v : data) for (const int &x : v) stats.freqs[x]++;
    stats.normalize_freqs(prob_scale);

    static size_t out_max_size = 32<<20; // 32MB
    uint8_t* out_buf = new uint8_t[out_max_size];

    // try rANS encode
    uint8_t *rans_begin;
    RansEncSymbol esyms[256];
    RansDecSymbol dsyms[256];

    for (int i=0; i < 256; i++) {
        RansEncSymbolInit(&esyms[i], stats.cum_freqs[i], stats.freqs[i], prob_bits);
        RansDecSymbolInit(&dsyms[i], stats.cum_freqs[i], stats.freqs[i]);
    }

    // ---- regular rANS encode/decode. Typical usage.
	int sumsz = 0, squsz = 0;

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

		int head_bytes = (int) (out_buf + out_max_size - rans_begin);
		// printf("rANS: %d bytes\n", head_bytes);

		int x = (((S1_LOBITS(logn) + 1) << logn) + 7) / 8;

		x += head_bytes;
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

	vector<vector<int>> data(num_sig, vector<int>(1 << logn));
	for (int i = 0; i < num_sig; i++) {
		for (int u = 0; u < (1<<logn); u++) {
			uint16_t w = (uint16_t) sigs[i][u];
			w ^= -(w >> 15);
			data[i][u] = w >> S1_LOBITS(logn);
			assert(data[i][u] < 256);
		}
	}

	report_ANS_size(logn, data);
	return 0;
}
