#include <cassert>
#include <climits>
#include <cstdio>
#include <vector>

#include <mutex>
#include <thread>

extern "C" {
	#define restrict
	#include "../inner.h"
}
using namespace std;

constexpr int MAXN = 1024;
constexpr int n_repetitions = 1 << 14;
constexpr int spread_keygen = 1 << 6;

double max_sizes[4][11][11] = {};

void FPC_MUL(double &a, double &b, double c, double d, double e, double f) {
	double x = c * e - d * f;
	double y = c * f + d * e;
	a = x;
	b = y;
}

void _FFT(double *f, unsigned logn, int type)
{
	unsigned u;
	size_t t, n, hn, m;

	n = (size_t)1 << logn;
	hn = n >> 1;
	t = hn;

	for (size_t i = 0; i < n; i++) {
		double x = fabs(f[i]);
		if (x > max_sizes[type][logn][0])
			max_sizes[type][logn][0] = x;
	}

	for (u = 1, m = 2; u < logn; u ++, m <<= 1) {
		size_t ht, hm, i1, j1;

		ht = t >> 1;
		hm = m >> 1;
		for (i1 = 0, j1 = 0; i1 < hm; i1 ++, j1 += t) {
			size_t j, j2;

			j2 = j1 + ht;
			double s_re, s_im;

			s_re = fpr_gm_tab[((m + i1) << 1) + 0].v;
			s_im = fpr_gm_tab[((m + i1) << 1) + 1].v;
			for (j = j1; j < j2; j ++) {
				double x_re, x_im, y_re, y_im;

				x_re = f[j];
				x_im = f[j + hn];
				y_re = f[j + ht];
				y_im = f[j + ht + hn];
				FPC_MUL(y_re, y_im, y_re, y_im, s_re, s_im);

				f[j] = x_re + y_re;
				f[j + hn] = x_im + y_im;

				f[j + ht] = x_re - y_re;
				f[j + ht + hn] = x_im - y_im;
			}
		}
		t = ht;

		for (size_t i = 0; i < n; i++) {
			double x = fabs(f[i]);
			if (x > max_sizes[type][logn][u])
				max_sizes[type][logn][u] = x;
		}

		for (j1 = 0; j1 < n; j1++) {
			f[j1] /= 2;
		}
	}
}

void _iFFT(double *f, unsigned logn)
{
	size_t u, n, hn, t, m;

	n = (size_t)1 << logn;
	t = 1;
	m = n;
	hn = n >> 1;

	for (size_t i = 0; i < n; i++) {
		double x = fabs(f[i]);
		if (x > max_sizes[3][logn][logn - 1])
			max_sizes[3][logn][logn - 1] = x;
	}

	for (u = logn; u > 1; u --) {
		size_t hm, dt, i1, j1;

		hm = m >> 1;
		dt = t << 1;
		for (i1 = 0, j1 = 0; j1 < hn; i1 ++, j1 += dt) {
			size_t j, j2;

			j2 = j1 + t;
			double s_re, s_im;

			s_re = fpr_gm_tab[((hm + i1) << 1)+0].v;
			s_im =-fpr_gm_tab[((hm + i1) << 1)+1].v;
			for (j = j1; j < j2; j ++) {
				double x_re, x_im, y_re, y_im;

				x_re = f[j];
				x_im = f[j + hn];
				y_re = f[j + t];
				y_im = f[j + t + hn];

				f[j] = x_re + y_re;
				f[j + hn] = x_im + y_im;

				x_re -= y_re;
				x_im -= y_im;
				FPC_MUL(f[j + t], f[j + t + hn], x_re, x_im, s_re, s_im);
			}
		}
		t = dt;
		m = hm;

		for (size_t i = 0; i < n; i++) {
			double x = fabs(f[i]);
			if (x > max_sizes[3][logn][u - 2])
				max_sizes[3][logn][u - 2] = x;
		}
	}
}

void measure_sizes(unsigned logn) {
	uint8_t b[48 * MAXN], h[MAXN / 4];
	int8_t f[MAXN], g[MAXN], F[MAXN], G[MAXN];
	int16_t iq00[MAXN], iq10[MAXN], s1[MAXN];
	size_t n, hn;
	double t0[MAXN], t1[MAXN], t2[MAXN];
	unsigned char seed[48];
	inner_shake256_context sc;

	n = MKN(logn);
	hn = n >> 1;

	// Initialize a RNG.
	Zf(get_seed)(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

	for (int rep = 0; rep < n_repetitions; rep++) {
		if (rep % spread_keygen == 0) {
			// Generate key pair.
			Zf(keygen)(&sc, f, g, F, G, iq00, iq10, logn, b);
		}

		// make a signature of a random message...
		inner_shake256_extract(&sc, h, sizeof h);

		// Compute the signature.
		while (!Zf(sign_dyn)(&sc, s1, f, g, F, G, h, logn, b)) {}

		for (size_t u = 0; u < n; u++) {
			t0[u] = (double) (((h[u>>3] >> (u & 7)) & 1) - 2 * s1[u]);
			t1[u] = (double) iq10[u];
			t2[u] = (double) iq00[u];
		}
		t2[0] = 0;

		_FFT(t0, logn, 0);
		_FFT(t1, logn, 1);
		_FFT(t2, logn, 2);
		for (size_t u = 0; u < hn; u++) {
			t2[u] += (double)iq00[0] / (1 << (logn - 1));
		}

		for (size_t u = 0; u < hn; u++) {
			double re = t0[u] * t1[u] - t0[u + hn] * t1[u + hn];
			double im = t0[u] * t1[u + hn] + t0[u + hn] * t1[u];
			t0[u] = re / t2[u];
			t0[u + hn] = im / t2[u];
		}

		for (size_t u = 0; u < n; u++) {
			printf("%.3f ", t0[u]);
		}
		printf("\n");

		_iFFT(t0, logn);
	}

	printf("logn = %d:\n", logn);
	for (int type = 0; type < 4; type++) {
		if (type == 0) printf("s0:  ");
		if (type == 1) printf("q10: ");
		if (type == 2) printf("q00: ");
		if (type == 3) printf("inv: ");

		double worst = 0.0;
		for (size_t u = 0; u < logn; u++) {
			printf("%.8f ", max_sizes[type][logn][u]);
			worst = max(worst, max_sizes[type][logn][u]);
		}
		printf(" -> %.8f (%d bits)\n", worst, (int) ceil(log2(worst)));
	}
	printf("\n");
}

int main() {
	for (size_t logn = 1; logn <= 10; logn++) {
		measure_sizes(logn);
	}
	return 0;
}
