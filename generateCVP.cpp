// g++ build/compress.o build/ffo.o build/fft.o build/fpr.o build/keygen.o build/rng.o build/sampler.o build/shake.o build/sign.o build/vrfy.o generateCVP.cpp -lgmp -lfplll
#include <assert.h>
#include <stdio.h>
// x86_64 specific:
#include <sys/time.h>

#include <fplll.h>
#include <vector>

extern "C" {
	#define restrict
	#include "inner.h"
}

using namespace std;
using namespace fplll;

// Simple randomness generator:
void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x)
		*x = ((unsigned char) rand());
}

// =============================================================================
// | TESTING CODE                                                              |
// =============================================================================
int valid_sigma(int logn, fpr sigma_sig) {
	return !fpr_lt(sigma_sig, fpr_sigma_min[logn])
		&& fpr_lt(sigma_sig, fpr_div(fpr_of(18205), fpr_of(10000)));
}

void output_poly(int8_t *p, int logn) {
	for (size_t u = 0; u < MKN(logn); u++) {
		if (u) printf(" ");
		printf("%d", p[u]);
	}
	printf("\n");
}

void norm2(int8_t *f, int8_t *g, size_t logn) {
	unsigned res = 0;
	for (size_t u = 0; u < MKN(logn); u++) {
		res += (unsigned) f[u] * f[u];
		res += (unsigned) g[u] * g[u];
	}
	printf("Norm: %u\n", res);
}

void report_bytes(fpr *q00, fpr *q10, size_t logn) {
	int16_t Q0[512], Q1[512];
	size_t n = MKN(logn);

	Zf(iFFT)(q00, logn);
	Zf(iFFT)(q10, logn);
	for (size_t u = 0; u < n; u++) {
		Q0[u] = fpr_rint(q00[u]);
		Q1[u] = fpr_rint(q10[u]);
	}
	Zf(FFT)(q00, logn);
	Zf(FFT)(q10, logn);

	for (size_t u = 0; u < n; u++)
		printf("%d,", Q1[u]);
	printf("\n");
	size_t nr_bytes = Zf(encode_pubkey)(nullptr, 0, Q0, Q1, logn);
	printf("#bytes(pubkey) = %u\n", nr_bytes);
}

static void get_gscoords(const Matrix<FP_NR<mpfr_t>> &matrix, const Matrix<FP_NR<mpfr_t>> &mu,
                         const Matrix<FP_NR<mpfr_t>> &r, const vector<FP_NR<mpfr_t>> &v,
                         vector<FP_NR<mpfr_t>> &vcoord)
{

  int n = matrix.get_rows(), m = matrix.get_cols();

  if (static_cast<int>(vcoord.size()) != n)
    vcoord.resize(n);
  FPLLL_DEBUG_CHECK(mu.get_rows() == n && mu.get_cols() == n && r.get_rows() == n &&
                    r.get_cols() == n && static_cast<int>(v.size()) == m);

  for (int i = 0; i < n; i++)
  {
    vcoord[i] = 0.0;
    for (int j = 0; j < m; j++)
      vcoord[i].addmul(v[j], matrix(i, j));
    for (int j = 0; j < i; j++)
      vcoord[i].submul(mu(i, j), vcoord[j]);
  }
  for (int i = 0; i < n; i++)
  {
    vcoord[i].div(vcoord[i], r(i, i));
  }
}


static void babai(const FP_mat<mpfr_t> &matrix, const Matrix<FP_NR<mpfr_t>> &mu,
                  const Matrix<FP_NR<mpfr_t>> &r, const vector<FP_NR<mpfr_t>> &target,
                  vector<FP_NR<mpfr_t>> &target_coord)
{

  int d = matrix.get_rows();
  get_gscoords(matrix, mu, r, target, target_coord);
  for (int i = d - 1; i >= 0; i--)
  {
    target_coord[i].rnd(target_coord[i]);
    for (int j = 0; j < i; j++)
      target_coord[j].submul(mu(i, j), target_coord[i]);
  }
}

int quick_CVP(ZZ_mat<mpz_t> &b, const vector<Z_NR<mpz_t>> &int_target, vector<Z_NR<mpz_t>> &sol_coord)
{
	// status = fplll::closest_vector(m, target, sol_coord, CVPM_PROVED, CVP_DEFAULT);

  // d = lattice dimension (note that it might decrease during preprocessing)
  int d = b.get_rows();
  // n = dimension of the space
  int n = b.get_cols();

  FPLLL_CHECK(d > 0 && n > 0, "closestVector: empty matrix");
  FPLLL_CHECK(d <= n, "closestVector: number of vectors > size of the vectors");

  // Sets the floating-point precision
  // Error bounds on GSO are valid if prec >= minprec
  double rho;
  int min_prec = gso_min_prec(rho, d, LLL_DEF_DELTA, LLL_DEF_ETA);
  int prec     = max(53, min_prec + 10);
  int old_prec = FP_NR<mpfr_t>::set_prec(prec);

  // Allocates space for vectors and matrices in constructors
  ZZ_mat<mpz_t> empty_mat;
  MatGSO<Z_NR<mpz_t>, FP_NR<mpfr_t>> gso(b, empty_mat, empty_mat, GSO_INT_GRAM);
  vector<FP_NR<mpfr_t>> target_coord;
  FP_NR<mpfr_t> max_dist;
  Z_NR<mpz_t> itmp1;

  // Computes the Gram-Schmidt orthogonalization in floating-point
  gso.update_gso();
  gen_zero_vect(sol_coord, d);

  /* Applies Babai's algorithm. Because we use fp, it might be necessary to
      do it several times (if ||target|| >> ||b_i||) */
  FP_mat<mpfr_t> float_matrix(d, n);
  vector<FP_NR<mpfr_t>> target(n), babai_sol;
  vector<Z_NR<mpz_t>> int_new_target = int_target;

  for (int i = 0; i < d; i++)
    for (int j = 0; j < n; j++)
      float_matrix(i, j).set_z(b(i, j));

  for (int loop_idx = 0;; loop_idx++) {
    if (loop_idx >= 0x100 && ((loop_idx & (loop_idx - 1)) == 0))
      FPLLL_INFO("warning: possible infinite loop in Babai's algorithm");

    for (int i = 0; i < n; i++) {
      target[i].set_z(int_new_target[i]);
    }
    babai(float_matrix, gso.get_mu_matrix(), gso.get_r_matrix(), target, babai_sol);
    int idx;
    for (idx = 0; idx < d && babai_sol[idx] >= -1 && babai_sol[idx] <= 1; idx++) { }
    if (idx == d)
      break;

    for (int i = 0; i < d; i++) {
      itmp1.set_f(babai_sol[i]);
      sol_coord[i].add(sol_coord[i], itmp1);
      for (int j = 0; j < n; j++)
        int_new_target[j].submul(itmp1, b(i, j));
    }
  }
  // FPLLL_TRACE("BabaiSol=" << sol_coord);
  get_gscoords(float_matrix, gso.get_mu_matrix(), gso.get_r_matrix(), target, target_coord);

  /* Computes a very large bound to make the algorithm work
      until the first solution is found */
  max_dist = 0.0;
  for (int i = 0; i < d; i++) {
    // get_r_exp(i, i) = r(i, i) because gso is initialized without GSO_ROW_EXPO
    max_dist.add(max_dist, gso.get_r_exp(i, i));
  }

  vector<int> max_indices(d, 1);
  FastErrorBoundedEvaluator evaluator(n, gso.get_mu_matrix(), gso.get_r_matrix(), EVALMODE_CV);

  printf("Start of enumeration.\n");
  fflush(stdout);
  // Main loop of the enumeration
  Enumeration<Z_NR<mpz_t>, FP_NR<mpfr_t>> enumobj(gso, evaluator, max_indices);
  enumobj.enumerate(0, d, max_dist, 0, target_coord);

  int result = RED_ENUM_FAILURE;
  if (!evaluator.empty()) {
    for (int i = 0; i < d; i++) {
      itmp1.set_f(evaluator.begin()->second[i]);
      sol_coord[i].add(sol_coord[i], itmp1);
    }
    result = RED_SUCCESS;
  }

  FP_NR<mpfr_t>::set_prec(old_prec);
  return result;
}

int main(int argc, char **argv) {
	uint8_t b[48 * 512];
	int8_t f[512], g[512], F[512], G[512];
	fpr q00[512], q10[512], q11[512];
	unsigned char shakeseed[48];
	inner_shake256_context sc;

	if (argc != 2) {
		printf("Usage: %s logn\n", argv[0]);
		exit(1);
	}
	size_t logn = atoi(argv[1]);
	if (logn < 1 || logn > 9) {
		printf("Error: %s only works with 1 <= logn <= 9\n", argv[0]);
		exit(1);
	}

	// set seed
	struct timeval tv;
	gettimeofday(&tv, NULL);
	int seed = 1000000 * tv.tv_sec + tv.tv_usec;
	printf("Seed: %d\n", seed);
	srand(seed);

	const fpr sigma_kg = fpr_div(fpr_of(1425), fpr_of(1000));
	assert(valid_sigma(logn, sigma_kg));

	// Initialize a RNG.
	randombytes(shakeseed, sizeof shakeseed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, shakeseed, sizeof shakeseed);
	inner_shake256_flip(&sc);

	// Generate key pair.
	Zf(keygen)(&sc, f, g, F, G, q00, q10, q11, fpr_inv(sigma_kg), logn, b);

	// Output secret key (f, g, F, G)
	output_poly(F, logn);
	output_poly(G, logn);

	norm2(f, g, logn);
	norm2(F, G, logn);
	report_bytes(q00, q10, logn);

	size_t n = MKN(logn);
	ZZ_mat<mpz_t> m(n, 2*n);
	vector<Z_NR<mpz_t>> target(2*n), sol_coord(2*n);

	for (size_t u = 0; u < n; u ++) {
		// (f,g) x^u
		for (size_t v = 0; v < u; v++) {
			m(u, v) = -f[n-u+v];
			m(u, n + v) = -g[n-u+v];
		}
		for (size_t v = u; v < n; v++) {
			m(u, v) = f[v-u];
			m(u, n + v) = g[v-u];
		}
	}
	for (size_t v = 0; v < n; v++) {
		target[v] = F[v];
		target[n + v] = G[v];
	}

const char *const RED_STATUS_STR[RED_STATUS_MAX] = {
    "success",
    "",
    "infinite number in GSO",
    "infinite loop in babai",
    "infinite loop in LLL",
    "error in SVP solver",
    "error in BKZ",
    "time limit exceeded in BKZ",
    "loops limit exceeded in BKZ",
    "error in HLLL",
    "increase of the norm",
    "error in weak size reduction",
    "Please see https://github.com/fplll/fplll/wiki/fplll-errors-FAQ for more information."
};


	int status;
	status = lll_reduction(m, LLL_DEF_DELTA, LLL_DEF_ETA, LM_PROVED, FT_DEFAULT, 0, LLL_DEFAULT);
	if (status != RED_SUCCESS) {
		cerr << "Status LLL: " << RED_STATUS_STR[status] << endl;
		return status;
	}

	status = quick_CVP(m, target, sol_coord);
	// status = fplll::closest_vector(m, target, sol_coord, CVPM_PROVED, CVP_DEFAULT);
	if (status != RED_SUCCESS) {
		cerr << "Status CVP: " << RED_STATUS_STR[status] << endl;
		return status;
	}

	for (size_t v = 0; v < n; v++) {
		cout << sol_coord[v] << ",";
	}
	cout << endl;

	int8_t newF[512], newG[512];
	memcpy(newF, F, sizeof newF);
	memcpy(newG, G, sizeof newG);

	fpr fft_k[512];
	for (size_t v = 0; v < n; v++) {
		fft_k[v] = fpr_of(sol_coord[v].get_si());
	}
	Zf(FFT)(fft_k, logn);
	Zf(poly_mul_autoadj_fft)(fft_k, q00, logn);
	Zf(poly_sub)(q10, fft_k, logn);

	for (size_t v = 0; v < n; v++) {
		int x = sol_coord[v].get_si();
		if (x == 0) continue;
		// subtract the CVP from (F,G)
		for (size_t i = 0; i < n; i++) {
			newF[i] -= x * m(v, i).get_si();
			newG[i] -= x * m(v, i + n).get_si();
		}
	}

	report_bytes(q00, q10, logn);

	output_poly(newF, logn);
	output_poly(newG, logn);

	norm2(newF, newG, logn);
	return 0;
}
