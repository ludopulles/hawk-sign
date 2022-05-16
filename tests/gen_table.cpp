/*
 * Code to be used for printing the precomputed tables for discrete Gaussian
 * sampling in keygen.c and sign.c.
 */
#include<bits/stdc++.h>
typedef long double FT;
const int LEN = 100;
const FT spk = 1.500, ssig = 1.278, p63 = powl(2, 63);
FT rho(FT x, FT sigma) { return expl(-(x * x) / (2.0 * sigma * sigma)); }

int main() {
	// Table for key generation:
	std::map<int, long double> table, csum;
	unsigned long long results[2][LEN] = {};

	for (int x = LEN; x >= -LEN; x--) {
		table[x] = rho(x, spk);
		csum[x] = csum[x + 1] + table[x];
	}

	results[0][0] = llroundl(p63 * table[0] / csum[-LEN]); // P(X == 0)
	for (int k = LEN; --k > 0; ) {
		// result[k] ~ P(X >= k+1 | X >= 1)
		results[0][k] = llroundl(p63 * csum[k + 1] / csum[1]);
	}

	int len = 0;
	while (results[0][len] != 0) len++;
	printf("static const uint64_t gauss_keygen[%d] = {\n", len);
	for (int x = 0; x < len; x++) {
		printf("\t%19lluu,\n", results[0][x]);
	}
	printf("};\n\n");

	// Table for signing:
	for (int coset = 0; coset < 2; coset++) {
		for (int x = LEN; x >= -LEN; x--) {
			table[x] = rho(x - 0.5 * coset, ssig);
			csum[x] = csum[x + 1] + table[x];
		}
		if (coset == 0) {
			// result[0] ~ P(X == 0)
			results[0][0] = llroundl(p63 * table[0] / csum[-LEN]);
		} else {
			// result[1] ~ (P(X == 1) + P(X == 0)) = P(X == 1 | X >= 1)
			results[1][1] = llroundl(p63 * table[1] / csum[1]);
		}

		for (int k = LEN; --k > coset; ) {
			// result[k] ~ P(X >= k + 1 | X >= 1 + coset)
			results[coset][k] = llroundl(p63 * csum[k + 1] / csum[1 + coset]);
		}
	}

	len = 0;
	while (results[0][len] != 0 || results[1][len + 1] != 0) len++;
	printf("static const uint64_t gauss_sign[%d] = {\n", 2 * len);
	for (int k = 0; k < len; k++)
		printf("\t%19lluu, %19lluu,\n", results[0][k], results[1][k + 1]);
	printf("};\n\n");

	/*
	 * Generate l2bounds using:
	 *     l2bound(logn) = floor( (verif_margin * 2 sigma_sig)^2 * 2n ),
	 * where verif_margin = 1.1.
	 */
	printf("const uint32_t Zf(l2bound)[10] = {\n\t0u /* unused */");
	for (unsigned logn = 1; logn <= 9; logn++) {
		FT bound = powl(1.1 * 2.0 * ssig, 2) * powl(2, logn + 1);
		printf(", %lldu", (long long) floorl(bound));
	}
	printf("\n};\n");

	return 0;
}
