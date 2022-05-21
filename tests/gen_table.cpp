/*
 * Prints the precomputed table for discrete gaussian sampling in keygen.c
 * For generating the CDT for signing, use the sage script `renyi.sage`, it
 * requires higher precision.
 */
#include<bits/stdc++.h>
typedef long double FT;
const int LEN = 100;
const FT spk = 1.500, ssec = 1.425, scale = powl(2, 63);

FT rho(FT x, FT sigma) { return expl(-(x * x) / (2.0 * sigma * sigma)); }

int main() {
	std::map<int, FT> table;
	unsigned long long results[LEN] = {};

	FT weight = 0.0, sum = 0.0;
	for (int x = LEN; x >= -LEN; x--) {
		table[x] = rho(x, spk);
		weight += table[x];
	}

	for (int k = LEN; k -- > 0; ) {
		// results[k] ~ P( |X| >= k + 1)
		results[k] = floorl(scale * sum / weight);
		sum += table[k] + table[-k];
	}

	int len = 0;
	while (results[len] != 0) len++;
	printf("static const uint64_t gauss_keygen[%d] = {\n\t", len);
	for (int x = 0; x < len; x++) {
		printf("%lluu, ", results[x]);
	}
	printf("\n};\n\n");

	/*
	 * Generate l2bounds using:
	 *     l2bound(logn) = floor( (2 * sigma_sec)^2 * 2n ).
	 */
	printf("const uint32_t Zf(l2bound)[10] = {\n\t0u /* unused */");
	for (unsigned logn = 1; logn <= 9; logn++) {
		FT bound = powl(2.0 * ssec, 2) * powl(2, logn + 1);
		printf(", %lldu", (long long) floorl(bound));
	}
	printf("\n};\n");
	return 0;
}
