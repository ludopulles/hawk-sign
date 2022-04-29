#include<bits/stdc++.h>
int main() {
	long double spk = 1.5, p63 = powl(2, 63), table[100], csum[100];
	unsigned long long results[2][25] = {};
	for (int x = 100; x --> 0; ) {
		table[x] = expl(-0.5* x*x / spk / spk);
		csum[x] = table[x];
		if (x < 99) csum[x] += csum[x+1];
	}
	results[0][0] = llroundl(p63 * table[0] / (csum[0] + csum[1]));
	for (int x = 24; --x >= 1; ) {
		results[0][x] = llroundl(p63 * csum[1+x] / csum[1]);
	}

	int len = 0;
	while (results[0][len] != 0) len++;
	len++;
	printf("static const uint64_t gauss_1500[%d] = {\n", len);
	for (int x = 0; x < len; x++)
		printf("\t%19lluu,\n", results[0][x]);
	printf("};\n\n");

	long double ssig = 1.292, mu = 0;
	for (int i = 0; i < 2; i++, mu += 0.5) {
		for (int x = 100; x --> 0; ) {
			table[x] = expl(-0.5* (x-mu)*(x-mu) / ssig / ssig);
			csum[x] = table[x];
			if (x < 99) csum[x] += csum[x+1];
		}
		results[i][0] = llroundl(p63 * table[i] / (csum[i] + csum[1]));
		for (int x = 19; --x >= 1; )
			results[i][x] = llroundl(p63 * csum[1+i+x] / csum[1+i]);
	}
	printf("static const uint64_t gauss_1292[26] = {\n");
	for (int x = 0; x < 13; x++)
		printf("\t%19lluu, %19lluu,\n", results[0][x], results[1][x]);
	printf("};\n");

}

