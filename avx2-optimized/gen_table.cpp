#include<bits/stdc++.h>
int main() {
	long double sigma = 1.292, mu = 0, p63 = powl(2, 63), table[100], csum[100];
	unsigned long long results[2][20] = {};
	for (int i = 0; i < 2; i++, mu += 0.5) {
		for (int x = 100; x --> 0; ) {
			table[x] = expl(-0.5* (x-mu)*(x-mu) / sigma / sigma);
			csum[x] = table[x];
			if (x < 99) csum[x] += csum[x+1];
		}
		results[i][0] = llroundl(p63 * table[i] / (csum[i] + csum[1]));
		for (int x = 19; --x >= 1; )
			results[i][x] = llroundl(p63 * csum[1+i+x] / csum[1+i]);
	}
	for (int x = 0; x < 20; x++)
		printf("\t%19lluu, %19lluu\n", results[0][x], results[1][x]);
}
