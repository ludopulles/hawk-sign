#include<bits/stdc++.h>
using namespace std;

typedef long double FT; 
constexpr FT sigma_sig = 1.292;

// Returns the error function for a normal distributed variable with standard deviation sigma,
// i.e. the probability that a sample X has value that does NOT lie in (-x, x).
FT erfcl(FT x, FT sigma) {
	return erfcl(x / (sigma * sqrtl(2)));
}

FT F(FT iq00, long long n) {
	FT val = 1.0 - powl(erfl(0.5 / sqrtl(2*iq00) / sigma_sig), n);
	if (val < 1e-50)
		return n * erfcl(0.5, sqrt(iq00) * sigma_sig);
	return val;
}

signed main() {
	long long n;
	cin >> n;
	FT target = powl(0.5, 96);

	FT iq00_lo = 0.00001, iq00_hi = 1.0;
	while (iq00_hi - iq00_lo > 1e-10) {
		FT mid = (iq00_lo + iq00_hi) * 0.5;
		if (F(mid, n) > target)
			iq00_hi = mid;
		else
			iq00_lo = mid;
	}

	cout << iq00_lo << ": " << F(iq00_lo, n) << endl;
	cout << iq00_hi << ": " << F(iq00_hi, n) << endl;

	FT x;
	while (cin >> x) {
		cout << x << ": " << F(x, n) << endl;
	}
    return 0;
}
