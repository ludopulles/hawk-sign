#include<bits/stdc++.h>

extern "C" {
	#ifndef restrict
		#define restrict
	#endif
	#include "inner.h"
}

using namespace std;
typedef long long ll;

void randombytes(unsigned char *x, unsigned long long xlen) {
	for (; xlen -- > 0; ++x) *x = ((unsigned char) rand());
}

ll time_diff(const struct timeval *begin, const struct timeval *end) {
	return 1000000LL * (end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec);
}

// Protocol parameters:
constexpr size_t logn = 9;
const fpr isigma_sig = fpr_div(fpr_of(1000), fpr_of(1292));
constexpr int num_samples = 10 * 1000 * 1000;

int main() {
	unsigned char seed[48];
	inner_shake256_context sc;
	struct timeval t0, t1;

	// Initialize a RNG.
	randombytes(seed, sizeof seed);
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);
	
	sampler_context spc;
	spc.sigma_min = fpr_sigma_min[logn];
	Zf(prng_init)(&spc.p, &sc); // Use a fast PRNG for gaussian sampling.

	fpr middle;
	ll us;

	middle = fpr_zero;
	// heat up computer:
	for (int rep = 0; rep < num_samples/5; rep++)
		Zf(sampler)((void *)&spc, middle, isigma_sig);
	middle = fpr_half(fpr_one);
	for (int rep = 0; rep < num_samples/5; rep++)
		Zf(sampler)((void *)&spc, middle, isigma_sig);

	// Do the real test

	middle = fpr_half(fpr_one);
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++)
		Zf(sampler)((void *)&spc, middle, isigma_sig);
	gettimeofday(&t1, NULL);
	us = time_diff(&t0, &t1);
	printf("Time sampling around center = 1/2: %lld us\n", us);

	middle = fpr_zero;
	gettimeofday(&t0, NULL);
	for (int rep = 0; rep < num_samples; rep++)
		Zf(sampler)((void *)&spc, middle, isigma_sig);
	gettimeofday(&t1, NULL);
	us = time_diff(&t0, &t1);
	printf("Time sampling around center = 0:   %lld us\n", us);

	// Conclusion: both sample with the same time consumption
	return 0;
}
