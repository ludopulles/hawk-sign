// TODO: remove since it is not used.

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_big(uint32_t *x, size_t xlen)
{
	uint32_t *y = malloc(xlen * sizeof(uint32_t));
	memcpy(y, x, xlen * sizeof(uint32_t));

	uint8_t *out = malloc(xlen * 31);
	size_t nout = 0, ylen = xlen;

	size_t sgn = y[ylen-1]>>30;
	if (sgn) {
		for (size_t u = 0; u < ylen; u++)
			y[u] = (1U<<31) - 1 - y[u];
		y[0]++;
		printf("-");
	}

	while (ylen > 0) {
		uint64_t C = 0;
		for (size_t p = ylen; p -- > 0; ) {
			C = (C << 31) | y[p];
			y[p] = (C / 10);
			// new carry:
			C %= 10;
		}
		// now C is the remainder
		out[nout++] = (int8_t) C;
		while (ylen > 0 && y[ylen - 1] == 0) --ylen;
	}
	while (nout -- > 0)
		printf("%d", out[nout]);
	free(y);
	free(out);
}
