/* parity.c -- libFuzzer harness for the erasure code (geometry + pattern) */
#include <stdlib.h>
#include "vivi/parity.h"

int LLVMFuzzerTestOneInput(const uint8_t *fdata, size_t size)
{
	if (size < 8) return 0;
	size_t n = 1 + (size_t)(fdata[0] & 15);   /* 1..16 data shards */
	size_t m = (size_t)(fdata[1] & 7);        /* 0..7 parity shards */
	size_t len = (size - 2) / (n + m);
	if (len == 0) return 0;

	const uint8_t *data[16];
	for (size_t i = 0; i < n; i++) data[i] = fdata + 2 + i * len;

	uint8_t *shards[23];
	if (!vivi_parity_encode(shards, (const uint8_t *const *)data, n, m, len, nullptr))
		return 0;

	uint8_t present[23];
	size_t pos = 2 + n * len;
	for (size_t k = 0; k < n + m; k++)
		present[k] = (pos < size) ? (uint8_t)(fdata[pos++] & 1u) : 1u;

	uint8_t *out[16];
	if (vivi_parity_decode(out, (const uint8_t *const *)shards, present, n, m, len, nullptr))
		vivi_parity_release(out, n);
	vivi_parity_release(shards, n + m);
	return 0;
}
