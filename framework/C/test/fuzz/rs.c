#include <stdlib.h>
#include <string.h>
#include "vivi/rs.h"

static vivi_bytes g_seed;
static int g_seed_done;

static void seed_cleanup(void)
{
	vivi_bytes_free(&g_seed);
}

static void seed_once(void)
{
	if (g_seed_done) return;
	g_seed_done = 1;
	static const uint8_t payload[32] = {
		0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
		0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff,
		0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef,
		0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10
	};
	(void)rs_encode(&g_seed, payload, 32, 8, nullptr);
	atexit(seed_cleanup);
}

static void decode_one(const uint8_t *data, size_t size)
{
	size_t n = size < 255 ? size : 255;
	if (n <= 32) return;
	uint8_t buf[255];
	memcpy(buf, data, n);
	size_t corrected = 0;
	(void)rs_decode(buf, n, 32, &corrected, nullptr);
}

static void encode_one(const uint8_t *data, size_t size)
{
	size_t k = size < 200 ? size : 200;
	size_t m = size ? (size_t)(data[0] % 64) + 1 : 8;
	if (k + m > 255) m = 255 - k;
	if (m < 1 || m > 255) return;
	vivi_bytes cw = { 0 };
	if (rs_encode(&cw, data, k, m, nullptr)) vivi_bytes_free(&cw);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
	seed_once();
	if (g_seed.data) {
		uint8_t buf[255];
		memcpy(buf, g_seed.data, g_seed.len);
		size_t corrected = 0;
		(void)rs_decode(buf, g_seed.len, 32, &corrected, nullptr);
	}
	decode_one(data, size);
	encode_one(data, size);
	return 0;
}
