/* parity.c -- systematic Reed-Solomon over GF(256), Lagrange form */
#include "vivi/parity.h"
#include <string.h>

/* absorbs messages when a caller passes err == nullptr */
static const char *par_err_sink;

/* GF(256) with the classic 0x11D polynomial and generator 2 */
#define GF_POLY 0x11Du

static uint8_t gf_exp[512];
static uint8_t gf_log[256];
static int gf_ready;

static void gf_init(void)
{
	if (gf_ready) return;
	uint32_t x = 1;
	for (int i = 0; i < 255; i++) {
		gf_exp[i] = (uint8_t)x;
		gf_log[x] = (uint8_t)i;
		x <<= 1;
		if (x & 0x100u) x ^= GF_POLY;
	}
	for (int i = 255; i < 512; i++) gf_exp[i] = gf_exp[i - 255];
	gf_ready = 1;
}

static uint8_t gf_mul(uint8_t a, uint8_t b)
{
	if (!a || !b) return 0;
	int t = (int)gf_log[a] + (int)gf_log[b];
	if (t >= 255) t -= 255;
	return gf_exp[t];
}

static uint8_t gf_inv(uint8_t a)
{
	return gf_exp[255 - (int)gf_log[a]];
}

/* Lagrange basis coefficients at x over the points xs[0..n-1] */
static void lagrange_coefs(uint8_t *coef, const uint8_t *xs, size_t n, uint8_t x)
{
	for (size_t i = 0; i < n; i++) {
		uint8_t num = 1, den = 1;
		for (size_t k = 0; k < n; k++) {
			if (k == i) continue;
			num = gf_mul(num, (uint8_t)(x ^ xs[k]));
			den = gf_mul(den, (uint8_t)(xs[i] ^ xs[k]));
		}
		coef[i] = gf_mul(num, gf_inv(den));
	}
}

bool vivi_parity_encode(uint8_t **shards, const uint8_t *const *data,
	size_t n, size_t m, size_t len, const char **err)
{
	if (!err) err = &par_err_sink;
	if (!shards || !data) { *err = "expected buffers"; return false; }
	if (n < 1 || len < 1 || n + m > 255) { *err = "invalid shard geometry"; return false; }
	for (size_t k = 0; k < n + m; k++) shards[k] = nullptr;
	gf_init();
	uint8_t xs[255];
	for (size_t i = 0; i < n; i++) xs[i] = (uint8_t)(i + 1);
	for (size_t i = 0; i < n; i++) {
		shards[i] = vivi_alloc(len);
		if (!shards[i]) goto oom;
		memcpy(shards[i], data[i], len);
	}
	for (size_t j = 0; j < m; j++) {
		uint8_t coef[255];
		lagrange_coefs(coef, xs, n, (uint8_t)(n + j + 1));
		shards[n + j] = vivi_zalloc(len, 1);
		if (!shards[n + j]) goto oom;
		for (size_t b = 0; b < len; b++) {
			uint8_t acc = 0;
			for (size_t i = 0; i < n; i++)
				acc ^= gf_mul(coef[i], data[i][b]);
			shards[n + j][b] = acc;
		}
	}
	return true;
oom:
	vivi_parity_release(shards, n + m);
	*err = "out of memory";
	return false;
}

bool vivi_parity_decode(uint8_t **out, const uint8_t *const *shards,
	const uint8_t *present, size_t n, size_t m, size_t len, const char **err)
{
	if (!err) err = &par_err_sink;
	if (!out || !shards || !present) { *err = "expected buffers"; return false; }
	if (n < 1 || len < 1 || n + m > 255) { *err = "invalid shard geometry"; return false; }
	for (size_t i = 0; i < n; i++) out[i] = nullptr;
	gf_init();
	uint8_t px[255];
	size_t np = 0;
	for (size_t k = 0; k < n + m && np < n; k++) {
		if (!present[k]) continue;
		if (!shards[k]) { *err = "present shard is null"; return false; }
		px[np++] = (uint8_t)(k + 1);
	}
	if (np < n) { *err = "too few shards"; return false; }
	for (size_t i = 0; i < n; i++) {
		out[i] = vivi_alloc(len);
		if (!out[i]) goto oom;
		if (present[i]) {
			memcpy(out[i], shards[i], len);
			continue;
		}
		uint8_t coef[255];
		lagrange_coefs(coef, px, n, (uint8_t)(i + 1));
		for (size_t b = 0; b < len; b++) {
			uint8_t acc = 0;
			size_t k = 0;
			for (size_t s = 0; s < n + m && k < n; s++) {
				if (!present[s]) continue;
				acc ^= gf_mul(coef[k], shards[s][b]);
				k++;
			}
			out[i][b] = acc;
		}
	}
	return true;
oom:
	vivi_parity_release(out, n);
	*err = "out of memory";
	return false;
}

void vivi_parity_release(uint8_t **shards, size_t count)
{
	if (!shards) return;
	for (size_t k = 0; k < count; k++) {
		vivi_dealloc(shards[k]);
		shards[k] = nullptr;
	}
}
