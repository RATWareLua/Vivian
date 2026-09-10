/* chr.c -- libFuzzer harness for chromosome parsing and reading */
#include <stdlib.h>
#include "vivi/organism.h"

static vivi_bytes g_seed, g_seed_par;
static int g_seed_done;

static void seed_cleanup(void)
{
	vivi_bytes_free(&g_seed);
	vivi_bytes_free(&g_seed_par);
}

static void seed_once(void)
{
	if (g_seed_done) return;
	g_seed_done = 1;
	static const uint8_t payload[] = "chromosome fuzz seed payload";
	chr_opts co = { 64, 0, 3, 4, 0, 0 };
	(void)chr_encode(&g_seed, 3, payload, sizeof(payload) - 1, &co, nullptr);
	chr_opts par = { 64, 0, 3, 4, 0, 2 };
	(void)chr_encode(&g_seed_par, 4, payload, sizeof(payload) - 1, &par, nullptr);
	atexit(seed_cleanup);
}

static void parse_one(const uint8_t *data, size_t size)
{
	chr_record rec;
	if (!chr_parse(&rec, data, size, -1, nullptr)) return;
	vivi_bytes out = { 0 };
	if (chr_read(&out, &rec, nullptr)) vivi_bytes_free(&out);
	vivi_bytes gen = { 0 };
	if (chr_set_generation(&gen, data, size, 7, nullptr)) vivi_bytes_free(&gen);
	chr_record_free(&rec);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
	seed_once();
	if (g_seed.data) parse_one(g_seed.data, g_seed.len);
	if (g_seed_par.data) parse_one(g_seed_par.data, g_seed_par.len);
	parse_one(data, size);
	return 0;
}
