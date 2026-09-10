/* viv1.c -- libFuzzer harness for the VIV1 container parser */
#include <stdlib.h>
#include "vivi/organism.h"

static vivi_bytes g_seed, g_seed14;
static int g_seed_done;

static void seed_cleanup(void)
{
	vivi_bytes_free(&g_seed);
	vivi_bytes_free(&g_seed14);
}

static void seed_once(void)
{
	if (g_seed_done) return;
	g_seed_done = 1;
	static const uint8_t d0[] = "fuzz genome alpha";
	static const uint8_t d1[] = "fuzz genome beta";
	int ids[2] = { 0, 1 };
	chr_opts co = { 1024, 0, 3, 4, 0 };
	chr_opts opts[2] = { co, co };
	const uint8_t *datas[2] = { d0, d1 };
	size_t lens[2] = { sizeof(d0) - 1, sizeof(d1) - 1 };
	vivi_organism *o = nullptr;
	if (organism_new(&o, ids, opts, datas, lens, 2, 60, nullptr)) {
		(void)organism_serialize(&g_seed, o, nullptr);
		(void)organism_serialize14(&g_seed14, o, nullptr);
		organism_free(o);
	}
	atexit(seed_cleanup);
}

static void load_one(const uint8_t *data, size_t size)
{
	vivi_organism *o = nullptr;
	if (!organism_deserialize(&o, data, size, nullptr)) return;
	vivi_bytes *rd = nullptr;
	cell_report rep;
	if (organism_read(&rd, o, &rep, nullptr)) vivi_bytes_free_n(rd, o->nchr);
	organism_free(o);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
	seed_once();
	if (g_seed.data) load_one(g_seed.data, g_seed.len);
	if (g_seed14.data) load_one(g_seed14.data, g_seed14.len);
	load_one(data, size);
	return 0;
}
