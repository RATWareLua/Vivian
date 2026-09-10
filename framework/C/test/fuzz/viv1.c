/* viv1.c -- libFuzzer harness for the VIV1 container parser */
#include <stdlib.h>
#include "vivi/organism.h"

static vivi_bytes g_seed, g_seed14, g_seed14n, g_seed14nb;
static int g_seed_done;

static void seed_cleanup(void)
{
	vivi_bytes_free(&g_seed);
	vivi_bytes_free(&g_seed14);
	vivi_bytes_free(&g_seed14n);
	vivi_bytes_free(&g_seed14nb);
}

static void seed_once(void)
{
	if (g_seed_done) return;
	g_seed_done = 1;
	static const uint8_t d0[] = "fuzz genome alpha";
	static const uint8_t d1[] = "fuzz genome beta";
	int ids[2] = { 0, 1 };
	chr_opts co = { 1024, 0, 3, 4, 0, 0, 0 };
	chr_opts ban = { 1024, 0, 3, 4, 0, 2, 4242 };
	chr_opts opts[2] = { co, ban };
	const uint8_t *datas[2] = { d0, d1 };
	size_t lens[2] = { sizeof(d0) - 1, sizeof(d1) - 1 };
	vivi_organism *o = nullptr;
	if (organism_new(&o, ids, opts, datas, lens, 2, 60, nullptr)) {
		(void)organism_serialize(&g_seed, o, nullptr);
		(void)organism_serialize14(&g_seed14, o, nullptr);
		(void)organism_serialize14n(&g_seed14n, o, nullptr);
		(void)organism_serialize14nb(&g_seed14nb, o, nullptr);
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
	if (g_seed14n.data) load_one(g_seed14n.data, g_seed14n.len);
	if (g_seed14nb.data) load_one(g_seed14nb.data, g_seed14nb.len);
	load_one(data, size);
	return 0;
}
