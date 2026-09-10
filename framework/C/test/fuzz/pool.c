/* pool.c -- libFuzzer harness for pool building and primer amplification */
#include <stdlib.h>
#include "vivi/organism.h"
#include "vivi/pool.h"

static vivi_bytes g_seed;
static int g_seed_done;

static void seed_cleanup(void) { vivi_bytes_free(&g_seed); }

static void seed_once(void)
{
	if (g_seed_done) return;
	g_seed_done = 1;
	static const uint8_t payload[] = "pool fuzz seed payload, several genes long";
	chr_opts co = { 16, 0, 3, 4, 0, 0, 0 };
	(void)chr_encode(&g_seed, 1, payload, sizeof(payload) - 1, &co, nullptr);
	atexit(seed_cleanup);
}

static void amplify_one(const uint8_t *data, size_t size, const vivi_amp_opts *ao)
{
	vivi_pool pool = { 0 };
	if (!vivi_pool_add_chromosome(&pool, data, size, nullptr)) {
		vivi_pool_free(&pool);
		return;
	}
	if (pool.count) {
		int id = pool.ids[0];
		vivi_amp_result ar;
		if (vivi_pool_amplify(&ar, &pool, id, ao, nullptr))
			vivi_bytes_free(&ar.read.strand);
	}
	vivi_pool_free(&pool);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
	seed_once();
	vivi_amp_opts ao = { 0.05, 0.01, (uint32_t)size + 1u,
		{ 0.01, 0.001, 0.001, 0.01, 2u }, 0.01 };
	if (g_seed.data) amplify_one(g_seed.data, g_seed.len, &ao);
	amplify_one(data, size, &ao);
	return 0;
}
