/* pool.c -- oligo pool and primer-based random access */
#include "vivi/pool.h"
#include "vivi/genome.h"
#include <string.h>

/* absorbs messages when a caller passes err == nullptr */
static const char *pool_err_sink;

void vivi_pool_free(vivi_pool *p)
{
	if (!p) return;
	for (size_t i = 0; i < p->count; i++) vivi_dealloc(p->strands[i]);
	vivi_dealloc(p->ids);
	vivi_dealloc(p->strands);
	vivi_dealloc(p->lens);
	p->ids = nullptr;
	p->strands = nullptr;
	p->lens = nullptr;
	p->count = 0;
	p->cap = 0;
}

static bool pool_reserve(vivi_pool *p, size_t need)
{
	if (p->cap >= need) return true;
	size_t cap = p->cap ? p->cap * 2 : 8;
	if (cap < need) cap = need;
	int *ids = vivi_alloc(cap * sizeof(int));
	uint8_t **strands = vivi_alloc(cap * sizeof(uint8_t *));
	size_t *lens = vivi_alloc(cap * sizeof(size_t));
	if (!ids || !strands || !lens) {
		vivi_dealloc(ids);
		vivi_dealloc(strands);
		vivi_dealloc(lens);
		return false;
	}
	if (p->count) {
		memcpy(ids, p->ids, p->count * sizeof(int));
		memcpy(strands, p->strands, p->count * sizeof(uint8_t *));
		memcpy(lens, p->lens, p->count * sizeof(size_t));
	}
	vivi_dealloc(p->ids);
	vivi_dealloc(p->strands);
	vivi_dealloc(p->lens);
	p->ids = ids;
	p->strands = strands;
	p->lens = lens;
	p->cap = cap;
	return true;
}

bool vivi_pool_add_gene(vivi_pool *p, int id, const uint8_t *strand, size_t slen,
	const char **err)
{
	if (!err) err = &pool_err_sink;
	if (!p) { *err = "expected pool"; return false; }
	if (!strand) { *err = "expected string"; return false; }
	if (!pool_reserve(p, p->count + 1)) { *err = "out of memory"; return false; }
	uint8_t *copy = vivi_alloc(slen ? slen : 1);
	if (!copy) { *err = "out of memory"; return false; }
	memcpy(copy, strand, slen);
	p->ids[p->count] = id;
	p->strands[p->count] = copy;
	p->lens[p->count] = slen;
	p->count++;
	return true;
}

bool vivi_pool_add_chromosome(vivi_pool *p, const uint8_t *strand, size_t slen,
	const char **err)
{
	if (!err) err = &pool_err_sink;
	genome_scan_result res;
	if (!genome_gene_scan(&res, strand, slen, err)) return false;
	bool ok = true;
	for (size_t i = 0; i < res.count && ok; i++) {
		const genome_gene *g = &res.genes[i];
		if (!g->crc_ok) continue;   /* only intact molecules enter the pool */
		ok = vivi_pool_add_gene(p, g->id, strand + g->offset, g->size, err);
	}
	genome_scan_free(&res);
	return ok;
}

bool vivi_pool_amplify(vivi_amp_result *out, const vivi_pool *pool, int id,
	const vivi_amp_opts *opts, const char **err)
{
	if (!err) err = &pool_err_sink;
	memset(out, 0, sizeof(*out));
	out->id = -1;
	if (!pool || pool->count == 0) { *err = "empty pool"; return false; }
	vivi_amp_opts def = { 0.0, 0.0, 0, { 0.0, 0.0, 0.0, 0.0, 0 } };
	if (!opts) opts = &def;
	if (!(opts->p_access >= 0.0 && opts->p_access <= 1.0)
		|| !(opts->p_cross >= 0.0 && opts->p_cross <= 1.0)) {
		*err = "invalid probability";
		return false;
	}
	size_t target = SIZE_MAX;
	for (size_t i = 0; i < pool->count; i++)
		if (pool->ids[i] == id) { target = i; break; }
	if (target == SIZE_MAX) { *err = "target not in pool"; return false; }
	vivi_prng rng;
	vivi_prng_init(&rng, opts->seed);
	if (vivi_prng_chance(&rng, opts->p_access)) {
		out->read.dropped = 1;
		return true;
	}
	size_t pick = target;
	if (pool->count > 1 && vivi_prng_chance(&rng, opts->p_cross)) {
		pick = (size_t)(vivi_prng_next(&rng) % pool->count);
		if (pool->ids[pick] == id) pick = (pick + 1) % pool->count;
	}
	vivi_channel_opts ch = opts->ch;
	ch.seed = opts->seed + (uint32_t)pick * 0x9E37'79B9u + 1u;
	if (!vivi_channel_read(&out->read, pool->strands[pick], pool->lens[pick], &ch, err))
		return false;
	out->id = pool->ids[pick];
	return true;
}
