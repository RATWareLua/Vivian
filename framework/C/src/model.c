/* model.c -- population and selection layer over organisms */
#include "vivi/model.h"
#include <string.h>

struct vivi_population {
	int *ids;
	chr_opts *opts;
	size_t nchr, size, max_gen;
	vivi_organism **members;
	double *fit;
};

vivi_population *population_new(const int *ids, const chr_opts *opts,
	const uint8_t *const *seed_datas, const size_t *seed_lens, size_t nchr,
	size_t pop_size, int max_gen, uint32_t diversity, uint32_t seed, const char **err)
{
	static VIVI_THREAD_LOCAL const char *sink;
	if (!err) err = &sink;
	if (!ids || !opts || !seed_datas || !seed_lens || nchr == 0 || pop_size == 0) {
		*err = "invalid population spec";
		return nullptr;
	}
	vivi_population *pop = vivi_zalloc(1, sizeof(*pop));
	if (!pop) { *err = "out of memory"; return nullptr; }
	pop->nchr = nchr;
	pop->size = pop_size;
	pop->max_gen = (max_gen > 0) ? (size_t)max_gen : 60;
	pop->ids = vivi_alloc(nchr * sizeof(int));
	pop->opts = vivi_alloc(nchr * sizeof(chr_opts));
	pop->members = vivi_zalloc(pop_size, sizeof(vivi_organism *));
	pop->fit = vivi_alloc(pop_size * sizeof(double));
	if (!pop->ids || !pop->opts || !pop->members || !pop->fit) {
		population_free(pop);
		*err = "out of memory";
		return nullptr;
	}
	for (size_t i = 0; i < nchr; i++) {
		pop->ids[i] = ids[i];
		pop->opts[i] = opts[i];
	}
	for (size_t m = 0; m < pop_size; m++) {
		vivi_organism *o = nullptr;
		if (!organism_new(&o, ids, opts, seed_datas, seed_lens, nchr, (int)pop->max_gen, err)) {
			population_free(pop);
			return nullptr;
		}
		if (m && diversity) {
			org_mut_result mut = { 0 };
			if (organism_mutate(&mut, o, (int)diversity, seed + (uint32_t)m, err))
				vivi_bytes_free(&mut.data);
		}
		pop->members[m] = o;
		pop->fit[m] = 0.0;
	}
	return pop;
}

void population_free(vivi_population *pop)
{
	if (!pop) return;
	if (pop->members) {
		for (size_t i = 0; i < pop->size; i++) organism_free(pop->members[i]);
		vivi_dealloc(pop->members);
	}
	vivi_dealloc(pop->fit);
	vivi_dealloc(pop->ids);
	vivi_dealloc(pop->opts);
	vivi_dealloc(pop);
}

size_t population_size(const vivi_population *pop) { return pop ? pop->size : 0; }

const vivi_organism *population_at(const vivi_population *pop, size_t i)
{
	return (pop && i < pop->size) ? pop->members[i] : nullptr;
}

double population_fitness_at(const vivi_population *pop, size_t i)
{
	return (pop && i < pop->size) ? pop->fit[i] : -1e308;
}

double population_best(const vivi_population *pop, size_t *index)
{
	if (!pop || pop->size == 0) {
		if (index) *index = 0;
		return -1e308;
	}
	size_t best = 0;
	for (size_t i = 1; i < pop->size; i++)
		if (pop->fit[i] > pop->fit[best]) best = i;
	if (index) *index = best;
	return pop->fit[best];
}

bool population_evaluate(vivi_population *pop, vivi_fitness_fn fn, void *ctx,
	const char **err)
{
	static VIVI_THREAD_LOCAL const char *sink;
	if (!err) err = &sink;
	if (!pop || !fn) { *err = "expected population and fitness"; return false; }
	for (size_t m = 0; m < pop->size; m++) {
		vivi_bytes *data = nullptr;
		cell_report rep;
		if (!organism_peek(&data, pop->members[m], &rep, err)) {
			if (vivi_error_code(*err) == VIVI_ERR_OOM) return false;
			pop->fit[m] = -1e308;
			continue;
		}
		pop->fit[m] = fn(data, pop->nchr, ctx);
		vivi_bytes_free_n(data, pop->nchr);
	}
	return true;
}

bool population_step(vivi_population *pop, vivi_fitness_fn fn, void *ctx,
	uint32_t seed, const char **err)
{
	static VIVI_THREAD_LOCAL const char *sink;
	if (!err) err = &sink;
	if (!pop || !fn) { *err = "expected population and fitness"; return false; }
	if (!population_evaluate(pop, fn, ctx, err)) return false;
	size_t n = pop->size;
	size_t *idx = vivi_alloc(n * sizeof(size_t));
	vivi_organism **nm = vivi_zalloc(n, sizeof(vivi_organism *));
	double *nf = vivi_alloc(n * sizeof(double));
	uint8_t *kept = vivi_zalloc(n, 1);
	if (!idx || !nm || !nf || !kept) {
		vivi_dealloc(idx); vivi_dealloc(nm); vivi_dealloc(nf); vivi_dealloc(kept);
		*err = "out of memory";
		return false;
	}
	for (size_t i = 0; i < n; i++) idx[i] = i;
	for (size_t i = 1; i < n; i++) {
		size_t key = idx[i];
		size_t j = i;
		while (j > 0 && pop->fit[idx[j - 1]] < pop->fit[key]) {
			idx[j] = idx[j - 1];
			j--;
		}
		idx[j] = key;
	}
	size_t keep = n / 2;
	if (keep == 0) keep = 1;
	for (size_t i = 0; i < keep; i++) {
		kept[idx[i]] = 1;
		nm[i] = pop->members[idx[i]];
		nf[i] = pop->fit[idx[i]];
	}
	vivi_prng rng;
	vivi_prng_init(&rng, seed ? seed : 1);
	for (size_t i = keep; i < n; i++) {
		vivi_organism *pa = nm[vivi_prng_next(&rng) % keep];
		vivi_organism *pb = nm[vivi_prng_next(&rng) % keep];
		vivi_organism *child = nullptr;
		if (!organism_cross(&child, pa, pb, vivi_prng_next(&rng), err)
			&& !organism_cross(&child, pa, pa, vivi_prng_next(&rng), err)) {
			for (size_t j = keep; j < i; j++) organism_free(nm[j]);
			vivi_dealloc(idx); vivi_dealloc(nm); vivi_dealloc(nf); vivi_dealloc(kept);
			return false;
		}
		org_mut_result mut = { 0 };
		if (organism_mutate(&mut, child, 1, vivi_prng_next(&rng), err))
			vivi_bytes_free(&mut.data);
		nm[i] = child;
		nf[i] = 0.0;
	}
	for (size_t i = 0; i < n; i++)
		if (!kept[i]) organism_free(pop->members[i]);
	vivi_dealloc(pop->members);
	vivi_dealloc(pop->fit);
	pop->members = nm;
	pop->fit = nf;
	vivi_dealloc(idx);
	vivi_dealloc(kept);
	return population_evaluate(pop, fn, ctx, err);
}

bool population_evolve(vivi_population *pop, vivi_fitness_fn fn, void *ctx,
	int generations, uint32_t seed, const char **err)
{
	static VIVI_THREAD_LOCAL const char *sink;
	if (!err) err = &sink;
	if (!pop || !fn) { *err = "expected population and fitness"; return false; }
	if (generations <= 0) return population_evaluate(pop, fn, ctx, err);
	for (int g = 0; g < generations; g++)
		if (!population_step(pop, fn, ctx, seed + (uint32_t)g + 1u, err)) return false;
	return true;
}
