/* model.h -- population and selection layer over organisms (research).
 *
 * A genome is raw data; a fitness function scores one genome. The
 * population keeps a fixed number of organisms, evaluates them with
 * organism_peek (read-only), keeps the fittest half and breeds the rest
 * with organism_cross plus a point mutation. Deterministic per seed.
 */
#ifndef MODEL_H
#define MODEL_H

#include "vivi.h"
#include "organism.h"

typedef struct vivi_population vivi_population;

/* fitness of one genome: data[i] holds chromosome i's raw bytes */
typedef double (*vivi_fitness_fn)(const vivi_bytes *data, size_t nchr, void *ctx);

/* every member starts from the same seed genome; members after the first
 * are diversified with `diversity` raw bit-flips */
[[nodiscard]] vivi_population *population_new(const int *ids, const chr_opts *opts,
	const uint8_t *const *seed_datas, const size_t *seed_lens, size_t nchr,
	size_t pop_size, int max_gen, uint32_t diversity, uint32_t seed, const char **err);

void population_free(vivi_population *pop);
size_t population_size(const vivi_population *pop);
const vivi_organism *population_at(const vivi_population *pop, size_t i);
double population_fitness_at(const vivi_population *pop, size_t i);

[[nodiscard]] bool population_evaluate(vivi_population *pop, vivi_fitness_fn fn, void *ctx,
	const char **err);
[[nodiscard]] bool population_step(vivi_population *pop, vivi_fitness_fn fn, void *ctx,
	uint32_t seed, const char **err);
[[nodiscard]] bool population_evolve(vivi_population *pop, vivi_fitness_fn fn, void *ctx,
	int generations, uint32_t seed, const char **err);

/* best stored fitness after the last evaluate/step; index may be NULL */
double population_best(const vivi_population *pop, size_t *index);

#endif /* MODEL_H */
