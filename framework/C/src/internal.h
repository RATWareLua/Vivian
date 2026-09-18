/* internal.h -- private layouts of the model objects.
 *
 * NOT a public header: the stable API in include/vivi/ never exposes these
 * fields; external code reads and writes them through the accessors in
 * cell.h / organism.h. */
#ifndef VIVI_INTERNAL_H
#define VIVI_INTERNAL_H

#include "vivi/vivi.h"
#include "vivi/chromosome.h"

struct vivi_cell {
	int chr_id, generation, max_gen;
	int stem, dead;
	const struct vivi_cell *stem_source;
	uint8_t *hom[2];
	size_t hlen[2];
};

struct vivi_organism {
	int *chr_ids;
	chr_opts *chr_opts;
	size_t nchr;
	uint8_t **hom[2];
	size_t *hlen[2];
	int generation, max_gen;
	int stem, dead;
	const struct vivi_organism *stem_source;
};

#endif /* VIVI_INTERNAL_H */
