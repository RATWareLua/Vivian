/* organism.h -- the organism layer (C23 port of framework/organism.lua).
 *
 * An organism is a diploid genome of N chromosomes: homologs are
 * tables of chromosome strands. Data-level evolution API: mutate()
 * flips random bits in the RAW data of a random chromosome and
 * re-encodes it (the mutation becomes a legitimate part of the
 * genome); cross(parentA, parentB) builds a per-gene mosaic child.
 * Serialization: the versioned "VIV1" container.
 */
#ifndef ORGANISM_H
#define ORGANISM_H

#include "vivi.h"
#include "cell.h"

typedef struct vivi_organism {
	int *chr_ids;
	chr_opts *chr_opts;
	size_t nchr;
	uint8_t **hom[2];   /* per homolog: array of nchr strand pointers */
	size_t *hlen[2];    /* per homolog: array of strand lengths */
	int generation, max_gen;
	int stem, dead;
	const struct vivi_organism *stem_source;
} vivi_organism;

typedef struct {
	int chr;            /* chromosome id */
	vivi_bytes data;     /* mutated raw data (owned) */
} org_mut_result;

[[nodiscard]] bool organism_new(vivi_organism **out, const int *ids, const chr_opts *opts,
	const uint8_t *const *datas, const size_t *lens, size_t nchr,
	int max_gen, const char **err);
[[nodiscard]] bool organism_stem(vivi_organism **out, const int *ids, const chr_opts *opts,
	const uint8_t *const *datas, const size_t *lens, size_t nchr,
	int max_gen, const char **err);
void organism_attach_stem(vivi_organism *o, const vivi_organism *stem);
bool organism_renew(vivi_organism *o, const vivi_organism *stem, const char **err);
void organism_kill(vivi_organism *o);
void organism_free(vivi_organism *o);

/* reads every chromosome through the repair machinery; on success
 * out[i] holds the raw data of chr_ids[i] */
[[nodiscard]] bool organism_read(vivi_bytes **out, vivi_organism *o, cell_report *rep, const char **err);
cell_report organism_checkpoint(vivi_organism *o);
[[nodiscard]] bool organism_replicate(vivi_organism *o, const char **err);
[[nodiscard]] bool organism_mitosis(vivi_organism **out, vivi_organism *o, const char **err);
[[nodiscard]] bool organism_damage(vivi_organism *o, int count, uint32_t seed, const char **err);
[[nodiscard]] bool organism_mutate(org_mut_result *out, vivi_organism *o, int count, uint32_t seed,
	const char **err);
[[nodiscard]] bool organism_cross(vivi_organism **out, vivi_organism *pa, vivi_organism *pb,
	uint32_t seed, const char **err);
[[nodiscard]] bool organism_serialize(vivi_bytes *out, vivi_organism *o, const char **err);
/* VIV14 writer: stores both homologs and max_gen. organism_serialize()
 * keeps writing VIV1, byte-compatible with the Lua reference. */
[[nodiscard]] bool organism_serialize14(vivi_bytes *out, vivi_organism *o, const char **err);
/* VIV14N writer: as VIV14 plus the per-chromosome parity count */
[[nodiscard]] bool organism_serialize14n(vivi_bytes *out, vivi_organism *o, const char **err);
/* container revision: 1 (VIV1), 14 (VIV14), 141 (VIV14N), 143 (reserved), 0 */
int organism_container_version(const uint8_t *s, size_t len);
/* reads VIV1, VIV14 and VIV14N containers */
[[nodiscard]] bool organism_deserialize(vivi_organism **out, const uint8_t *s, size_t len, const char **err);
cell_report organism_maintain(vivi_organism **organisms, size_t count,
	const vivi_organism *stem);

#endif /* ORGANISM_H */


