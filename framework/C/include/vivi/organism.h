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

/* Opaque: the layout lives in the implementation and is not part of the
 * stable API. Read it through the accessors below, or via organism_read /
 * organism_peek. */
typedef struct vivi_organism vivi_organism;

/* stable accessors: read an organism without depending on its layout */
size_t organism_count(const vivi_organism *o);
int organism_chr_id_at(const vivi_organism *o, size_t i);
const chr_opts *organism_chr_opts_at(const vivi_organism *o, size_t i);
int organism_generation(const vivi_organism *o);
int organism_max_generation(const vivi_organism *o);
bool organism_is_dead(const vivi_organism *o);
bool organism_is_stem(const vivi_organism *o);
bool organism_has_stem(const vivi_organism *o);

/* advanced research access (homolog 0 or 1). organism_strand returns a
 * borrowed pointer; organism_replace_strand takes ownership of `data` and
 * frees the previous strand (NULL, 0 clears). organism_set_strand_len only
 * rewrites the recorded length (no reallocation). */
const uint8_t *organism_strand(const vivi_organism *o, size_t homolog, size_t i, size_t *len);
uint8_t *organism_strand_mut(vivi_organism *o, size_t homolog, size_t i, size_t *len);
[[nodiscard]] bool organism_replace_strand(vivi_organism *o, size_t homolog, size_t i,
	uint8_t *data, size_t len);
[[nodiscard]] bool organism_set_strand_len(vivi_organism *o, size_t homolog, size_t i, size_t len);
void organism_set_max_generation(vivi_organism *o, int max_gen);

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
/* inspect without repairing, mutating or renewing; the organism keeps its
 * exact strands and size, and out[i] holds the data observed for chr i */
[[nodiscard]] bool organism_peek(vivi_bytes **out, const vivi_organism *o, cell_report *rep, const char **err);
[[nodiscard]] bool organism_checkpoint_checked(vivi_organism *o, cell_report *rep, const char **err);
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
/* VIV14NB4NSH33 writer ("Banshee"): parity plus the per-chromosome
 * on-strand primer barcode stored in the header */
[[nodiscard]] bool organism_serialize14nb(vivi_bytes *out, vivi_organism *o, const char **err);
/* container revision: 1 (VIV1), 14 (VIV14), 141 (VIV14N), 143 (Banshee), 0 */
int organism_container_version(const uint8_t *s, size_t len);
/* reads VIV1, VIV14, VIV14N and VIV14NB4NSH33 containers */
[[nodiscard]] bool organism_deserialize(vivi_organism **out, const uint8_t *s, size_t len, const char **err);
cell_report organism_maintain(vivi_organism **organisms, size_t count,
	const vivi_organism *stem);

#endif /* ORGANISM_H */


