/* cell.h -- the diploid cell and repair (C23 port of framework/cell.lua).
 *
 * Every chromosome is kept as a HOMOLOG PAIR. checkpoint() repairs a
 * damaged gene by splicing the healthy homolog at the same locus
 * (homologous recombination; damage never shifts the layout). Structural
 * damage (telomeres/centromere) is spliced from the healthy homolog at
 * fixed positions. A locus damaged in BOTH homologs is dead -- the cell
 * is regenerated from a stem cell (germ line).
 */
#ifndef CELL_H
#define CELL_H

#include "vivi.h"
#include "chromosome.h"

typedef struct {
	int repaired, dead, structural, anomaly;
	int renewed;   /* set by read/maintain when a stem cell was used */
	int failed;
} cell_report;

/* Opaque: the layout lives in the implementation and is not part of the
 * stable API. Read it through the accessors below, or via cell_read. */
typedef struct vivi_cell vivi_cell;

/* stable accessors: read a cell without depending on its layout */
int cell_chr_id(const vivi_cell *c);
int cell_generation(const vivi_cell *c);
int cell_max_generation(const vivi_cell *c);
void cell_set_max_generation(vivi_cell *c, int max_gen);
bool cell_is_dead(const vivi_cell *c);
bool cell_is_stem(const vivi_cell *c);

/* advanced strand access (research layer): homolog 0 or 1. cell_strand
 * returns a borrowed pointer; cell_replace_strand takes ownership of
 * `data` and frees the previous strand (pass NULL, 0 to clear). */
const uint8_t *cell_strand(const vivi_cell *c, size_t homolog, size_t *len);
uint8_t *cell_strand_mut(vivi_cell *c, size_t homolog, size_t *len);
[[nodiscard]] bool cell_replace_strand(vivi_cell *c, size_t homolog, uint8_t *data, size_t len);

/* low-level machinery (also used by organism.c) */
[[nodiscard]] bool cell_damage_strand(vivi_bytes *out, const uint8_t *strand, size_t slen, int count,
	uint32_t seed, const char **err);
[[nodiscard]] bool cell_repair_homolog(vivi_bytes *out, const uint8_t *dst, size_t dlen,
	const uint8_t *src, size_t slen, cell_report *rep, const char **err);
[[nodiscard]] bool cell_checkpoint_pair(vivi_bytes *a_out, vivi_bytes *b_out,
	const uint8_t *h1, size_t l1, const uint8_t *h2, size_t l2,
	cell_report *rep, const char **err);

/* cell lifecycle; all produced strands/chromosomes are malloc'd */
[[nodiscard]] bool cell_new(vivi_cell **out, int chr_id, const uint8_t *data, size_t len,
	const chr_opts *opts, const char **err);
[[nodiscard]] bool cell_stem(vivi_cell **out, int chr_id, const uint8_t *data, size_t len,
	const chr_opts *opts, const char **err);
void cell_attach_stem(vivi_cell *c, const vivi_cell *stem);
[[nodiscard]] bool cell_renew(vivi_cell *c, const vivi_cell *stem, const char **err);
void cell_kill(vivi_cell *c);
void cell_free(vivi_cell *c);

[[nodiscard]] bool cell_read(vivi_bytes *out, vivi_cell *c, cell_report *rep, const char **err);
/* inspect without repairing, mutating or renewing: reads whichever homolog
 * is intact and reports what is observed */
[[nodiscard]] bool cell_peek(vivi_bytes *out, const vivi_cell *c, cell_report *rep, const char **err);
[[nodiscard]] bool cell_checkpoint_checked(vivi_cell *c, cell_report *rep, const char **err);
cell_report cell_checkpoint(vivi_cell *c);
[[nodiscard]] bool cell_replicate(vivi_cell *c, const char **err);
[[nodiscard]] bool cell_mitosis(vivi_cell **out, vivi_cell *c, const char **err);
[[nodiscard]] bool cell_damage(vivi_cell *c, int count, uint32_t seed, const char **err);
cell_report cell_maintain(vivi_cell **cells, size_t count, const vivi_cell *stem);

#endif /* CELL_H */


