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
} cell_report;

typedef struct vivi_cell {
	int chr_id, generation, max_gen;
	int stem, dead;
	const struct vivi_cell *stem_source;
	uint8_t *hom[2];
	size_t hlen[2];
} vivi_cell;

/* low-level machinery (also used by organism.c) */
[[nodiscard]] [[nodiscard]] bool cell_damage_strand(vivi_bytes *out, const uint8_t *strand, size_t slen, int count,
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
cell_report cell_checkpoint(vivi_cell *c);
[[nodiscard]] bool cell_replicate(vivi_cell *c, const char **err);
[[nodiscard]] bool cell_mitosis(vivi_cell **out, vivi_cell *c, const char **err);
[[nodiscard]] bool cell_damage(vivi_cell *c, int count, uint32_t seed, const char **err);
cell_report cell_maintain(vivi_cell **cells, size_t count, const vivi_cell *stem);

#endif /* CELL_H */


