/* chromosome.h -- the chromosome layer (C23 port of framework/chromosome.lua).
 *
 * [TELOMERE_L][CENTROMERE][gene 0]..[gene n-1][TELOMERE_R]
 * Telomeres: TTAGGG-style repeats (digits 3,3,0,2,2,2). Centromere:
 * marker + codon-coded header (id, flags, gene count, generation)
 * protected by its own Chaskey tag. Genes are found by scanning.
 */
#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include "vivi.h"
#include "genome.h"

typedef struct {
	int gene_raw;   /* 0 = default 1024; else >= 16 */
	int codon;      /* payload mode: 0 = dense, 1 = codon */
	int h;          /* dense gene homopolymer limit, 0 = default 3 */
	int units;      /* telomere repeat units, 0 = default 4 */
	int flags;      /* user byte */
} chr_opts;

typedef struct {
	int id, flags, ngenes, generation;
	int telomere_ok, cen_ok;
	size_t telomere_bytes;
	genome_gene *genes;
	size_t gene_count;
} chr_record;

[[nodiscard]] bool chr_encode(vivi_bytes *out, int chr_id, const uint8_t *data, size_t len,
	const chr_opts *opts, const char **err);
/* units_hint < 0 = autodetect; record genes point into malloc'd copies */
[[nodiscard]] bool chr_parse(chr_record *out, const uint8_t *strand, size_t slen, int units_hint,
	const char **err);
void chr_record_free(chr_record *r);
[[nodiscard]] bool chr_read(vivi_bytes *out, const chr_record *rec, const char **err);
[[nodiscard]] bool chr_set_generation(vivi_bytes *out, const uint8_t *strand, size_t slen,
	int generation, const char **err);

#endif /* CHROMOSOME_H */


