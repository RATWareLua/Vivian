/* chromosome.h -- the chromosome layer (C23 port of framework/chromosome.lua).
 *
 * [TELOMERE_L][CENTROMERE][gene 0]..[gene n-1][TELOMERE_R]
 * Telomeres: TTAGGG-style repeats (digits 3,3,0,2,2,2). Centromere:
 * marker + codon-coded header (id, flags, gene count, generation)
 * protected by its own Chaskey tag. Genes are found by scanning.
 * VIV14NB4NSH33 ("Banshee") adds a primer site on each side:
 * [TELO][PRIMER_L][CENTROMERE][genes][PRIMER_R][TELO], so a molecule
 * carries its own random-access address (chr_opts.primer).
 */
#ifndef CHROMOSOME_H
#define CHROMOSOME_H

#include "vivi.h"
#include "genome.h"

#define CENBYTES 18   /* VIV1/VIV14 centromere: marker (3) + 20 codons (15) */
#define CEN2BYTES 27  /* VIV14N centromere: marker (3) + 32 codons (24) */
#define PRIMERBYTES 15 /* VIV14NB4NSH33 primer site: marker (3) + 16 codons (12) */

typedef struct {
	int gene_raw;   /* 0 = default 1024; else >= 16 */
	int codon;      /* payload mode: 0 = dense, 1 = codon */
	int h;          /* dense gene homopolymer limit, 0 = default 3 */
	int units;      /* telomere repeat units, 0 = default 4 */
	int flags;      /* user byte */
	int parity;     /* VIV14N: parity genes appended (0..16, 0 = off) */
	int primer;     /* VIV14NB4NSH33: on-strand barcode 1..65535 (0 = off) */
} chr_opts;

typedef struct {
	int id, flags, ngenes, generation;
	int parity;          /* VIV14N: number of parity genes */
	size_t rawlen;       /* VIV14N: true payload length */
	int cen_version;     /* 1 = VIV1/VIV14 centromere, 2 = VIV14N */
	int primer;          /* VIV14NB4NSH33: barcode, 0 = no primer sites */
	int primer_ok;       /* both sites valid and carrying the same barcode */
	int primer_bytes;    /* PRIMERBYTES when the left site anchors the layout */
	int telomere_ok, cen_ok;
	size_t telomere_bytes;
	genome_gene *genes;
	size_t gene_count;
} chr_record;

/* physical random access: a molecule is amplifiable when its primer pair
 * survived (strands without primer sites are pool-addressable as before) */
static inline bool chr_amplifiable(const chr_record *rec)
{
	return rec->primer == 0 || rec->primer_ok;
}

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


