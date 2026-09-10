/* genome.h -- the gene layer (C23 port of framework/genome.lua).
 *
 * A gene is a self-delimiting data segment framed by promoter/
 * terminator markers, with a codon-coded header and a Chaskey-12
 * integrity tag over the raw data:
 *   [PROM 12 bases][header 12 codons = 6 bytes][payload][tag 8 codons
 *   = 4 bytes][TERM 12 bases]
 *
 * Codon code: codon = 3 bases (b1, b2, wobble); value = b1*4+b2 (16
 * classes, 4 bits). The wobble base carries no data -- any corruption
 * of it leaves the value unchanged (wobble != b2, capping homopolymer
 * runs at 3 inside codon sections).
 */
#ifndef GENOME_H
#define GENOME_H

#include "vivi.h"
#include "dna.h"

/* Chaskey-12 (ISO/IEC 29192-6), official MAC mode, tag truncated to 32
 * bits. key may be NULL for the built-in fixed (public) default key. */
void genome_chaskey_tag(uint32_t out[4], const uint8_t *s, size_t len,
	const uint32_t key[4]);
uint32_t genome_chaskey32(const uint8_t *s, size_t len);

/* 16-bit values <-> codons (count must be a multiple of 4) */
[[nodiscard]] bool genome_pack_codons(vivi_bytes *out, const int *values, size_t nvalues, const char **err);
[[nodiscard]] bool genome_unpack_codons(int *values, size_t nvalues, const uint8_t *packed, size_t plen,
	const char **err);

[[nodiscard]] bool genome_values_from_bytes(int *values, size_t nvalues, const uint8_t *s, size_t len);
[[nodiscard]] bool genome_bytes_from_values(vivi_bytes *out, const int *values, size_t nvalues, size_t nbytes,
	const char **err);

/* mode: codon_mode != 0 -> codon payload; h is the dense payload's
 * homopolymer limit (3..12) */
[[nodiscard]] bool genome_gene_encode(vivi_bytes *out, int id, int usertype, int codon_mode, int h,
	const uint8_t *data, size_t len, const char **err);

typedef struct {
	int id, type;
	int codon;        /* payload mode: 0 = dense, 1 = codon */
	int rawlen, packedlen;
	size_t offset;    /* 0-based promoter byte index in the strand */
	size_t size;      /* total gene bytes: 21 + packedlen */
	int crc_ok;
	uint32_t tag;     /* expected tag (decoded from the codon section) */
	uint8_t *data;    /* decoded raw data (NULL if undecodable) */
	size_t data_len;
} genome_gene;

typedef struct {
	genome_gene *genes;
	size_t count;
} genome_scan_result;

void genome_scan_free(genome_scan_result *r);
[[nodiscard]] bool genome_gene_scan(genome_scan_result *out, const uint8_t *strand, size_t slen,
	const char **err);
/* id < 0 = first valid gene */
[[nodiscard]] bool genome_gene_read(genome_gene *out, const uint8_t *strand, size_t slen,
	int id, const char **err);

/* rewrite a single base (0-based index), in place */
void genome_set_base(uint8_t *strand, size_t n, size_t idx, int digit);


#endif /* GENOME_H */



