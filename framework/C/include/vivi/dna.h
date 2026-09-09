/* dna.h -- the quaternary strand codec (C23 port of framework/dna.lua).
 *
 * Binary data <-> packed nucleotide strand over A/C/G/T: one base per
 * 2 bits, 4 bases per byte (00=A 01=C 10=G 11=T, high bits first), so a
 * strand costs ~2% more memory than the payload.
 *
 * Pipeline (encode):
 *   32-bit big-endian payload length header -> header+payload whitened
 *   with xorshift32 keystreams (header: fixed constant seed; payload:
 *   seed derived from the length; deterministic, NOT encryption) ->
 *   homopolymer-constrained coding: 2 bits per base with a 1-bit
 *   fallback alphabet over the opposite base pair when the run reaches
 *   h -> trailing fill/padding bases to keep GC content in
 *   [0.5-gc_eps, 0.5+gc_eps].
 *
 * Byte-exact compatible with the Lua implementation.
 */
#ifndef DNA_H
#define DNA_H

#include "vivi.h"

typedef struct {
	int h;          /* 0 = default 3; otherwise 1..12 */
	double gc_eps;  /* 0 = default 0.05; otherwise [0.005, 0.5) */
} dna_opts;

[[nodiscard]] bool dna_encode(vivi_bytes *out, const uint8_t *data, size_t len,
	const dna_opts *opts, const char **err);
[[nodiscard]] bool dna_decode(vivi_bytes *out, const uint8_t *strand, size_t slen,
	const char **err);

double dna_gc_content(const uint8_t *s, size_t n);
size_t dna_max_homopolymer(const uint8_t *s, size_t n);

[[nodiscard]] bool dna_complement(vivi_bytes *out, const uint8_t *s, size_t n, const char **err);
[[nodiscard]] bool dna_reverse_complement(vivi_bytes *out, const uint8_t *s, size_t n,
	const char **err);

/* nbases < 0 = all bases of the strand */
[[nodiscard]] bool dna_to_ascii(vivi_bytes *out, const uint8_t *s, size_t n, long nbases,
	const char **err);
[[nodiscard]] bool dna_from_ascii(vivi_bytes *out, const char *s, size_t n, const char **err);

[[nodiscard]] bool dna_validate(const uint8_t *s, size_t n, int h, double eps, const char **err);

#endif /* DNA_H */

