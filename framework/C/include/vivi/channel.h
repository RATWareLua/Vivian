/* channel.h -- synthesis/sequencing error channel (research layer).
 *
 * The channel sits on top of the format: it takes a packed strand
 * (4 bases per byte, high bits first, exactly what the codec produces)
 * and returns one "read" of it -- a base sequence with substitutions,
 * insertions and deletions applied, or nothing at all when the read is
 * dropped. Everything is deterministic per seed, so an experiment is
 * reproducible from its parameters alone.
 *
 * The damaged strand is repacked into whole bytes; `bases` reports the
 * exact base count (0..3 trailing filler bases are zero/A). A dropped
 * read comes back as `dropped = 1` with an empty strand.
 *
 * Draw order per original base (pinned for reproducibility):
 *   1. delete?        then skip the base
 *   2. substitute?    then draw the replacement digit (never the same one)
 *   3. insert?        then draw the inserted digit (appended after)
 */
#ifndef CHANNEL_H
#define CHANNEL_H

#include "vivi.h"

typedef struct {
	double p_sub;   /* substitution probability per original base, [0, 1] */
	double p_ins;   /* insertion probability after each original base, [0, 1] */
	double p_del;   /* deletion probability per original base, [0, 1] */
	double p_drop;  /* probability that the whole read is lost, [0, 1] */
	uint32_t seed;  /* rng seed; 0 behaves as 1 */
	double p_sub_gc;    /* EXTRA substitution probability when the ORIGINAL base is G or C */
	double p_sub_hp;    /* EXTRA substitution probability when the ORIGINAL base equals the previous original base (homopolymer) */
	double p_trunc;     /* probability the read is truncated at a random base */
	double p_burst;     /* per-base probability of STARTING a correlated burst of substitutions */
	uint32_t burst_len; /* number of original bases substituted at the start of a burst (0 = 8) */
} vivi_channel_opts;

typedef struct {
	vivi_bytes strand;  /* damaged packed strand ({ NULL, 0 } when dropped) */
	size_t bases;       /* exact base count of the read */
	int dropped;        /* 1 = the read was lost */
} vivi_read;

/* opts may be NULL (an identity channel). */
[[nodiscard]] bool vivi_channel_read(vivi_read *out, const uint8_t *strand, size_t slen,
	const vivi_channel_opts *opts, const char **err);

/* Coverage: read the same molecule `coverage` times and majority-vote
 * every base, the way a real pipeline sequences many copies of one oligo.
 * Substitution-only: with coverage > 1 every surviving read must keep the
 * original base count, so p_ins and p_del must stay 0. A per-base tie
 * picks the lowest digit, so the result is fully deterministic. Reads are
 * seeded `ch.seed + r` for r = 0..coverage-1; coverage <= 1 is exactly
 * vivi_channel_read with `ch`, so single-read streams are unchanged. */
typedef struct {
	vivi_channel_opts ch;  /* per-read channel; seed is the base read seed */
	uint32_t coverage;     /* independent reads, 1 = a plain read */
} vivi_consensus_opts;

[[nodiscard]] bool vivi_consensus_read(vivi_read *out, const uint8_t *strand, size_t slen,
	const vivi_consensus_opts *opts, const char **err);

#endif /* CHANNEL_H */
