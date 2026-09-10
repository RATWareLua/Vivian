/* pool.h -- oligo pool and primer-based random access (research layer).
 *
 * Random-access DNA storage keeps each oligo separate in a pool and
 * amplifies one address with a specific primer pair. This module models
 * that: a vivi_pool stores gene strands (each a self-contained molecule),
 * and vivi_pool_amplify() returns the damaged amplicon of one gene id,
 * including primer dropout and off-target cross-talk. From
 * VIV14NB4NSH33 the primer sites live on the strand itself; p_primer
 * models their dropout on top of the pool-level p_access.
 */
#ifndef POOL_H
#define POOL_H

#include "vivi.h"
#include "channel.h"

typedef struct {
	int *ids;
	uint8_t **strands;
	size_t *lens;
	size_t count, cap;
} vivi_pool;

void vivi_pool_free(vivi_pool *p);

/* vivi_pool_add_chromosome copies every valid gene of a chromosome;
 * vivi_pool_add_gene adds one molecule under an explicit id */
[[nodiscard]] bool vivi_pool_add_chromosome(vivi_pool *p, const uint8_t *strand, size_t slen,
	const char **err);
[[nodiscard]] bool vivi_pool_add_gene(vivi_pool *p, int id, const uint8_t *strand, size_t slen,
	const char **err);

typedef struct {
	double p_access;       /* probability the primer fails: no product */
	double p_cross;        /* probability the product is an off-target member */
	uint32_t seed;         /* rng seed; 0 behaves as 1 */
	vivi_channel_opts ch;  /* amplification/sequencing errors on the amplicon */
	double p_primer;       /* VIV14NB4NSH33: on-strand primer-site dropout */
} vivi_amp_opts;

typedef struct {
	vivi_read read;        /* damaged amplicon ({ NULL, 0 } when dropped) */
	int id;                /* gene id that was amplified, -1 when dropped */
} vivi_amp_result;

/* Amplifies one gene id from the pool. A primer failure is reported as
 * dropped = 1 with id = -1. An off-target product is reported by id
 * differing from the requested one: the caller must check it. The amplicon
 * channel seed is derived from the amplification seed and the member. */
[[nodiscard]] bool vivi_pool_amplify(vivi_amp_result *out, const vivi_pool *pool, int id,
	const vivi_amp_opts *opts, const char **err);

#endif /* POOL_H */
