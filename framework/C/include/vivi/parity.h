/* parity.h -- systematic Reed-Solomon erasure coding (research layer).
 *
 * Splits a payload-agnostic byte stream into n data shards and m parity
 * shards over GF(256): the code is systematic (shards 0..n-1 are the data
 * verbatim) and MDS, so any n of the n+m shards reconstruct the data.
 * Used to model an outer erasure code across gene molecules: with m parity
 * molecules the library survives any m lost/damaged genes.
 *
 * Shards must all be the same length (pad shorter data at the caller's
 * side). n + m <= 255 because the code evaluates at the 255 nonzero field
 * elements.
 */
#ifndef PARITY_H
#define PARITY_H

#include "vivi.h"

/* Encodes n data shards of `len` bytes into n + m shards. The caller owns
 * the `shards` array; on success every entry holds an allocated buffer, on
 * failure the function releases whatever it allocated and nulls the
 * entries. Release with vivi_parity_release(). */
[[nodiscard]] bool vivi_parity_encode(uint8_t **shards, const uint8_t *const *data,
	size_t n, size_t m, size_t len, const char **err);

/* Reconstructs the n data shards. `shards` holds n + m read buffers,
 * `present[k]` is nonzero when shard k is usable (at least n must be).
 * `out` is caller-owned; on success every entry holds a fresh buffer.
 * Fails with "too few shards" when fewer than n are present. */
[[nodiscard]] bool vivi_parity_decode(uint8_t **out, const uint8_t *const *shards,
	const uint8_t *present, size_t n, size_t m, size_t len, const char **err);

/* Frees `count` shard buffers (null entries are ignored). */
void vivi_parity_release(uint8_t **shards, size_t count);

#endif /* PARITY_H */
