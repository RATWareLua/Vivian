/* rs.h -- systematic Reed-Solomon over GF(256), byte-oriented (research layer).
 *
 * A short inner code for one molecule: k data bytes plus m parity bytes
 * (n = k + m <= 255), correcting up to floor(m/2) unknown byte errors.
 * Unlike vivi/parity.h, which recovers whole shards whose positions are
 * known, this decoder locates the errors itself.
 */
#ifndef RS_H
#define RS_H

#include "vivi.h"

#define RS_MAX_N 255

/* data may be NULL when k == 0. out receives n = k + m bytes. */
[[nodiscard]] bool rs_encode(vivi_bytes *out, const uint8_t *data, size_t k, size_t m,
	const char **err);

/* cw holds n = k + m bytes and is corrected in place; the first k bytes
 * are the data. *corrected receives the number of byte errors fixed. */
[[nodiscard]] bool rs_decode(uint8_t *cw, size_t n, size_t k, size_t *corrected,
	const char **err);

#endif /* RS_H */
