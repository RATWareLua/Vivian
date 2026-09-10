/* channel.c -- deterministic synthesis/sequencing error channel */
#include "vivi/channel.h"
#include <string.h>

/* absorbs messages when a caller passes err == nullptr */
static const char *channel_err_sink;

static int digit_at(const uint8_t *s, size_t i)
{
	return (int)((s[i / 4] >> (2 * (3 - i % 4))) & 3u);
}

bool vivi_channel_read(vivi_read *out, const uint8_t *strand, size_t slen,
	const vivi_channel_opts *opts, const char **err)
{
	if (!err) err = &channel_err_sink;
	memset(out, 0, sizeof(*out));
	if (!strand) { *err = "expected string"; return false; }
	vivi_channel_opts def = { 0.0, 0.0, 0.0, 0.0, 0 };
	if (!opts) opts = &def;
	const double ps[4] = { opts->p_sub, opts->p_ins, opts->p_del, opts->p_drop };
	for (int i = 0; i < 4; i++)
		if (!(ps[i] >= 0.0 && ps[i] <= 1.0)) {
			*err = "invalid probability";
			return false;
		}
	if (slen > (((size_t)-1) - 4) / 8) { *err = "strand too long"; return false; }
	vivi_prng rng;
	vivi_prng_init(&rng, opts->seed);
	if (	vivi_prng_chance(&rng, opts->p_drop)) {
		out->dropped = 1;
		return true;
	}
	size_t n = slen * 4;
	uint8_t *digits = vivi_alloc(2 * n + 4);
	if (!digits) { *err = "out of memory"; return false; }
	size_t m = 0;
	for (size_t i = 0; i < n; i++) {
		int d = digit_at(strand, i);
		if (!	vivi_prng_chance(&rng, opts->p_del)) {
			if (	vivi_prng_chance(&rng, opts->p_sub)) {
				int k = (int)(	vivi_prng_next(&rng) % 3u);
				d = (d + 1 + k) % 4;   /* never the original base */
			}
			digits[m++] = (uint8_t)d;
		}
		if (	vivi_prng_chance(&rng, opts->p_ins))
			digits[m++] = (uint8_t)(	vivi_prng_next(&rng) & 3u);
	}
	size_t bytes = (m + 3) / 4;
	if (bytes) {
		uint8_t *packed = vivi_zalloc(bytes, 1);
		if (!packed) {
			vivi_dealloc(digits);
			*err = "out of memory";
			return false;
		}
		for (size_t i = 0; i < m; i++)
			packed[i / 4] |= (uint8_t)(digits[i] << (2 * (3 - i % 4)));
		out->strand.data = packed;
		out->strand.len = bytes;
	}
	out->bases = m;
	vivi_dealloc(digits);
	return true;
}
