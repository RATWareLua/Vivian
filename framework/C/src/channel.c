/* channel.c -- deterministic synthesis/sequencing error channel */
#include "vivi/channel.h"
#include <string.h>

/* absorbs messages when a caller passes err == nullptr */
static VIVI_THREAD_LOCAL const char *channel_err_sink;

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
	vivi_channel_opts def = { 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0 };
	if (!opts) opts = &def;
	const double ps[8] = { opts->p_sub, opts->p_ins, opts->p_del, opts->p_drop,
		opts->p_sub_gc, opts->p_sub_hp, opts->p_trunc, opts->p_burst };
	for (int i = 0; i < 8; i++)
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
	size_t limit = 0;
	if (opts->p_trunc > 0.0 && n >= 2 && vivi_prng_chance(&rng, opts->p_trunc))
		limit = 1 + (size_t)(vivi_prng_next(&rng) % (uint32_t)(n - 1));
	uint8_t *digits = vivi_alloc(2 * n + 4);
	if (!digits) { *err = "out of memory"; return false; }
	size_t m = 0;
	uint32_t burst_rem = 0;
	for (size_t i = 0; i < n; i++) {
		if (limit && m >= limit) break;
		int d = digit_at(strand, i);
		int bursted = 0;
		if (opts->p_burst > 0.0) {
			if (burst_rem > 0) {
				bursted = 1;
				burst_rem--;
			} else if (vivi_prng_chance(&rng, opts->p_burst)) {
				burst_rem = opts->burst_len ? opts->burst_len : 8;
				bursted = 1;
			}
		}
		if (!	vivi_prng_chance(&rng, opts->p_del)) {
			double p_sub = opts->p_sub;
			if (bursted) {
				p_sub = 1.0;
			} else {
				if (opts->p_sub_gc > 0.0 && (d == 1 || d == 2))
					p_sub = 1.0 - (1.0 - p_sub) * (1.0 - opts->p_sub_gc);
				if (opts->p_sub_hp > 0.0 && i > 0 && d == digit_at(strand, i - 1))
					p_sub = 1.0 - (1.0 - p_sub) * (1.0 - opts->p_sub_hp);
			}
			if (	vivi_prng_chance(&rng, p_sub)) {
				int k = (int)(	vivi_prng_next(&rng) % 3u);
				d = (d + 1 + k) % 4;   /* never the original base */
			}
			if (!limit || m < limit) digits[m++] = (uint8_t)d;
		}
		if (	vivi_prng_chance(&rng, opts->p_ins)) {
			uint8_t ins = (uint8_t)(	vivi_prng_next(&rng) & 3u);
			if (!limit || m < limit) digits[m++] = ins;
		}
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

bool vivi_consensus_read(vivi_read *out, const uint8_t *strand, size_t slen,
	const vivi_consensus_opts *opts, const char **err)
{
	if (!err) err = &channel_err_sink;
	if (!out) { *err = "expected output"; return false; }
	if (!opts) return vivi_channel_read(out, strand, slen, nullptr, err);
	if (opts->coverage <= 1) return vivi_channel_read(out, strand, slen, &opts->ch, err);
	if (!strand) { *err = "expected string"; return false; }
	if (slen > ((size_t)-1) / 4) { *err = "strand too long"; return false; }
	size_t n = slen * 4;
	if (n == 0) return vivi_channel_read(out, strand, slen, &opts->ch, err);
	uint32_t *counts = vivi_zalloc(n, 4 * sizeof(uint32_t));
	if (!counts) { *err = "out of memory"; return false; }
	size_t survivors = 0;
	for (uint32_t r = 0; r < opts->coverage; r++) {
		vivi_channel_opts ch = opts->ch;
		ch.seed = opts->ch.seed + r;
		vivi_read rd;
		if (!vivi_channel_read(&rd, strand, slen, &ch, err)) {
			vivi_dealloc(counts);
			return false;
		}
		if (rd.dropped) { vivi_bytes_free(&rd.strand); continue; }
		if (rd.bases != n) {
			vivi_bytes_free(&rd.strand);
			vivi_dealloc(counts);
			*err = "consensus requires equal-length reads (no indels)";
			return false;
		}
		for (size_t i = 0; i < n; i++)
			counts[i * 4 + (size_t)digit_at(rd.strand.data, i)]++;
		vivi_bytes_free(&rd.strand);
		survivors++;
	}
	memset(out, 0, sizeof(*out));
	if (survivors == 0) {
		vivi_dealloc(counts);
		out->dropped = 1;
		return true;
	}
	size_t bytes = (n + 3) / 4;
	uint8_t *packed = vivi_zalloc(bytes, 1);
	if (!packed) {
		vivi_dealloc(counts);
		*err = "out of memory";
		return false;
	}
	for (size_t i = 0; i < n; i++) {
		size_t best = 0;
		for (size_t d = 1; d < 4; d++)
			if (counts[i * 4 + d] > counts[i * 4 + best]) best = d;
		packed[i / 4] |= (uint8_t)(best << (2 * (3 - i % 4)));
	}
	out->strand.data = packed;
	out->strand.len = bytes;
	out->bases = n;
	vivi_dealloc(counts);
	return true;
}
