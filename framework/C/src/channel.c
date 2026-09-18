/* channel.c -- deterministic synthesis/sequencing error channel */
#include "vivi/channel.h"
#include <string.h>

/* absorbs messages when a caller passes err == nullptr */
static VIVI_THREAD_LOCAL const char *channel_err_sink;

static int digit_at(const uint8_t *s, size_t i)
{
	return (int)((s[i / 4] >> (2 * (3 - i % 4))) & 3u);
}

static bool channel_read_impl(vivi_read *out, const uint8_t *strand, size_t slen,
	const vivi_channel_opts *opts, int soft, const char **err)
{
	memset(out, 0, sizeof(*out));
	if (!strand) { *err = "expected string"; return false; }
	vivi_channel_opts def = { 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0.0 };
	if (!opts) opts = &def;
	const double ps[9] = { opts->p_sub, opts->p_ins, opts->p_del, opts->p_drop,
		opts->p_sub_gc, opts->p_sub_hp, opts->p_trunc, opts->p_burst,
		opts->p_burst_del };
	for (int i = 0; i < 9; i++)
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
	uint8_t *q = soft ? vivi_alloc(2 * n + 4) : nullptr;
	if (soft && !q) { vivi_dealloc(digits); *err = "out of memory"; return false; }
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
			} else if (	vivi_prng_chance(&rng, opts->p_burst)) {
				burst_rem = opts->burst_len ? opts->burst_len : 8;
				bursted = 1;
			}
		}
		if (!	vivi_prng_chance(&rng, opts->p_del)) {
			double p_sub = opts->p_sub;
			int killed = 0;
			if (bursted) {
				if (opts->p_burst_del > 0.0
					&& 	vivi_prng_chance(&rng, opts->p_burst_del))
					killed = 1;
				else
					p_sub = 1.0;
			} else {
				if (opts->p_sub_gc > 0.0 && (d == 1 || d == 2))
					p_sub = 1.0 - (1.0 - p_sub) * (1.0 - opts->p_sub_gc);
				if (opts->p_sub_hp > 0.0 && i > 0 && d == digit_at(strand, i - 1))
					p_sub = 1.0 - (1.0 - p_sub) * (1.0 - opts->p_sub_hp);
			}
			if (!killed) {
				int sub = 	vivi_prng_chance(&rng, p_sub);
				if (sub) {
					int k = (int)(	vivi_prng_next(&rng) % 3u);
					d = (d + 1 + k) % 4;   /* never the original base */
				}
				if (!limit || m < limit) {
					digits[m] = (uint8_t)d;
					if (q) q[m] = (uint8_t)(sub ? VIVI_QUAL_LO : VIVI_QUAL_HI);
					m++;
				}
			}
		}
		if (	vivi_prng_chance(&rng, opts->p_ins)) {
			uint8_t ins = (uint8_t)(	vivi_prng_next(&rng) & 3u);
			if (!limit || m < limit) {
				digits[m] = ins;
				if (q) q[m] = VIVI_QUAL_LO;
				m++;
			}
		}
	}
	size_t bytes = (m + 3) / 4;
	if (bytes) {
		uint8_t *packed = vivi_zalloc(bytes, 1);
		if (!packed) {
			vivi_dealloc(digits);
			vivi_dealloc(q);
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
	if (soft) {
		out->qual = vivi_alloc(m ? m : 1);
		if (!out->qual) {
			vivi_bytes_free(&out->strand);
			vivi_dealloc(q);
			*err = "out of memory";
			return false;
		}
		if (m) memcpy(out->qual, q, m);
	}
	vivi_dealloc(q);
	return true;
}

bool vivi_channel_read(vivi_read *out, const uint8_t *strand, size_t slen,
	const vivi_channel_opts *opts, const char **err)
{
	if (!err) err = &channel_err_sink;
	if (!out) { *err = "expected output"; return false; }
	return channel_read_impl(out, strand, slen, opts, 0, err);
}

bool vivi_channel_read_soft(vivi_read *out, const uint8_t *strand, size_t slen,
	const vivi_channel_opts *opts, const char **err)
{
	if (!err) err = &channel_err_sink;
	if (!out) { *err = "expected output"; return false; }
	return channel_read_impl(out, strand, slen, opts, 1, err);
}

void vivi_read_free(vivi_read *out)
{
	if (!out) return;
	vivi_bytes_free(&out->strand);
	vivi_dealloc(out->qual);
	out->qual = nullptr;
	out->bases = 0;
	out->dropped = 0;
}

bool vivi_read_from_bytes(vivi_read *out, const uint8_t *strand, size_t slen,
	size_t bases, const char **err)
{
	if (!err) err = &channel_err_sink;
	if (out) *out = (vivi_read){ 0 };
	if (!out || (!strand && slen)) { *err = "expected buffers"; return false; }
	if (slen < (bases + 3) / 4) { *err = "strand too short for base count"; return false; }
	if (slen) {
		out->strand.data = vivi_alloc(slen);
		if (!out->strand.data) { *err = "out of memory"; return false; }
		memcpy(out->strand.data, strand, slen);
	}
	out->strand.len = slen;
	out->bases = bases;
	return true;
}

uint8_t vivi_read_base(const vivi_read *r, size_t i)
{
	if (!r || !r->strand.data || i >= r->bases) return 0;
	return (uint8_t)((r->strand.data[i / 4] >> (2 * (3 - i % 4))) & 3u);
}

bool vivi_read_alloc_quality(vivi_read *r, uint8_t fill, const char **err)
{
	if (!err) err = &channel_err_sink;
	if (!r) { *err = "expected read"; return false; }
	if (!r->qual) {
		r->qual = vivi_alloc(r->bases ? r->bases : 1);
		if (!r->qual) { *err = "out of memory"; return false; }
	}
	if (r->bases) memset(r->qual, fill, r->bases);
	return true;
}

bool vivi_read_set_quality(vivi_read *r, size_t i, uint8_t qual)
{
	if (!r || !r->qual || i >= r->bases) return false;
	r->qual[i] = qual;
	return true;
}

bool vivi_consensus_vote(vivi_read *out, const vivi_read *const *reads, size_t count,
	int soft, const char **err)
{
	if (!err) err = &channel_err_sink;
	if (!out) { *err = "expected output"; return false; }
	memset(out, 0, sizeof(*out));
	if (!reads && count) { *err = "expected reads"; return false; }
	size_t n = 0;
	int have = 0;
	for (size_t r = 0; r < count; r++) {
		if (!reads[r] || reads[r]->dropped) continue;
		if (!have) { n = reads[r]->bases; have = 1; }
		else if (reads[r]->bases != n) {
			*err = "consensus requires equal-length reads (no indels)";
			return false;
		}
	}
	if (!have) { out->dropped = 1; return true; }
	uint32_t *counts = vivi_zalloc(n ? n : 1, 4 * sizeof(uint32_t));
	if (!counts) { *err = "out of memory"; return false; }
	for (size_t r = 0; r < count; r++) {
		const vivi_read *rd = reads[r];
		if (!rd || rd->dropped) continue;
		if (rd->bases != n || (n && !rd->strand.data)) {
			vivi_dealloc(counts);
			*err = "consensus requires equal-length reads (no indels)";
			return false;
		}
		for (size_t i = 0; i < n; i++) {
			uint32_t w = 1;
			if (soft && rd->qual) w = (uint32_t)rd->qual[i] + 1u;
			counts[i * 4 + (size_t)vivi_read_base(rd, i)] += w;
		}
	}
	size_t bytes = (n + 3) / 4;
	uint8_t *packed = vivi_zalloc(bytes ? bytes : 1, 1);
	if (!packed) { vivi_dealloc(counts); *err = "out of memory"; return false; }
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

bool vivi_consensus_read(vivi_read *out, const uint8_t *strand, size_t slen,
	const vivi_consensus_opts *opts, const char **err)
{
	if (!err) err = &channel_err_sink;
	if (!out) { *err = "expected output"; return false; }
	if (!opts) return vivi_channel_read(out, strand, slen, nullptr, err);
	int soft = opts->soft != 0;
	if (opts->coverage <= 1)
		return soft ? vivi_channel_read_soft(out, strand, slen, &opts->ch, err)
			: vivi_channel_read(out, strand, slen, &opts->ch, err);
	if (!strand) { *err = "expected string"; return false; }
	if (slen > ((size_t)-1) / 4) { *err = "strand too long"; return false; }
	size_t n = slen * 4;
	if (n == 0)
		return soft ? vivi_channel_read_soft(out, strand, slen, &opts->ch, err)
			: vivi_channel_read(out, strand, slen, &opts->ch, err);
	uint32_t *counts = vivi_zalloc(n, 4 * sizeof(uint32_t));
	if (!counts) { *err = "out of memory"; return false; }
	size_t survivors = 0;
	for (uint32_t r = 0; r < opts->coverage; r++) {
		vivi_channel_opts ch = opts->ch;
		ch.seed = opts->ch.seed + r;
		vivi_read rd;
		if (!(soft ? vivi_channel_read_soft(&rd, strand, slen, &ch, err)
			: vivi_channel_read(&rd, strand, slen, &ch, err))) {
			vivi_dealloc(counts);
			return false;
		}
		if (rd.dropped) { vivi_read_free(&rd); continue; }
		if (rd.bases != n) {
			vivi_read_free(&rd);
			vivi_dealloc(counts);
			*err = "consensus requires equal-length reads (no indels)";
			return false;
		}
		for (size_t i = 0; i < n; i++) {
			uint32_t w = 1;
			if (soft && rd.qual) w = (uint32_t)rd.qual[i] + 1u;
			counts[i * 4 + (size_t)digit_at(rd.strand.data, i)] += w;
		}
		vivi_read_free(&rd);
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
