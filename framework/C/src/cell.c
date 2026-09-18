/* cell.c -- diploidy and repair, transliteration of cell.lua */
#include "vivi/cell.h"
#include "internal.h"
#include <string.h>

static bool cell_valid(const vivi_cell *c)
{
	return c && !c->dead && c->hom[0] && c->hom[1]
		&& c->hlen[0] && c->hlen[1];
}

static bool cell_out_of_memory(const char *err)
{
	return err && strcmp(err, "out of memory") == 0;
}

bool cell_damage_strand(vivi_bytes *out, const uint8_t *strand, size_t slen, int count,
	uint32_t seed, const char **err)
{
	const char *ignored_error = nullptr;
	if (!err) err = &ignored_error;
	if (!out) { *err = "expected output"; return false; }
	*out = (vivi_bytes){ 0 };
	if (!strand) { *err = "expected string"; return false; }
	if (slen > ((size_t)-1) / 4) { *err = "strand too long"; return false; }
	size_t n = slen * 4;
	if (n < 1) {
		out->data = vivi_alloc(1);
		if (!out->data) { *err = "out of memory"; return false; }
		out->len = 0;
		return true;
	}
	if (count <= 0) {
		out->data = vivi_alloc(slen);
		if (!out->data) { *err = "out of memory"; return false; }
		memcpy(out->data, strand, slen);
		out->len = slen;
		return true;
	}
	vivi_rng rng;
	vivi_rng_init(&rng, seed ? seed : 1);
	/* edits are applied straight onto a copy: each draw reads the ORIGINAL
	 * digit, so writing in order gives the same last-write-wins result as
	 * collecting them first -- with no slen*4 scratch array */
	uint8_t *res = vivi_alloc(slen);
	if (!res) { *err = "out of memory"; return false; }
	memcpy(res, strand, slen);
	for (int i = 0; i < count; i++) {
		size_t idx = vivi_rng_next(&rng) % n;
		size_t bytei = idx / 4;
		int j = (int)(idx % 4);
		int sh = 2 * (3 - j);
		int d = (int)((strand[bytei] >> sh) & 3u);
		uint32_t nd = (uint32_t)((d + 1 + (int)(vivi_rng_next(&rng) % 3)) % 4);
		res[bytei] = (uint8_t)((res[bytei] & ~(3u << sh)) | (nd << sh));
	}
	out->data = res;
	out->len = slen;
	return true;
}

bool cell_repair_homolog(vivi_bytes *out, const uint8_t *dst, size_t dlen,
	const uint8_t *src, size_t slen, cell_report *rep, const char **err)
{
	const char *ignored_error = nullptr;
	if (!err) err = &ignored_error;
	if (out) *out = (vivi_bytes){ 0 };
	if (rep) *rep = (cell_report){ 0 };
	if (!out || !rep || !dst || !src) { *err = "invalid repair arguments"; return false; }
	if (dlen != slen) { *err = "homolog length mismatch"; return false; }
	chr_record recS;
	if (!chr_parse(&recS, src, slen, -1, err)) return false;
	if (!recS.cen_ok) {
		rep->dead = 1;
		chr_record_free(&recS);
		/* return a copy of dst unchanged */
		out->data = vivi_alloc(dlen ? dlen : 1);
		if (!out->data) { *err = "out of memory"; chr_record_free(&recS); return false; }
		memcpy(out->data, dst, dlen);
		out->len = dlen;
		return true;
	}
	size_t tb = recS.telomere_bytes;
	uint8_t *cur = vivi_alloc(dlen ? dlen : 1);
	size_t clen = dlen;
	if (!cur) { *err = "out of memory"; chr_record_free(&recS); return false; }
	memcpy(cur, dst, dlen);
	chr_record recD;
	if (!chr_parse(&recD, cur, clen, -1, err)) {
		vivi_dealloc(cur);
		chr_record_free(&recS);
		return false;
	}
	if (!recD.cen_ok) {
		/* centromere destroyed in dst: the head region lives at fixed
		 * positions -- splice it wholesale from the healthy homolog,
		 * with the source layout (primer site + centromere revision) */
		size_t cen_bytes = (recS.cen_version == 2) ? CEN2BYTES : CENBYTES;
		size_t head = tb + (size_t)recS.primer_bytes + cen_bytes;
		uint8_t *nc = vivi_alloc(clen);
		if (!nc) {
			vivi_dealloc(cur);
			chr_record_free(&recD);
			chr_record_free(&recS);
			*err = "out of memory";
			return false;
		}
		memcpy(nc, src, head);
		memcpy(nc + head, cur + head, clen - head);
		vivi_dealloc(cur);
		cur = nc;
		rep->structural = 1;
		chr_record_free(&recD);
		if (!chr_parse(&recD, cur, clen, -1, err)) {
			vivi_dealloc(cur);
			chr_record_free(&recS);
			return false;
		}
	}
	if (!recD.telomere_ok && recS.telomere_ok) {
		int fixed = 0;
		if (memcmp(cur, src, tb) != 0) {
			memcpy(cur, src, tb);
			fixed = 1;
		}
		if (memcmp(cur + clen - tb, src + slen - tb, tb) != 0) {
			memcpy(cur + clen - tb, src + slen - tb, tb);
			fixed = 1;
		}
		if (fixed) {
			rep->structural = 1;
			chr_record_free(&recD);
			if (!chr_parse(&recD, cur, clen, -1, err)) {
				vivi_dealloc(cur);
				chr_record_free(&recS);
				return false;
			}
		}
	}
	if (!recD.cen_ok) {
		rep->dead = 1;
		out->data = cur;
		out->len = clen;
		chr_record_free(&recD);
		chr_record_free(&recS);
		return true;
	}
	/* gene-level repair: splice the healthy src gene at the same locus */
	for (size_t i = 0; i < recS.gene_count; i++) {
		const genome_gene *g = &recS.genes[i];
		const genome_gene *d = nullptr;
		for (size_t j = 0; j < recD.gene_count; j++)
			if (recD.genes[j].id == g->id) { d = &recD.genes[j]; break; }
		int need = 0;
		if (!g->crc_ok) {
			if (!d || !d->crc_ok) rep->dead++;
		} else if (!d) {
			need = 1;
		} else if (!d->crc_ok || d->data_len != g->data_len
			|| (g->data_len > 0 && memcmp(d->data, g->data, g->data_len) != 0)) {
			need = 1;
		}
		if (need && g->crc_ok) {
			/* the true locus always equals the template's offset */
			if (g->offset > dlen || g->size > dlen - g->offset) {
				rep->anomaly++;
				continue;
			}
			memcpy(cur + g->offset, src + g->offset, g->size);
			rep->repaired++;
		}
	}
	if (recD.gene_count > recS.gene_count)
		rep->anomaly += (int)(recD.gene_count - recS.gene_count);
	chr_record_free(&recD);
	chr_record_free(&recS);
	out->data = cur;
	out->len = clen;
	return true;
}

bool cell_checkpoint_pair(vivi_bytes *a_out, vivi_bytes *b_out,
	const uint8_t *h1, size_t l1, const uint8_t *h2, size_t l2,
	cell_report *rep, const char **err)
{
	const char *ignored_error = nullptr;
	if (!err) err = &ignored_error;
	memset(rep, 0, sizeof(*rep));
	vivi_bytes a2, b2, a3;
	cell_report r1, r2, r3;
	if (!cell_repair_homolog(&a2, h1, l1, h2, l2, &r1, err)) return false;
	if (!cell_repair_homolog(&b2, h2, l2, a2.data, a2.len, &r2, err)) {
		vivi_bytes_free(&a2);
		return false;
	}
	if (!cell_repair_homolog(&a3, a2.data, a2.len, b2.data, b2.len, &r3, err)) {
		vivi_bytes_free(&a2);
		vivi_bytes_free(&b2);
		return false;
	}
	rep->repaired = r1.repaired + r2.repaired + r3.repaired;
	rep->dead = r1.dead + r2.dead + r3.dead;
	rep->structural = r1.structural + r2.structural + r3.structural;
	rep->anomaly = r1.anomaly + r2.anomaly + r3.anomaly;
	*a_out = a3;   /* first homolog: healed against the healed second one */
	*b_out = b2;   /* second homolog: healed on the second pass */
	vivi_bytes_free(&a2);
	return true;
}

/* ---- cell lifecycle ---- */

static bool cell_alloc_from_strand(vivi_cell **out, int chr_id, const uint8_t *strand,
	size_t slen)
{
	vivi_cell *c = vivi_zalloc(1, sizeof(vivi_cell));
	if (!c) return false;
	for (int h = 0; h < 2; h++) {
		c->hom[h] = vivi_alloc(slen ? slen : 1);
		if (!c->hom[h]) {
			cell_free(c);
			return false;
		}
		memcpy(c->hom[h], strand, slen);
		c->hlen[h] = slen;
	}
	c->chr_id = chr_id;
	c->generation = 0;
	c->max_gen = 60;
	return *out = c, true;
}

bool cell_new(vivi_cell **out, int chr_id, const uint8_t *data, size_t len,
	const chr_opts *opts, const char **err)
{
	const char *ignored_error = nullptr;
	if (!err) err = &ignored_error;
	*out = nullptr;
	vivi_bytes strand;
	if (!chr_encode(&strand, chr_id, data, len, opts, err)) return false;
	if (!cell_alloc_from_strand(out, chr_id, strand.data, strand.len)) {
		vivi_bytes_free(&strand);
		*err = "out of memory";
		return false;
	}
	vivi_bytes_free(&strand);
	return true;
}

bool cell_stem(vivi_cell **out, int chr_id, const uint8_t *data, size_t len,
	const chr_opts *opts, const char **err)
{
	if (!cell_new(out, chr_id, data, len, opts, err)) return false;
	(*out)->stem = 1;
	return true;
}

void cell_attach_stem(vivi_cell *c, const vivi_cell *stem)
{
	c->stem_source = stem;
}

bool cell_renew(vivi_cell *c, const vivi_cell *stem, const char **err)
{
	const char *ignored_error = nullptr;
	if (!err) err = &ignored_error;
	if (!stem || !stem->hom[0]) { *err = "no stem cell"; return false; }
	chr_record rec;
	if (!chr_parse(&rec, stem->hom[0], stem->hlen[0], -1, err)) return false;
	if (!rec.cen_ok) {
		chr_record_free(&rec);
		*err = "stem cell damaged";
		return false;
	}
	chr_record_free(&rec);
	/* stage both homologs first: a failed allocation leaves c untouched */
	uint8_t *n0 = vivi_alloc(stem->hlen[0] ? stem->hlen[0] : 1);
	if (!n0) { *err = "out of memory"; return false; }
	uint8_t *n1 = vivi_alloc(stem->hlen[1] ? stem->hlen[1] : 1);
	if (!n1) {
		vivi_dealloc(n0);
		*err = "out of memory";
		return false;
	}
	memcpy(n0, stem->hom[0], stem->hlen[0]);
	memcpy(n1, stem->hom[1], stem->hlen[1]);
	vivi_dealloc(c->hom[0]);
	vivi_dealloc(c->hom[1]);
	c->hom[0] = n0;
	c->hlen[0] = stem->hlen[0];
	c->hom[1] = n1;
	c->hlen[1] = stem->hlen[1];
	c->chr_id = stem->chr_id;
	c->generation = 0;
	c->dead = 0;
	return true;
}

void cell_kill(vivi_cell *c)
{
	c->dead = 1;
	for (int h = 0; h < 2; h++) {
		vivi_dealloc(c->hom[h]);
		c->hom[h] = nullptr;
		c->hlen[h] = 0;
	}
}

void cell_free(vivi_cell *c)
{
	if (!c) return;
	for (int h = 0; h < 2; h++) vivi_dealloc(c->hom[h]);
	vivi_dealloc(c);
}

int cell_chr_id(const vivi_cell *c) { return c ? c->chr_id : 0; }
int cell_generation(const vivi_cell *c) { return c ? c->generation : 0; }
int cell_max_generation(const vivi_cell *c) { return c ? c->max_gen : 0; }
void cell_set_max_generation(vivi_cell *c, int max_gen) { if (c) c->max_gen = max_gen; }
bool cell_is_dead(const vivi_cell *c) { return c ? c->dead != 0 : true; }
bool cell_is_stem(const vivi_cell *c) { return c ? c->stem != 0 : false; }

const uint8_t *cell_strand(const vivi_cell *c, size_t homolog, size_t *len)
{
	if (!c || homolog > 1) { if (len) *len = 0; return nullptr; }
	if (len) *len = c->hlen[homolog];
	return c->hom[homolog];
}

uint8_t *cell_strand_mut(vivi_cell *c, size_t homolog, size_t *len)
{
	if (!c || homolog > 1) { if (len) *len = 0; return nullptr; }
	if (len) *len = c->hlen[homolog];
	return c->hom[homolog];
}

bool cell_replace_strand(vivi_cell *c, size_t homolog, uint8_t *data, size_t len)
{
	if (!c || homolog > 1) return false;
	vivi_dealloc(c->hom[homolog]);
	c->hom[homolog] = data;
	c->hlen[homolog] = len;
	return true;
}

bool cell_checkpoint_checked(vivi_cell *c, cell_report *rep, const char **err)
{
	const char *ignored_error = nullptr;
	if (!err) err = &ignored_error;
	if (!rep) { *err = "expected report"; return false; }
	*rep = (cell_report){ 0 };
	if (!cell_valid(c)) {
		rep->failed = 1;
		*err = "invalid cell";
		return false;
	}
	vivi_bytes a, b;
	if (!cell_checkpoint_pair(&a, &b, c->hom[0], c->hlen[0], c->hom[1], c->hlen[1],
		rep, err)) {
		*rep = (cell_report){ .failed = 1 };
		return false;
	}
	vivi_dealloc(c->hom[0]);
	vivi_dealloc(c->hom[1]);
	c->hom[0] = a.data;
	c->hlen[0] = a.len;
	c->hom[1] = b.data;
	c->hlen[1] = b.len;
	return true;
}

cell_report cell_checkpoint(vivi_cell *c)
{
	cell_report rep;
	(void)cell_checkpoint_checked(c, &rep, nullptr);
	return rep;
}

/* strict read: checkpoint + both homologs; no stem fallback */
static int cell_read_checked(vivi_bytes *out, vivi_cell *c, cell_report *rep, const char **err)
{
	if (!cell_checkpoint_checked(c, rep, err)) {
		if (cell_out_of_memory(*err)) { rep->failed = 1; return 0; }
		/* structural damage is a model event, not an execution failure:
		 * fall through so the caller can still try the raw homologs and,
		 * failing those, renew from the stem niche */
		*rep = (cell_report){ 0 };
	}
	for (int h = 0; h < 2; h++) {
		chr_record rec;
		*err = nullptr;
		if (!chr_parse(&rec, c->hom[h], c->hlen[h], -1, err)) {
			if (cell_out_of_memory(*err)) { rep->failed = 1; return 0; }
			continue;
		}
		bool ok = rec.cen_ok && chr_read(out, &rec, err);
		chr_record_free(&rec);
		if (ok) return 1;
		if (cell_out_of_memory(*err)) { rep->failed = 1; return 0; }
	}
	return 0;
}

bool cell_read(vivi_bytes *out, vivi_cell *c, cell_report *rep, const char **err)
{
	const char *ignored_error = nullptr;
	if (!err) err = &ignored_error;
	if (out) *out = (vivi_bytes){ 0 };
	if (rep) *rep = (cell_report){ 0 };
	if (!out || !rep || !c) {
		if (rep) rep->failed = 1;
		*err = "invalid read arguments";
		return false;
	}
	if (!c->dead && c->hom[0]) {
		if (cell_read_checked(out, c, rep, err)) return true;
		if (rep->failed) return false;
	}
	if (c->stem_source) {
		if (!cell_renew(c, c->stem_source, err)) { rep->failed = 1; return false; }
		if (cell_read_checked(out, c, rep, err)) {
			rep->renewed = 1;
			return true;
		}
		if (rep->failed) return false;
	}
	*err = "cell is dead";
	return false;
}

bool cell_peek(vivi_bytes *out, const vivi_cell *c, cell_report *rep, const char **err)
{
	const char *ignored_error = nullptr;
	if (!err) err = &ignored_error;
	if (out) *out = (vivi_bytes){ 0 };
	if (rep) *rep = (cell_report){ 0 };
	if (!out || !rep || !c) {
		if (rep) rep->failed = 1;
		*err = "invalid read arguments";
		return false;
	}
	if (c->dead || !c->hom[0] || !c->hom[1] || !c->hlen[0] || !c->hlen[1]) {
		*err = "cell is dead";
		return false;
	}
	for (int h = 0; h < 2; h++) {
		chr_record rec;
		const char *e = nullptr;
		if (!chr_parse(&rec, c->hom[h], c->hlen[h], -1, &e)) {
			if (cell_out_of_memory(e)) { rep->failed = 1; *err = "out of memory"; return false; }
			continue;
		}
		bool ok = rec.cen_ok && chr_read(out, &rec, &e);
		chr_record_free(&rec);
		if (ok) return true;
		if (cell_out_of_memory(e)) { rep->failed = 1; *err = "out of memory"; return false; }
	}
	*err = "cell damaged";
	return false;
}

static bool cell_next(vivi_cell **out, const vivi_cell *c, const char **err)
{
	if (!cell_valid(c)) { *err = "invalid cell"; return false; }
	if (c->generation < 0 || c->generation >= 65535) {
		*err = "generation out of range";
		return false;
	}
	if (c->generation >= c->max_gen) { *err = "senescent"; return false; }
	vivi_bytes a, b, s;
	cell_report rep;
	if (!cell_checkpoint_pair(&a, &b, c->hom[0], c->hlen[0], c->hom[1], c->hlen[1],
		&rep, err)) return false;
	bool ok = chr_set_generation(&s, a.data, a.len, c->generation + 1, err);
	if (!ok && !cell_out_of_memory(*err))
		ok = chr_set_generation(&s, b.data, b.len, c->generation + 1, err);
	vivi_bytes_free(&a);
	vivi_bytes_free(&b);
	if (!ok) return false;
	ok = cell_alloc_from_strand(out, c->chr_id, s.data, s.len);
	vivi_bytes_free(&s);
	if (!ok) { *err = "out of memory"; return false; }
	(*out)->generation = c->generation + 1;
	(*out)->max_gen = c->max_gen;
	(*out)->stem_source = c->stem_source;
	return true;
}

static void cell_commit(vivi_cell *c, vivi_cell *next)
{
	for (int h = 0; h < 2; h++) {
		vivi_dealloc(c->hom[h]);
		c->hom[h] = next->hom[h];
		c->hlen[h] = next->hlen[h];
		next->hom[h] = nullptr;
	}
	c->generation = next->generation;
}

bool cell_replicate(vivi_cell *c, const char **err)
{
	const char *ignored_error = nullptr;
	if (!err) err = &ignored_error;
	vivi_cell *next = nullptr;
	if (!cell_next(&next, c, err)) return false;
	cell_commit(c, next);
	cell_free(next);
	return true;
}

bool cell_mitosis(vivi_cell **out, vivi_cell *c, const char **err)
{
	const char *ignored_error = nullptr;
	if (!err) err = &ignored_error;
	if (!out) { *err = "expected output"; return false; }
	*out = nullptr;
	vivi_cell *next = nullptr, *daughter = nullptr;
	if (!cell_next(&next, c, err)) return false;
	if (!cell_alloc_from_strand(&daughter, next->chr_id, next->hom[0], next->hlen[0])) {
		cell_free(next);
		*err = "out of memory";
		return false;
	}
	daughter->generation = next->generation;
	daughter->max_gen = next->max_gen;
	daughter->stem_source = next->stem_source;
	cell_commit(c, next);
	cell_free(next);
	*out = daughter;
	return true;
}

bool cell_damage(vivi_cell *c, int count, uint32_t seed, const char **err)
{
	const char *ignored_error = nullptr;
	if (!err) err = &ignored_error;
	if (!cell_valid(c)) { *err = "invalid cell"; return false; }
	vivi_bytes a, b;
	if (!cell_damage_strand(&a, c->hom[0], c->hlen[0], count, seed, err))
		return false;
	if (!cell_damage_strand(&b, c->hom[1], c->hlen[1], count, seed + 1, err)) {
		vivi_bytes_free(&a);
		return false;
	}
	vivi_dealloc(c->hom[0]);
	vivi_dealloc(c->hom[1]);
	c->hom[0] = a.data;
	c->hlen[0] = a.len;
	c->hom[1] = b.data;
	c->hlen[1] = b.len;
	return true;
}

cell_report cell_maintain(vivi_cell **cells, size_t count, const vivi_cell *stem)
{
	cell_report report = { 0 };
	if (!cells && count) { report.failed = 1; return report; }
	for (size_t i = 0; i < count; i++) {
		vivi_cell *c = cells[i];
		if (!c) { report.failed++; continue; }
		const char *err = nullptr;
		cell_report rep = { 0 };
		vivi_bytes data = { 0 };
		bool ok = false;
		if (!c->dead && c->hom[0]) ok = cell_read(&data, c, &rep, &err);
		report.repaired += rep.repaired + rep.structural;
		report.renewed += rep.renewed;
		vivi_bytes_free(&data);
		if (rep.failed) { report.failed++; continue; }
		if (ok && (c->generation < c->max_gen || !stem)) continue;
		if (stem) {
			if (cell_renew(c, stem, &err)) report.renewed++;
			else report.failed++;
		} else if (!ok) {
			cell_kill(c);
			report.dead++;
		}
	}
	return report;
}




