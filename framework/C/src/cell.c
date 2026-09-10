/* cell.c -- diploidy and repair, transliteration of cell.lua */
#include "vivi/cell.h"
#include <string.h>

/* absorbs messages when a caller passes err == nullptr */
static const char *cell_err_sink;

bool cell_damage_strand(vivi_bytes *out, const uint8_t *strand, size_t slen, int count,
	uint32_t seed, const char **err)
{
	if (!err) err = &cell_err_sink;
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
	/* collect edits against the ORIGINAL strand (last write wins),
	 * then apply them in one pass -- mirrors the Lua implementation */
	uint8_t *edits = vivi_zalloc(slen * 4, 1); /* per digit: digit+1, 0 = no edit */
	if (!edits) { *err = "out of memory"; return false; }
	for (int i = 0; i < count; i++) {
		size_t idx = vivi_rng_next(&rng) % n;
		size_t bytei = idx / 4;
		int sh = (int)(2 * (3 - idx % 4));
		int d = (int)((strand[bytei] >> sh) & 3u);
		/* stored as digit+1; 0 means "no edit" (any new digit is valid) */
		edits[bytei * 4 + (idx % 4)] =
			(uint8_t)(((d + 1 + (int)(vivi_rng_next(&rng) % 3)) % 4) + 1);
	}
	uint8_t *res = vivi_alloc(slen);
	if (!res) {
		vivi_dealloc(edits);
		*err = "out of memory";
		return false;
	}
	memcpy(res, strand, slen);
	for (size_t b = 0; b < slen; b++) {
		for (int j = 0; j < 4; j++) {
			uint8_t nd = edits[b * 4 + j];
			if (nd) {
				int sh = 2 * (3 - j);
				res[b] = (uint8_t)((res[b] & ~(3u << sh)) | ((uint32_t)(nd - 1) << sh));
			}
		}
	}
	vivi_dealloc(edits);
	out->data = res;
	out->len = slen;
	return true;
}

bool cell_repair_homolog(vivi_bytes *out, const uint8_t *dst, size_t dlen,
	const uint8_t *src, size_t slen, cell_report *rep, const char **err)
{
	if (!err) err = &cell_err_sink;
	memset(out, 0, sizeof(*out));
	memset(rep, 0, sizeof(*rep));
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
	if (!err) err = &cell_err_sink;
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
	if (!err) err = &cell_err_sink;
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
	if (!err) err = &cell_err_sink;
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

cell_report cell_checkpoint(vivi_cell *c)
{
	cell_report rep;
	memset(&rep, 0, sizeof(rep));
	if (c->dead || !c->hom[0]) return rep;
	vivi_bytes a, b;
	const char *err;
	if (!cell_checkpoint_pair(&a, &b, c->hom[0], c->hlen[0], c->hom[1], c->hlen[1],
		&rep, &err))
		return rep;
	vivi_dealloc(c->hom[0]);
	vivi_dealloc(c->hom[1]);
	c->hom[0] = a.data;
	c->hlen[0] = a.len;
	c->hom[1] = b.data;
	c->hlen[1] = b.len;
	return rep;
}

/* strict read: checkpoint + both homologs; no stem fallback */
static int cell_read_checked(vivi_bytes *out, vivi_cell *c, cell_report *rep, const char **err)
{
	*rep = cell_checkpoint(c);
	for (int h = 0; h < 2; h++) {
		chr_record rec;
		if (!chr_parse(&rec, c->hom[h], c->hlen[h], -1, err)) continue;
		if (rec.cen_ok && chr_read(out, &rec, err)) {
			chr_record_free(&rec);
			return 1;
		}
		chr_record_free(&rec);
	}
	return 0;
}

bool cell_read(vivi_bytes *out, vivi_cell *c, cell_report *rep, const char **err)
{
	if (!err) err = &cell_err_sink;
	memset(out, 0, sizeof(*out));
	memset(rep, 0, sizeof(*rep));
	if (c->dead || !c->hom[0]) {
		if (c->stem_source && cell_renew(c, c->stem_source, err)) {
			cell_report r2;
			if (cell_read_checked(out, c, &r2, err)) {
				rep->renewed = 1;
				return true;
			}
		}
		*err = "cell is dead";
		return false;
	}
	if (cell_read_checked(out, c, rep, err)) return true;
	if (c->stem_source && cell_renew(c, c->stem_source, err)) {
		cell_report r2;
		if (cell_read_checked(out, c, &r2, err)) {
			rep->renewed = 1;
			return true;
		}
	}
	*err = "cell is dead";
	return false;
}

bool cell_replicate(vivi_cell *c, const char **err)
{
	if (!err) err = &cell_err_sink;
	if (c->dead || !c->hom[0]) { *err = "cell is dead"; return false; }
	cell_checkpoint(c);
	if (c->generation + 1 > c->max_gen) { *err = "senescent"; return false; }
	vivi_bytes s;
	if (!chr_set_generation(&s, c->hom[0], c->hlen[0], c->generation + 1, err)
		&& !chr_set_generation(&s, c->hom[1], c->hlen[1], c->generation + 1, err)) {
		return false;
	}
	size_t slen = s.len;
	/* stage both copies first: a failed allocation leaves c at the old generation */
	uint8_t *n0 = vivi_alloc(slen ? slen : 1);
	uint8_t *n1 = n0 ? vivi_alloc(slen ? slen : 1) : nullptr;
	if (!n0 || !n1) {
		vivi_dealloc(n0);
		vivi_dealloc(n1);
		vivi_bytes_free(&s);
		*err = "out of memory";
		return false;
	}
	memcpy(n0, s.data, slen);
	memcpy(n1, s.data, slen);
	vivi_bytes_free(&s);
	vivi_dealloc(c->hom[0]);
	vivi_dealloc(c->hom[1]);
	c->hom[0] = n0;
	c->hom[1] = n1;
	c->hlen[0] = c->hlen[1] = slen;
	c->generation++;
	return true;
}

bool cell_mitosis(vivi_cell **out, vivi_cell *c, const char **err)
{
	if (!err) err = &cell_err_sink;
	*out = nullptr;
	if (c->dead || !c->hom[0]) { *err = "cell is dead"; return false; }
	cell_checkpoint(c);
	if (c->generation + 1 > c->max_gen) { *err = "senescent"; return false; }
	vivi_bytes s;
	if (!chr_set_generation(&s, c->hom[0], c->hlen[0], c->generation + 1, err)
		&& !chr_set_generation(&s, c->hom[1], c->hlen[1], c->generation + 1, err)) {
		return false;
	}
	size_t slen = s.len;
	/* stage mother and daughter copies before mutating anything */
	uint8_t *m0 = vivi_alloc(slen ? slen : 1);
	uint8_t *m1 = m0 ? vivi_alloc(slen ? slen : 1) : nullptr;
	uint8_t *d0 = m1 ? vivi_alloc(slen ? slen : 1) : nullptr;
	uint8_t *d1 = d0 ? vivi_alloc(slen ? slen : 1) : nullptr;
	if (!m0 || !m1 || !d0 || !d1) {
		vivi_dealloc(m0);
		vivi_dealloc(m1);
		vivi_dealloc(d0);
		vivi_dealloc(d1);
		vivi_bytes_free(&s);
		*err = "out of memory";
		return false;
	}
	memcpy(m0, s.data, slen);
	memcpy(m1, s.data, slen);
	memcpy(d0, s.data, slen);
	memcpy(d1, s.data, slen);
	vivi_bytes_free(&s);
	vivi_cell *daughter = vivi_zalloc(1, sizeof(vivi_cell));
	if (!daughter) {
		vivi_dealloc(m0);
		vivi_dealloc(m1);
		vivi_dealloc(d0);
		vivi_dealloc(d1);
		*err = "out of memory";
		return false;
	}
	/* the mother ages too and carries the new generation (mirrors Lua) */
	vivi_dealloc(c->hom[0]);
	vivi_dealloc(c->hom[1]);
	c->hom[0] = m0;
	c->hlen[0] = slen;
	c->hom[1] = m1;
	c->hlen[1] = slen;
	c->generation++;
	daughter->chr_id = c->chr_id;
	daughter->hom[0] = d0;
	daughter->hlen[0] = slen;
	daughter->hom[1] = d1;
	daughter->hlen[1] = slen;
	daughter->generation = c->generation;
	daughter->max_gen = c->max_gen;
	daughter->stem_source = c->stem_source;
	*out = daughter;
	return true;
}

bool cell_damage(vivi_cell *c, int count, uint32_t seed, const char **err)
{
	if (!err) err = &cell_err_sink;
	if (c->dead || !c->hom[0]) { *err = "cell is dead"; return false; }
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
	cell_report report;
	memset(&report, 0, sizeof(report));
	for (size_t i = 0; i < count; i++) {
		vivi_cell *c = cells[i];
		if (c->dead || !c->hom[0]) {
			const char *err;
			if (stem && cell_renew(c, stem, &err)) report.renewed++;
			else {
				cell_kill(c);
				report.dead++;
			}
		} else {
			cell_report rep = cell_checkpoint(c);
			report.repaired += rep.repaired + rep.structural;
			vivi_bytes data;
			const char *err;
			int ok = cell_read(&data, c, &rep, &err);
			if (ok) {
				vivi_bytes_free(&data);
				if (c->generation >= c->max_gen && stem) {
					(void)cell_renew(c, stem, &err);
					report.renewed++;
				}
			} else if (stem && cell_renew(c, stem, &err)) {
				report.renewed++;
			} else {
				cell_kill(c);
				report.dead++;
			}
		}
	}
	return report;
}




