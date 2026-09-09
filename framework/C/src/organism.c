/* organism.c -- organisms, transliteration of organism.lua */
#include "vivi/organism.h"
#include <string.h>

static const char ORG_MAGIC[4] = { 'G', 'S', 'Y', '1' };

static void organism_free_arrays(vivi_organism *o)
{
	for (int h = 0; h < 2; h++) {
		if (o->hom[h]) {
			for (size_t i = 0; i < o->nchr; i++) vivi_dealloc(o->hom[h][i]);
			vivi_dealloc(o->hom[h]);
			vivi_dealloc(o->hlen[h]);
			o->hom[h] = nullptr;
			o->hlen[h] = nullptr;
		}
	}
}

void organism_free(vivi_organism *o)
{
	if (!o) return;
	organism_free_arrays(o);
	vivi_dealloc(o->chr_ids);
	vivi_dealloc(o->chr_opts);
	vivi_dealloc(o);
}

static int organism_alloc_from_genome(vivi_organism **out, const int *ids,
	const chr_opts *opts, uint8_t **strands, const size_t *lens, size_t nchr,
	int max_gen)
{
	vivi_organism *o = vivi_zalloc(1, sizeof(vivi_organism));
	if (!o) return 0;
	o->nchr = nchr;
	o->chr_ids = vivi_alloc(nchr * sizeof(int));
	o->chr_opts = vivi_alloc(nchr * sizeof(chr_opts));
	for (int h = 0; h < 2; h++) {
		o->hom[h] = vivi_alloc(nchr * sizeof(uint8_t *));
		o->hlen[h] = vivi_alloc(nchr * sizeof(size_t));
	}
	if (!o->chr_ids || !o->chr_opts || !o->hom[0] || !o->hom[1]
		|| !o->hlen[0] || !o->hlen[1]) {
		organism_free(o);
		return 0;
	}
	for (size_t i = 0; i < nchr; i++) {
		o->chr_ids[i] = ids[i];
		o->chr_opts[i] = opts[i];
		for (int h = 0; h < 2; h++) {
			o->hom[h][i] = vivi_alloc(lens[i] ? lens[i] : 1);
			if (!o->hom[h][i]) {
				organism_free(o);
				return 0;
			}
			memcpy(o->hom[h][i], strands[i], lens[i]);
			o->hlen[h][i] = lens[i];
		}
	}
	o->generation = 0;
	o->max_gen = max_gen;
	*out = o;
	return 1;
}

bool organism_new(vivi_organism **out, const int *ids, const chr_opts *opts,
	const uint8_t *const *datas, const size_t *lens, size_t nchr,
	int max_gen, const char **err)
{
	*out = nullptr;
	if (!ids || nchr == 0) { *err = "expected non-empty specs array"; return false; }
	if (nchr > 255) { *err = "too many chromosomes (max 255)"; return false; }
	uint8_t **strands = vivi_alloc(nchr * sizeof(uint8_t *));
	size_t *slens = vivi_alloc(nchr * sizeof(size_t));
	if (!strands || !slens) {
		vivi_dealloc(strands);
		vivi_dealloc(slens);
		*err = "out of memory";
		return false;
	}
	int ok = 1;
	for (size_t i = 0; i < nchr && ok; i++) {
		if (ids[i] < 0 || ids[i] > 255) { *err = "invalid chr_id"; ok = 0; break; }
		for (size_t j = 0; j < i; j++) {
			if (ids[j] == ids[i]) { *err = "duplicate chr_id"; ok = 0; break; }
		}
		if (!ok) break;
		vivi_bytes strand;
		if (!chr_encode(&strand, ids[i], datas[i], lens[i], &opts[i], err)) {
			ok = 0;
			break;
		}
		strands[i] = strand.data;
		slens[i] = strand.len;
	}
	if (ok && !organism_alloc_from_genome(out, ids, opts, strands, slens, nchr,
		max_gen ? max_gen : 60)) {
		*err = "out of memory";
		ok = 0;
	}
	for (size_t i = 0; i < nchr; i++) vivi_dealloc(strands[i]);
	vivi_dealloc(strands);
	vivi_dealloc(slens);
	return ok;
}

bool organism_stem(vivi_organism **out, const int *ids, const chr_opts *opts,
	const uint8_t *const *datas, const size_t *lens, size_t nchr,
	int max_gen, const char **err)
{
	if (!organism_new(out, ids, opts, datas, lens, nchr, max_gen, err)) return false;
	(*out)->stem = 1;
	return true;
}

void organism_attach_stem(vivi_organism *o, const vivi_organism *stem)
{
	o->stem_source = stem;
}

bool organism_renew(vivi_organism *o, const vivi_organism *stem, const char **err)
{
	if (!stem || !stem->hom[0]) { *err = "no stem organism"; return false; }
	size_t nchr = stem->nchr;
	uint8_t **g1 = vivi_zalloc(nchr, sizeof(uint8_t *));
	uint8_t **g2 = vivi_zalloc(nchr, sizeof(uint8_t *));
	size_t *l1 = vivi_zalloc(nchr, sizeof(size_t));
	size_t *l2 = vivi_zalloc(nchr, sizeof(size_t));
	if (!g1 || !g2 || !l1 || !l2) {
		vivi_dealloc(g1); vivi_dealloc(g2); vivi_dealloc(l1); vivi_dealloc(l2);
		*err = "out of memory";
		return false;
	}
	int fail = 0;
	for (size_t i = 0; i < nchr && !fail; i++) {
		l1[i] = l2[i] = stem->hlen[0][i];
		g1[i] = vivi_alloc(l1[i] ? l1[i] : 1);
		g2[i] = vivi_alloc(l2[i] ? l2[i] : 1);
		if (!g1[i] || !g2[i]) fail = 1;
		else {
			memcpy(g1[i], stem->hom[0][i], l1[i]);
			memcpy(g2[i], stem->hom[1][i], l2[i]);
		}
	}
	if (fail) {
		for (size_t i = 0; i < nchr; i++) { vivi_dealloc(g1[i]); vivi_dealloc(g2[i]); }
		vivi_dealloc(g1); vivi_dealloc(g2); vivi_dealloc(l1); vivi_dealloc(l2);
		*err = "out of memory";
		return false;
	}
	organism_free_arrays(o);
	o->hom[0] = g1;
	o->hom[1] = g2;
	o->hlen[0] = l1;
	o->hlen[1] = l2;
	o->nchr = nchr;
	int *ids = vivi_zalloc(nchr, sizeof(int));
	chr_opts *op = vivi_alloc(nchr * sizeof(chr_opts));
	if (!ids || !op) {
		vivi_dealloc(ids);
		vivi_dealloc(op);
		*err = "out of memory";
		return false;
	}
	for (size_t i = 0; i < nchr; i++) {
		ids[i] = stem->chr_ids[i];
		op[i] = stem->chr_opts[i];
	}
	vivi_dealloc(o->chr_ids);
	vivi_dealloc(o->chr_opts);
	o->chr_ids = ids;
	o->chr_opts = op;
	o->generation = 0;
	o->dead = 0;
	return true;
}

void organism_kill(vivi_organism *o)
{
	o->dead = 1;
	organism_free_arrays(o);
}

cell_report organism_checkpoint(vivi_organism *o)
{
	cell_report rep;
	memset(&rep, 0, sizeof(rep));
	if (o->dead || !o->hom[0]) return rep;
	for (size_t i = 0; i < o->nchr; i++) {
		vivi_bytes a, b;
		cell_report r;
		const char *err;
		if (!cell_checkpoint_pair(&a, &b, o->hom[0][i], o->hlen[0][i],
			o->hom[1][i], o->hlen[1][i], &r, &err))
			continue;
		vivi_dealloc(o->hom[0][i]);
		vivi_dealloc(o->hom[1][i]);
		o->hom[0][i] = a.data;
		o->hlen[0][i] = a.len;
		o->hom[1][i] = b.data;
		o->hlen[1][i] = b.len;
		rep.repaired += r.repaired;
		rep.dead += r.dead;
		rep.structural += r.structural;
		rep.anomaly += r.anomaly;
	}
	return rep;
}

/* checkpoint + strict read (every chromosome must be readable);
 * out is an array of o->nchr vivi_bytes; no stem fallback */
static void org_karyotype_clear(vivi_bytes *arr, size_t n)
{
	for (size_t i = 0; i < n; i++) vivi_bytes_free(&arr[i]);
	vivi_dealloc(arr);
}

static int organism_read_checked(vivi_bytes *out, vivi_organism *o, cell_report *rep,
	int *broken_chr, const char **err)
{
	*rep = organism_checkpoint(o);
	int all_ok = 1;
	for (size_t i = 0; i < o->nchr; i++) {
		vivi_bytes_free(&out[i]); /* safe on zeroed entries; needed on renew retry */
		out[i] = (vivi_bytes){ 0 };
		int good = 0;
		for (int h = 0; h < 2; h++) {
			chr_record rec;
			if (!chr_parse(&rec, o->hom[h][i], o->hlen[h][i], -1, err)) continue;
			if (rec.cen_ok && chr_read(&out[i], &rec, err)) good = 1;
			chr_record_free(&rec);
			if (good) break;
		}
		if (!good) {
			*broken_chr = o->chr_ids[i];
			all_ok = 0;
		}
	}
	return all_ok;
}

bool organism_read(vivi_bytes **out, vivi_organism *o, cell_report *rep, const char **err)
{
	memset(rep, 0, sizeof(*rep));
	vivi_bytes *tmp = vivi_zalloc(o->nchr ? o->nchr : 1, sizeof(vivi_bytes));
	if (!tmp) { *err = "out of memory"; return false; }
	int broken;
	cell_report r;
	if (o->dead || !o->hom[0]) {
		if (o->stem_source && organism_renew(o, o->stem_source, err)) {
			if (organism_read_checked(tmp, o, &r, &broken, err)) {
				*out = tmp;

				rep->renewed = 1;
				return true;
			}
		}
		org_karyotype_clear(tmp, o->nchr ? o->nchr : 1);
		*err = "organism dead";
		return false;
	}
	if (organism_read_checked(tmp, o, &r, &broken, err)) {
		*out = tmp;

		*rep = r;
		return true;
	}
	if (o->stem_source && organism_renew(o, o->stem_source, err)) {
		if (organism_read_checked(tmp, o, &r, &broken, err)) {
			*out = tmp;

			rep->renewed = 1;
			return true;
		}
	}
	org_karyotype_clear(tmp, o->nchr ? o->nchr : 1);
	*err = "chromosome dead";
	return false;
}

static int org_set_generation_all(vivi_organism *o, int generation, const char **err)
{
	for (size_t i = 0; i < o->nchr; i++) {
		vivi_bytes s;
		if (!chr_set_generation(&s, o->hom[0][i], o->hlen[0][i], generation, err)
			&& !chr_set_generation(&s, o->hom[1][i], o->hlen[1][i], generation, err)) {
			return 0;
		}
		vivi_dealloc(o->hom[0][i]);
		vivi_dealloc(o->hom[1][i]);
		o->hom[0][i] = vivi_alloc(s.len ? s.len : 1);
		o->hom[1][i] = vivi_alloc(s.len ? s.len : 1);
		if (!o->hom[0][i] || !o->hom[1][i]) {
			vivi_bytes_free(&s);
			*err = "out of memory";
			return 0;
		}
		memcpy(o->hom[0][i], s.data, s.len);
		memcpy(o->hom[1][i], s.data, s.len);
		o->hlen[0][i] = o->hlen[1][i] = s.len;
		vivi_bytes_free(&s);
	}
	return 1;
}

bool organism_replicate(vivi_organism *o, const char **err)
{
	if (o->dead || !o->hom[0]) { *err = "organism dead"; return false; }
	organism_checkpoint(o);
	if (o->generation + 1 > o->max_gen) { *err = "senescent"; return false; }
	o->generation++;
	return org_set_generation_all(o, o->generation, err);
}

bool organism_mitosis(vivi_organism **out, vivi_organism *o, const char **err)
{
	*out = nullptr;
	if (o->dead || !o->hom[0]) { *err = "organism dead"; return false; }
	organism_checkpoint(o);
	if (o->generation + 1 > o->max_gen) { *err = "senescent"; return false; }
	o->generation++;
	if (!org_set_generation_all(o, o->generation, err)) return false;
	vivi_organism *d = vivi_zalloc(1, sizeof(vivi_organism));
	if (!d) { *err = "out of memory"; return false; }
	d->nchr = o->nchr;
	d->chr_ids = vivi_alloc(o->nchr * sizeof(int));
	d->chr_opts = vivi_alloc(o->nchr * sizeof(chr_opts));
	for (int h = 0; h < 2; h++) {
		d->hom[h] = vivi_alloc(o->nchr * sizeof(uint8_t *));
		d->hlen[h] = vivi_alloc(o->nchr * sizeof(size_t));
	}
	if (!d->chr_ids || !d->chr_opts || !d->hom[0] || !d->hom[1]
		|| !d->hlen[0] || !d->hlen[1]) {
		organism_free(d);
		*err = "out of memory";
		return false;
	}
	for (size_t i = 0; i < o->nchr; i++) {
		d->chr_ids[i] = o->chr_ids[i];
		d->chr_opts[i] = o->chr_opts[i];
		for (int h = 0; h < 2; h++) {
			d->hom[h][i] = vivi_alloc(o->hlen[h][i] ? o->hlen[h][i] : 1);
			if (!d->hom[h][i]) {
				organism_free(d);
				*err = "out of memory";
				return false;
			}
			memcpy(d->hom[h][i], o->hom[h][i], o->hlen[h][i]);
			d->hlen[h][i] = o->hlen[h][i];
		}
	}
	d->generation = o->generation;
	d->max_gen = o->max_gen;
	d->stem_source = o->stem_source;
	*out = d;
	return true;
}

bool organism_damage(vivi_organism *o, int count, uint32_t seed, const char **err)
{
	if (o->dead || !o->hom[0]) { *err = "organism dead"; return false; }
	for (size_t i = 0; i < o->nchr; i++) {
		vivi_bytes a, b;
		if (!cell_damage_strand(&a, o->hom[0][i], o->hlen[0][i], count,
			seed + (uint32_t)i + 1, err))
			return false;
		if (!cell_damage_strand(&b, o->hom[1][i], o->hlen[1][i], count,
			seed + 1000u + (uint32_t)i + 1, err)) {
			vivi_bytes_free(&a);
			return false;
		}
		vivi_dealloc(o->hom[0][i]);
		vivi_dealloc(o->hom[1][i]);
		o->hom[0][i] = a.data;
		o->hlen[0][i] = a.len;
		o->hom[1][i] = b.data;
		o->hlen[1][i] = b.len;
	}
	return true;
}

bool organism_mutate(org_mut_result *out, vivi_organism *o, int count, uint32_t seed,
	const char **err)
{
	memset(out, 0, sizeof(*out));
	if (o->dead || !o->hom[0]) { *err = "organism dead"; return false; }
	vivi_rng rng;
	vivi_rng_init(&rng, seed);
	size_t ci = vivi_rng_next(&rng) % o->nchr;
	int cid = o->chr_ids[ci];
	vivi_bytes data = { 0 };
	for (int h = 0; h < 2; h++) {
		chr_record rec;
		if (!chr_parse(&rec, o->hom[h][ci], o->hlen[h][ci], -1, err)) continue;
		if (rec.cen_ok && chr_read(&data, &rec, err)) {
			chr_record_free(&rec);
			break;
		}
		chr_record_free(&rec);
	}
	if (!data.data) { *err = "chromosome unreadable"; return false; }
	if (data.len == 0 || count <= 0) {
		out->chr = cid;
		out->data = data; /* unchanged */
		return true;
	}
	uint8_t *t = vivi_alloc(data.len);
	if (!t) {
		vivi_bytes_free(&data);
		*err = "out of memory";
		return false;
	}
	memcpy(t, data.data, data.len);
	for (int i = 0; i < count; i++) {
		size_t bi = vivi_rng_next(&rng) % data.len;
		int bit = (int)(vivi_rng_next(&rng) % 8);
		t[bi] = (uint8_t)(t[bi] ^ (1u << bit));
	}
	size_t mlen = data.len;
	vivi_bytes_free(&data);
	vivi_bytes strand;
	if (!chr_encode(&strand, cid, t, mlen, &o->chr_opts[ci], err)) {
		vivi_dealloc(t);
		return false;
	}
	vivi_dealloc(o->hom[0][ci]);
	vivi_dealloc(o->hom[1][ci]);
	o->hom[0][ci] = vivi_alloc(strand.len ? strand.len : 1);
	o->hom[1][ci] = vivi_alloc(strand.len ? strand.len : 1);
	if (!o->hom[0][ci] || !o->hom[1][ci]) {
		vivi_bytes_free(&strand);
		vivi_dealloc(t);
		*err = "out of memory";
		return false;
	}
	memcpy(o->hom[0][ci], strand.data, strand.len);
	memcpy(o->hom[1][ci], strand.data, strand.len);
	o->hlen[0][ci] = o->hlen[1][ci] = strand.len;
	vivi_bytes_free(&strand);
	out->chr = cid;
	out->data.data = t;
	out->data.len = mlen;
	return true;
}

bool organism_cross(vivi_organism **out, vivi_organism *pa, vivi_organism *pb,
	uint32_t seed, const char **err)
{
	*out = nullptr;
	if (!pa || !pb || !pa->hom[0] || !pb->hom[0]) { *err = "expected organisms"; return false; }
	if (pa->nchr != pb->nchr) { *err = "incompatible genomes"; return false; }
	for (size_t i = 0; i < pa->nchr; i++) {
		if (pa->chr_ids[i] != pb->chr_ids[i]) { *err = "incompatible genomes"; return false; }
		const chr_opts *a = &pa->chr_opts[i], *b = &pb->chr_opts[i];
		if ((a->gene_raw ? a->gene_raw : 1024) != (b->gene_raw ? b->gene_raw : 1024)
			|| a->codon != b->codon
			|| (a->h ? a->h : 3) != (b->h ? b->h : 3)
			|| (a->units ? a->units : 4) != (b->units ? b->units : 4)
			|| a->flags != b->flags) {
			*err = "incompatible chromosome opts";
			return false;
		}
	}
	vivi_rng rng;
	vivi_rng_init(&rng, seed ? seed : 1);
	size_t nchr = pa->nchr;
	int *ids = vivi_zalloc(nchr, sizeof(int));
	chr_opts *opts_t = vivi_zalloc(nchr, sizeof(chr_opts));
	uint8_t **g1 = vivi_zalloc(nchr, sizeof(uint8_t *));
	uint8_t **g2 = vivi_zalloc(nchr, sizeof(uint8_t *));
	size_t *l1 = vivi_zalloc(nchr, sizeof(size_t));
	size_t *l2 = vivi_zalloc(nchr, sizeof(size_t));
	if (!ids || !opts_t || !g1 || !g2 || !l1 || !l2) {
		vivi_dealloc(ids); vivi_dealloc(opts_t); vivi_dealloc(g1); vivi_dealloc(g2); vivi_dealloc(l1); vivi_dealloc(l2);
		*err = "out of memory";
		return false;
	}
	int ok = 1;
	for (size_t i = 0; i < nchr && ok; i++) {
		/* per parent: healthy gene data, first valid per id across homologs */
		vivi_bytes *ga = nullptr, *gb = nullptr;
		int ngenesA = -1, ngenesB = -1;
		for (int h = 0; h < 2; h++) {
			chr_record rec;
			if (!chr_parse(&rec, pa->hom[h][i], pa->hlen[h][i], -1, err)) continue;
			if (rec.cen_ok && ngenesA < 0) ngenesA = rec.ngenes;
			chr_record_free(&rec);
		}
		for (int h = 0; h < 2; h++) {
			chr_record recB;
			if (!chr_parse(&recB, pb->hom[h][i], pb->hlen[h][i], -1, err)) continue;
			if (recB.cen_ok && ngenesB < 0) ngenesB = recB.ngenes;
			chr_record_free(&recB);
		}
		if (ngenesA < 0 || ngenesB < 0) { *err = "parent chromosome broken"; ok = 0; }
		else if (ngenesA != ngenesB) {
			*err = "incompatible gene counts on chr";
			ok = 0;
		}
		if (ok && ngenesA > 0) {
			ga = vivi_zalloc((size_t)ngenesA, sizeof(vivi_bytes));
			gb = vivi_zalloc((size_t)ngenesA, sizeof(vivi_bytes));
			if (!ga || !gb) ok = 0;
		}
		if (ok) {
			for (int hh = 0; hh < 2; hh++) {
				const vivi_organism *parent = hh ? pb : pa;
				chr_record rec;
				if (!chr_parse(&rec, parent->hom[hh][i], parent->hlen[hh][i], -1, err))
					continue;
				if (rec.cen_ok) {
					for (size_t kk = 0; kk < rec.gene_count; kk++) {
						const genome_gene *g = &rec.genes[kk];
						vivi_bytes *slot = hh ? &gb[g->id] : &ga[g->id];
						if (g->crc_ok && g->id < ngenesA && !slot->data) {
							slot->data = vivi_alloc(g->data_len ? g->data_len : 1);
							if (slot->data) {
								memcpy(slot->data, g->data, g->data_len);
								slot->len = g->data_len;
							}
						}
					}
				}
				chr_record_free(&rec);
			}
			/* one balanced recombinant mosaic, shared by both homologs */
			vivi_buf parts = { 0 };
			for (int gid = 0; gid < ngenesA && ok; gid++) {
				int useB = (int)(vivi_rng_next(&rng) % 2);
				const vivi_bytes *d;
				if (useB) d = gb[gid].data ? &gb[gid] : &ga[gid];
				else d = ga[gid].data ? &ga[gid] : &gb[gid];
				if (!d->data) {
					*err = "parent gene missing";
					ok = 0;
					break;
				}
				if (!vivi_buf_append(&parts, d->data, d->len)) {
					*err = "out of memory";
					ok = 0;
					break;
				}
			}
			if (ok) {
				vivi_bytes strand;
				if (!chr_encode(&strand, pa->chr_ids[i], parts.p, parts.len,
					&pa->chr_opts[i], err)) {
					ok = 0;
				} else {
					l1[i] = l2[i] = strand.len;
					g1[i] = vivi_alloc(strand.len ? strand.len : 1);
					if (!g1[i]) ok = 0;
					else {
						memcpy(g1[i], strand.data, strand.len);
						/* the child is homozygous: both homologs are the
						 * same balanced mosaic (fresh independent copies) */
						g2[i] = vivi_alloc(strand.len ? strand.len : 1);
						if (!g2[i]) ok = 0;
						else memcpy(g2[i], strand.data, strand.len);
					}
					if (!ok) {
						vivi_dealloc(g1[i]);
						vivi_dealloc(strand.data);
						g1[i] = nullptr;
						g2[i] = nullptr;
					} else {
						strand.data = nullptr;
					}
				}
			}
			vivi_buf_free(&parts);
		}
		if (ga) {
			for (int gid = 0; gid < ngenesA; gid++) vivi_bytes_free(&ga[gid]);
			vivi_dealloc(ga);
		}
		if (gb) {
			for (int gid = 0; gid < ngenesB; gid++) vivi_bytes_free(&gb[gid]);
			vivi_dealloc(gb);
		}
		if (ok) {
			ids[i] = pa->chr_ids[i];
			opts_t[i] = pa->chr_opts[i];
		}
	}
	if (!ok) {
		for (size_t i = 0; i < nchr; i++) { vivi_dealloc(g1[i]); vivi_dealloc(g2[i]); }
		vivi_dealloc(ids);
		vivi_dealloc(opts_t);
		vivi_dealloc(g1);
		vivi_dealloc(g2);
		vivi_dealloc(l1);
		vivi_dealloc(l2);
		return false;
	}
	vivi_organism *o = vivi_zalloc(1, sizeof(vivi_organism));
	if (!o) {
		for (size_t i = 0; i < nchr; i++) { vivi_dealloc(g1[i]); vivi_dealloc(g2[i]); }
		vivi_dealloc(ids);
		vivi_dealloc(opts_t);
		vivi_dealloc(g1);
		vivi_dealloc(g2);
		vivi_dealloc(l1);
		vivi_dealloc(l2);
		*err = "out of memory";
		return false;
	}
	o->nchr = nchr;
	o->chr_ids = ids;
	o->chr_opts = opts_t;
	o->hom[0] = g1;
	o->hom[1] = g2;
	o->hlen[0] = l1;
	o->hlen[1] = l2;
	o->generation = 0;
	o->max_gen = pa->max_gen;
	*out = o;
	return true;
}

bool organism_serialize(vivi_bytes *out, vivi_organism *o, const char **err)
{
	*out = (vivi_bytes){ 0 };
	if (o->dead || !o->hom[0]) { *err = "organism dead"; return false; }
	vivi_buf b = { 0 };
	if (!vivi_buf_append(&b, ORG_MAGIC, 4)
		|| !vivi_buf_push(&b, (uint8_t)((o->generation >> 8) & 255))
		|| !vivi_buf_push(&b, (uint8_t)(o->generation & 255))
		|| !vivi_buf_push(&b, (uint8_t)o->nchr)) {
		vivi_buf_free(&b);
		*err = "out of memory";
		return false;
	}
	for (size_t i = 0; i < o->nchr; i++) {
		const chr_opts *opt = &o->chr_opts[i];
		int mode = opt->codon ? 1 : 0;
		int h = (opt->h ? opt->h : 3) - 1;
		int gene_raw = opt->gene_raw ? opt->gene_raw : 1024;
		int units = opt->units ? opt->units : 4;
		uint8_t head[10] = {
			(uint8_t)o->chr_ids[i],
			(uint8_t)((gene_raw >> 8) & 255), (uint8_t)(gene_raw & 255),
			(uint8_t)((mode | (h << 1)) & 255),
			(uint8_t)units, (uint8_t)opt->flags,
			(uint8_t)((o->hlen[0][i] >> 24) & 255), (uint8_t)((o->hlen[0][i] >> 16) & 255),
			(uint8_t)((o->hlen[0][i] >> 8) & 255), (uint8_t)(o->hlen[0][i] & 255)
		};
		if (!vivi_buf_append(&b, head, 10)
			|| !vivi_buf_append(&b, o->hom[0][i], o->hlen[0][i])) {
			vivi_buf_free(&b);
			*err = "out of memory";
			return false;
		}
	}
	out->data = b.p;
	out->len = b.len;
	vivi_buf_release(&b);
	return true;
}

bool organism_deserialize(vivi_organism **out, const uint8_t *s, size_t len, const char **err)
{
	*out = nullptr;
	if (!s || len < 8) { *err = "bad format"; return false; }
	if (memcmp(s, ORG_MAGIC, 4) != 0) { *err = "bad magic"; return false; }
	int generation = s[4] * 256 + s[5];
	int nchr = s[6];
	if (nchr < 1) { *err = "bad chromosome count"; return false; }
	size_t pos = 7;
	int *ids = vivi_alloc((size_t)nchr * sizeof(int));
	chr_opts *opts_t = vivi_alloc((size_t)nchr * sizeof(chr_opts));
	uint8_t **g = vivi_alloc((size_t)nchr * sizeof(uint8_t *));
	size_t *lens = vivi_alloc((size_t)nchr * sizeof(size_t));
	if (!ids || !opts_t || !g || !lens) {
		vivi_dealloc(ids); vivi_dealloc(opts_t); vivi_dealloc(g); vivi_dealloc(lens);
		*err = "out of memory";
		return false;
	}
	int ok = 1;
	for (int i = 0; i < nchr && ok; i++) {
		if (pos + 10 > len) { *err = "truncated"; ok = 0; break; }
		int cid = s[pos];
		int gene_raw = s[pos + 1] * 256 + s[pos + 2];
		int mode = s[pos + 3] & 1;
		int h = ((s[pos + 3] >> 1) & 15) + 1;
		int units = s[pos + 4];
		int flags = s[pos + 5];
		size_t slen = ((size_t)s[pos + 6] << 24) | ((size_t)s[pos + 7] << 16)
			| ((size_t)s[pos + 8] << 8) | (size_t)s[pos + 9];
		pos += 10;
		if (pos + slen > len) { *err = "truncated"; ok = 0; break; }
		if (gene_raw < 16 || h < 3 || h > 12 || units < 1) {
			*err = "corrupt chromosome options";
			ok = 0;
			break;
		}
		chr_record rec;
		if (!chr_parse(&rec, s + pos, slen, -1, err) || !rec.cen_ok || rec.id != cid) {
			chr_record_free(&rec);
			*err = "corrupt chromosome";
			ok = 0;
			break;
		}
		chr_record_free(&rec);
		ids[i] = cid;
		opts_t[i].gene_raw = gene_raw;
		opts_t[i].codon = mode;
		opts_t[i].h = h;
		opts_t[i].units = units;
		opts_t[i].flags = flags;
		g[i] = vivi_alloc(slen ? slen : 1);
		if (!g[i]) { *err = "out of memory"; ok = 0; break; }
		memcpy(g[i], s + pos, slen);
		lens[i] = slen;
		pos += slen;
	}
	if (ok && !organism_alloc_from_genome(out, ids, opts_t, g, lens, (size_t)nchr, 60)) {
		*err = "out of memory";
		ok = 0;
	}
	for (int i = 0; i < nchr; i++) vivi_dealloc(g[i]);
	vivi_dealloc(g);
	vivi_dealloc(lens);
	vivi_dealloc(ids);
	vivi_dealloc(opts_t);
	if (ok) (*out)->generation = generation;
	return ok;
}

cell_report organism_maintain(vivi_organism **organisms, size_t count,
	const vivi_organism *stem)
{
	cell_report report;
	memset(&report, 0, sizeof(report));
	for (size_t i = 0; i < count; i++) {
		vivi_organism *o = organisms[i];
		if (o->dead || !o->hom[0]) {
			const char *err;
			if (stem && organism_renew(o, stem, &err)) report.renewed++;
			else {
				organism_kill(o);
				report.dead++;
			}
		} else {
			vivi_bytes *data = vivi_zalloc(o->nchr ? o->nchr : 1, sizeof(vivi_bytes));
			if (!data) continue;
			cell_report rep;
			const char *err;
			int ok = organism_read(&data, o, &rep, &err);
			report.repaired += rep.repaired + rep.structural;
			if (ok) {
				for (size_t k = 0; k < o->nchr; k++) vivi_bytes_free(&data[k]);
				if (o->generation >= o->max_gen && stem) {
					(void)organism_renew(o, stem, &err);
					report.renewed++;
				}
			} else if (stem && organism_renew(o, stem, &err)) {
				report.renewed++;
			} else {
				organism_kill(o);
				report.dead++;
			}
			vivi_dealloc(data);
		}
	}
	return report;
}






