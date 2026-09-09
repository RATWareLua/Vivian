/* chromosome.c -- chromosomes, transliteration of chromosome.lua */
#include "vivi/chromosome.h"
#include <string.h>

static constexpr uint8_t TELUNIT[3] = { 0xF2, 0xE3, 0x0E }; /* 2x TTAGGG */
static constexpr uint8_t CMARK[3] = { 0x6D, 0x6D, 0x6D };   /* CTGA x3 */
#define CENBYTES 18   /* marker (3) + 20 codons (15) */

static void build_centromere(uint8_t out[CENBYTES], int chr_id, int flags, int ngenes,
	int generation)
{
	uint8_t raw[6] = { (uint8_t)chr_id, (uint8_t)flags,
		(uint8_t)((ngenes >> 8) & 255), (uint8_t)(ngenes & 255),
		(uint8_t)((generation >> 8) & 255), (uint8_t)(generation & 255) };
	uint32_t tag = genome_chaskey32(raw, 6);
	uint8_t block[10];
	memcpy(block, raw, 6);
	block[6] = (uint8_t)((tag >> 24) & 255);
	block[7] = (uint8_t)((tag >> 16) & 255);
	block[8] = (uint8_t)((tag >> 8) & 255);
	block[9] = (uint8_t)(tag & 255);
	int v[20];
	(void)genome_values_from_bytes(v, 20, block, 10);
	vivi_bytes packed;
	(void)genome_pack_codons(&packed, v, 20, nullptr);
	memcpy(out, CMARK, 3);
	memcpy(out + 3, packed.data, packed.len);
	vivi_bytes_free(&packed);
}

bool chr_encode(vivi_bytes *out, int chr_id, const uint8_t *data, size_t len,
	const chr_opts *opts, const char **err)
{
	*out = (vivi_bytes){ 0 };
	chr_opts def;
	memset(&def, 0, sizeof(def));
	if (!opts) opts = &def;
	int flags = opts->flags;
	if (flags < 0 || flags > 255) { *err = "invalid flags (byte 0..255)"; return false; }
	int units = opts->units ? opts->units : 4;
	if (units < 1) { *err = "invalid units (integer >= 1)"; return false; }
	int gene_raw = opts->gene_raw ? opts->gene_raw : 1024;
	if (gene_raw < 16) { *err = "invalid gene_raw (integer >= 16)"; return false; }
	if (chr_id < 0 || chr_id > 255) { *err = "invalid chr_id (byte 0..255)"; return false; }

	vivi_buf o = { 0 };
	size_t tb = (size_t)units * 3;
	if (!vivi_buf_reserve(&o, tb * 2 + CENBYTES + len + len / 16 + 64)) {
		*err = "out of memory";
		return false;
	}
	for (int u = 0; u < units; u++)
		if (!vivi_buf_append(&o, TELUNIT, 3)) { *err = "out of memory"; return false; }
	uint8_t cen[CENBYTES];
	build_centromere(cen, chr_id, flags, 0, 0); /* ngenes patched after packing genes */
	size_t cen_pos = o.len;
	if (!vivi_buf_append(&o, cen, CENBYTES)) { *err = "out of memory"; return false; }

	/* pack genes */
	size_t ngenes = 0;
	vivi_bytes *gene_strands = nullptr;
	size_t gene_cap = 0;
	for (size_t off = 0; off < len; off += (size_t)gene_raw) {
		if (ngenes > 65535) {
			*err = "too many genes (max 65535)";
			goto fail;
		}
		size_t clen = (len - off < (size_t)gene_raw) ? (len - off) : (size_t)gene_raw;
		vivi_bytes g;
		if (!genome_gene_encode(&g, (int)ngenes, 0, opts->codon, opts->h ? opts->h : 3,
			data + off, clen, err)) {
			goto fail;
		}
		if (ngenes == gene_cap) {
			size_t ncap = gene_cap ? gene_cap * 2 : 8;
			vivi_bytes *ng = vivi_alloc(ncap * sizeof(vivi_bytes));
			if (!ng) {
				vivi_bytes_free(&g);
				*err = "out of memory";
				goto fail;
			}
			if (gene_cap) memcpy(ng, gene_strands, gene_cap * sizeof(vivi_bytes));
			vivi_dealloc(gene_strands);
			gene_strands = ng;
			gene_cap = ncap;
		}
		gene_strands[ngenes++] = g;
	}
	/* rebuild the centromere with the real gene count and append genes */
	build_centromere(cen, chr_id, flags, (int)ngenes, 0);
	for (size_t i = 0; i < (size_t)CENBYTES; i++) o.p[cen_pos + i] = cen[i];
	for (size_t i = 0; i < ngenes; i++)
		if (!vivi_buf_append(&o, gene_strands[i].data, gene_strands[i].len)) {
			*err = "out of memory";
			goto fail;
		}
	for (int u = 0; u < units; u++)
		if (!vivi_buf_append(&o, TELUNIT, 3)) { *err = "out of memory"; return false; }
	for (size_t i = 0; i < ngenes; i++) vivi_bytes_free(&gene_strands[i]);
	vivi_dealloc(gene_strands);
	out->data = o.p;
	out->len = o.len;
	vivi_buf_release(&o);
	return true;
fail:
	for (size_t i = 0; i < ngenes; i++) vivi_bytes_free(&gene_strands[i]);
	vivi_dealloc(gene_strands);
	vivi_buf_free(&o);
	return false;
}

static int cen_valid_at(const uint8_t *strand, size_t slen, size_t p)
{
	if (p + CENBYTES > slen) return 0;
	if (memcmp(strand + p, CMARK, 3) != 0) return 0;
	int v[20];
	(void)genome_unpack_codons(v, 20, strand + p + 3, slen - p - 3, nullptr);
	uint8_t raw[6];
	for (int i = 0; i < 6; i++) raw[i] = (uint8_t)(v[i * 2] * 16 + v[i * 2 + 1]);
	uint32_t expected = ((uint32_t)(v[12] * 16 + v[13]) << 24)
		| ((uint32_t)(v[14] * 16 + v[15]) << 16)
		| ((uint32_t)(v[16] * 16 + v[17]) << 8)
		| (uint32_t)(v[18] * 16 + v[19]);
	return genome_chaskey32(raw, 6) == expected;
}

static size_t find_bytes(const uint8_t *hay, size_t hlen, const uint8_t *needle, size_t nlen,
	size_t from)
{
	if (nlen == 0 || hlen < nlen) return hlen;
	for (size_t i = from; i + nlen <= hlen; i++)
		if (memcmp(hay + i, needle, nlen) == 0) return i;
	return hlen;
}

bool chr_parse(chr_record *out, const uint8_t *strand, size_t slen, int units_hint,
	const char **err)
{
	memset(out, 0, sizeof(*out));
	if (!strand) { *err = "expected string"; return false; }
	/* auto-detect the telomere size from the centromere marker position;
	 * the candidate is validated by the centromere tag at that offset */
	int units = (units_hint > 0) ? units_hint : 4;
	size_t pc = find_bytes(strand, slen, CMARK, 3, 0);
	if (pc > 0 && pc % 3 == 0 && cen_valid_at(strand, slen, pc))
		units = (int)(pc / 3);
	size_t tb = (size_t)units * 3;
	if (slen < 2 * tb + CENBYTES) { *err = "chromosome too short"; return false; }
	int telo_ok = (memcmp(strand, TELUNIT, 3) == 0);
	for (size_t u = 1; telo_ok && u < (size_t)units; u++)
		if (memcmp(strand + u * 3, TELUNIT, 3) != 0) telo_ok = 0;
	for (size_t u = 0; telo_ok && u < (size_t)units; u++)
		if (memcmp(strand + slen - tb + u * 3, TELUNIT, 3) != 0) telo_ok = 0;
	int cen_ok = 0;
	int id = 0, flags = 0, ngenes = 0, generation = 0;
	size_t cenp = tb;
	if (cenp + CENBYTES <= slen) {
		int v[20];
		if (memcmp(strand + cenp, CMARK, 3) == 0) {
			(void)genome_unpack_codons(v, 20, strand + cenp + 3, slen - cenp - 3, nullptr);
			uint8_t raw[6];
			for (int i = 0; i < 6; i++) raw[i] = (uint8_t)(v[i * 2] * 16 + v[i * 2 + 1]);
			uint32_t expected = ((uint32_t)(v[12] * 16 + v[13]) << 24)
				| ((uint32_t)(v[14] * 16 + v[15]) << 16)
				| ((uint32_t)(v[16] * 16 + v[17]) << 8)
				| (uint32_t)(v[18] * 16 + v[19]);
			if (genome_chaskey32(raw, 6) == expected) {
				cen_ok = 1;
				id = raw[0];
				flags = raw[1];
				ngenes = raw[2] * 256 + raw[3];
				generation = raw[4] * 256 + raw[5];
			}
		}
	}
	genome_gene *genes = nullptr;
	size_t gene_count = 0;
	if (cen_ok) {
		size_t interior_off = tb + CENBYTES;
		size_t interior_len = slen - interior_off - tb;
		genome_scan_result res;
		if (!genome_gene_scan(&res, strand + interior_off, interior_len, err)) {
			chr_record_free(out);
			return false;
		}
		gene_count = res.count;
		genes = res.genes;
		/* convert to absolute offsets */
		for (size_t i = 0; i < gene_count; i++)
			genes[i].offset += interior_off;
	}
	out->id = id;
	out->flags = flags;
	out->ngenes = ngenes;
	out->generation = generation;
	out->telomere_ok = telo_ok;
	out->cen_ok = cen_ok;
	out->telomere_bytes = tb;
	out->genes = genes;
	out->gene_count = gene_count;
	return true;
}

void chr_record_free(chr_record *r)
{
	if (!r) return;
	if (r->genes) {
		for (size_t i = 0; i < r->gene_count; i++) vivi_dealloc(r->genes[i].data);
		vivi_dealloc(r->genes);
	}
	r->genes = nullptr;
	r->gene_count = 0;
}

bool chr_read(vivi_bytes *out, const chr_record *rec, const char **err)
{
	*out = (vivi_bytes){ 0 };
	if (!rec || !rec->cen_ok) { *err = "centromere damaged"; return false; }
	if ((int)rec->gene_count != rec->ngenes) { *err = "gene count mismatch"; return false; }
	int seen[256];
	memset(seen, -1, sizeof(seen));
	for (size_t i = 0; i < rec->gene_count; i++) {
		const genome_gene *g = &rec->genes[i];
		if (seen[g->id] >= 0) { *err = "duplicate gene id"; return false; }
		seen[g->id] = (int)i;
	}
	for (int gid = 0; gid < rec->ngenes; gid++) {
		int idx = seen[gid];
		if (idx < 0) { *err = "gene missing"; return false; }
		const genome_gene *g = &rec->genes[idx];
		if (!g->crc_ok) { *err = "gene damaged"; return false; }
	}
	size_t total = 0;
	for (int gid = 0; gid < rec->ngenes; gid++)
		total += rec->genes[seen[gid]].data_len;
	out->data = vivi_alloc(total ? total : 1);
	if (!out->data) { *err = "out of memory"; return false; }
	size_t pos = 0;
	for (int gid = 0; gid < rec->ngenes; gid++) {
		const genome_gene *g = &rec->genes[seen[gid]];
		memcpy(out->data + pos, g->data, g->data_len);
		pos += g->data_len;
	}
	out->len = total;
	return true;
}

bool chr_set_generation(vivi_bytes *out, const uint8_t *strand, size_t slen, int generation,
	const char **err)
{
	*out = (vivi_bytes){ 0 };
	chr_record rec;
	if (!chr_parse(&rec, strand, slen, -1, err)) return false;
	if (!rec.cen_ok) { *err = "centromere damaged"; return false; }
	if (generation < 0 || generation > 65535) {
		chr_record_free(&rec);
		*err = "invalid generation (0..65535)";
		return false;
	}
	size_t tb = rec.telomere_bytes;
	uint8_t cen[CENBYTES];
	build_centromere(cen, rec.id, rec.flags, rec.ngenes, generation);
	uint8_t *ns = vivi_alloc(slen);
	if (!ns) {
		chr_record_free(&rec);
		*err = "out of memory";
		return false;
	}
	memcpy(ns, strand, tb);
	memcpy(ns + tb, cen, CENBYTES);
	memcpy(ns + tb + CENBYTES, strand + tb + CENBYTES, slen - tb - CENBYTES);
	chr_record_free(&rec);
	out->data = ns;
	out->len = slen;
	return true;
}






