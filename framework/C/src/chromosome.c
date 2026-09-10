/* chromosome.c -- chromosomes, transliteration of chromosome.lua */
#include "vivi/chromosome.h"
#include "vivi/parity.h"
#include <string.h>

/* absorbs messages when a caller passes err == nullptr */
static const char *chr_err_sink;

static constexpr uint8_t TELUNIT[3] = { 0xF2, 0xE3, 0x0E }; /* 2x TTAGGG */
static constexpr uint8_t CMARK[3] = { 0x6D, 0x6D, 0x6D };   /* CTGA x3 */
static constexpr uint8_t PMARK[3] = { 0x55, 0x55, 0x55 };   /* AAAA x3 */

static bool build_centromere(uint8_t out[CENBYTES], int chr_id, int flags, int ngenes,
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
	vivi_bytes packed = { 0 };
	if (!genome_pack_codons(&packed, v, 20, nullptr)) return false;
	memcpy(out, CMARK, 3);
	memcpy(out + 3, packed.data, packed.len);
	vivi_bytes_free(&packed);
	return true;
}

/* VIV14N centromere: raw 12 bytes (id, flags, ngenes, generation, parity,
 * reserved, rawlen) + 32-bit tag = 16 bytes -> 32 codons -> 24 bytes */
static bool build_centromere2(uint8_t out[CEN2BYTES], int chr_id, int flags, int ngenes,
	int generation, int parity, size_t rawlen)
{
	uint8_t raw[12] = { (uint8_t)chr_id, (uint8_t)flags,
		(uint8_t)((ngenes >> 8) & 255), (uint8_t)(ngenes & 255),
		(uint8_t)((generation >> 8) & 255), (uint8_t)(generation & 255),
		(uint8_t)parity, 0,
		(uint8_t)((rawlen >> 24) & 255), (uint8_t)((rawlen >> 16) & 255),
		(uint8_t)((rawlen >> 8) & 255), (uint8_t)(rawlen & 255) };
	uint32_t tag = genome_chaskey32(raw, 12);
	uint8_t block[16];
	memcpy(block, raw, 12);
	block[12] = (uint8_t)((tag >> 24) & 255);
	block[13] = (uint8_t)((tag >> 16) & 255);
	block[14] = (uint8_t)((tag >> 8) & 255);
	block[15] = (uint8_t)(tag & 255);
	int v[32];
	(void)genome_values_from_bytes(v, 32, block, 16);
	vivi_bytes packed = { 0 };
	if (!genome_pack_codons(&packed, v, 32, nullptr)) return false;
	memcpy(out, CMARK, 3);
	memcpy(out + 3, packed.data, packed.len);
	vivi_bytes_free(&packed);
	return true;
}

/* VIV14NB4NSH33 primer site: raw 4 bytes (barcode BE, zero, zero) plus a
 * 32-bit tag = 8 bytes -> 16 codons -> 12 bytes after the marker */
static bool build_primer(uint8_t out[PRIMERBYTES], int barcode)
{
	uint8_t raw[4] = { (uint8_t)((barcode >> 8) & 255), (uint8_t)(barcode & 255), 0, 0 };
	uint32_t tag = genome_chaskey32(raw, 4);
	uint8_t block[8] = { raw[0], raw[1], raw[2], raw[3],
		(uint8_t)((tag >> 24) & 255), (uint8_t)((tag >> 16) & 255),
		(uint8_t)((tag >> 8) & 255), (uint8_t)(tag & 255) };
	int v[16];
	(void)genome_values_from_bytes(v, 16, block, 8);
	vivi_bytes packed = { 0 };
	if (!genome_pack_codons(&packed, v, 16, nullptr)) return false;
	memcpy(out, PMARK, 3);
	memcpy(out + 3, packed.data, packed.len);
	vivi_bytes_free(&packed);
	return true;
}

/* parses a primer site at p; 0 is not a valid barcode */
static bool parse_primer(const uint8_t *strand, size_t slen, size_t p, int *barcode)
{
	if (p + PRIMERBYTES > slen || memcmp(strand + p, PMARK, 3) != 0) return false;
	int v[16];
	if (!genome_unpack_codons(v, 16, strand + p + 3, slen - p - 3, nullptr)) return false;
	uint8_t raw[4];
	for (int i = 0; i < 4; i++) raw[i] = (uint8_t)(v[i * 2] * 16 + v[i * 2 + 1]);
	uint32_t expected = ((uint32_t)(v[8] * 16 + v[9]) << 24)
		| ((uint32_t)(v[10] * 16 + v[11]) << 16)
		| ((uint32_t)(v[12] * 16 + v[13]) << 8)
		| (uint32_t)(v[14] * 16 + v[15]);
	if (genome_chaskey32(raw, 4) != expected) return false;
	*barcode = raw[0] * 256 + raw[1];
	return true;
}

bool chr_encode(vivi_bytes *out, int chr_id, const uint8_t *data, size_t len,
	const chr_opts *opts, const char **err)
{
	if (!err) err = &chr_err_sink;
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
	int parity = opts->parity;
	if (parity < 0 || parity > 16) { *err = "invalid parity (integer 0..16)"; return false; }
	int primer = opts->primer;
	if (primer < 0 || primer > 65535) { *err = "invalid primer (integer 0..65535)"; return false; }

	vivi_buf o = { 0 };
	size_t ngenes = 0, npar = 0;
	vivi_bytes *gene_strands = nullptr;
	vivi_bytes par_strands[16];
	size_t gene_cap = 0;
	size_t tb = (size_t)units * 3;
	if (!vivi_buf_reserve(&o, tb * 2 + 2 * PRIMERBYTES + CEN2BYTES + len + len / 16 + 64)) {
		*err = "out of memory";
		return false;
	}
	for (int u = 0; u < units; u++)
		if (!vivi_buf_append(&o, TELUNIT, 3)) { *err = "out of memory"; goto fail; }
	/* VIV14NB4NSH33: forward primer site right after the left telomere */
	if (primer > 0) {
		uint8_t pl[PRIMERBYTES];
		if (!build_primer(pl, primer) || !vivi_buf_append(&o, pl, PRIMERBYTES)) {
			*err = "out of memory";
			goto fail;
		}
	}

	/* pack data genes */
	for (size_t off = 0; off < len; off += (size_t)gene_raw) {
		if (ngenes >= 256) {   /* gene ids are one byte: 0..255 */
			*err = "too many genes (max 256)";
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
	if (parity > 0) {
		if (ngenes == 0) { *err = "parity needs a non-empty payload"; goto fail; }
		if (ngenes + (size_t)parity > 255) {
			*err = "too many genes with parity (max 255)";
			goto fail;
		}
	}

	/* centromere (VIV14N revision when parity > 0) */
	if (parity == 0) {
		uint8_t cen[CENBYTES];
		if (!build_centromere(cen, chr_id, flags, (int)ngenes, 0)
			|| !vivi_buf_append(&o, cen, CENBYTES)) {
			*err = "out of memory";
			goto fail;
		}
	} else {
		/* zero-padded data shards -> RS parity genes */
		uint8_t *padded = vivi_zalloc(ngenes * (size_t)gene_raw, 1);
		if (!padded) { *err = "out of memory"; goto fail; }
		for (size_t i = 0; i < ngenes; i++) {
			size_t off = i * (size_t)gene_raw;
			size_t clen = (len - off < (size_t)gene_raw) ? (len - off) : (size_t)gene_raw;
			memcpy(padded + off, data + off, clen);
		}
		const uint8_t *shard_in[255];
		for (size_t i = 0; i < ngenes; i++) shard_in[i] = padded + i * (size_t)gene_raw;
		uint8_t *shards[255];
		if (!vivi_parity_encode(shards, shard_in, ngenes, (size_t)parity,
			(size_t)gene_raw, err)) {
			vivi_dealloc(padded);
			goto fail;
		}
		uint8_t cen[CEN2BYTES];
		if (!build_centromere2(cen, chr_id, flags, (int)(ngenes + (size_t)parity), 0,
			parity, len)) {
			*err = "out of memory";
			vivi_parity_release(shards, ngenes + (size_t)parity);
			vivi_dealloc(padded);
			goto fail;
		}
		int ok = vivi_buf_append(&o, cen, CEN2BYTES);
		for (size_t j = 0; ok && j < (size_t)parity; j++) {
			if (!genome_gene_encode(&par_strands[npar], (int)(ngenes + j), 2, opts->codon,
				opts->h ? opts->h : 3, shards[ngenes + j], (size_t)gene_raw, err)) {
				ok = 0;
				break;
			}
			npar++;
		}
		vivi_parity_release(shards, ngenes + (size_t)parity);
		vivi_dealloc(padded);
		if (!ok) { *err = "out of memory"; goto fail; }
	}
	for (size_t i = 0; i < ngenes; i++)
		if (!vivi_buf_append(&o, gene_strands[i].data, gene_strands[i].len)) {
			*err = "out of memory";
			goto fail;
		}
	for (size_t j = 0; j < npar; j++)
		if (!vivi_buf_append(&o, par_strands[j].data, par_strands[j].len)) {
			*err = "out of memory";
			goto fail;
		}
	/* VIV14NB4NSH33: reverse primer site right before the right telomere */
	if (primer > 0) {
		uint8_t pr[PRIMERBYTES];
		if (!build_primer(pr, primer) || !vivi_buf_append(&o, pr, PRIMERBYTES)) {
			*err = "out of memory";
			goto fail;
		}
	}
	for (int u = 0; u < units; u++)
		if (!vivi_buf_append(&o, TELUNIT, 3)) { *err = "out of memory"; goto fail; }
	for (size_t i = 0; i < ngenes; i++) vivi_bytes_free(&gene_strands[i]);
	for (size_t j = 0; j < npar; j++) vivi_bytes_free(&par_strands[j]);
	vivi_dealloc(gene_strands);
	out->data = o.p;
	out->len = o.len;
	vivi_buf_release(&o);
	return true;
fail:
	for (size_t i = 0; i < ngenes; i++) vivi_bytes_free(&gene_strands[i]);
	for (size_t j = 0; j < npar; j++) vivi_bytes_free(&par_strands[j]);
	vivi_dealloc(gene_strands);
	vivi_buf_free(&o);
	return false;
}

typedef struct {
	int id, flags, ngenes, generation, parity;
	size_t rawlen;
} cen_fields;

/* parses either centromere revision at p; 2 = VIV14N, 1 = VIV1/VIV14, 0 = no */
static int cen_parse(const uint8_t *strand, size_t slen, size_t p, cen_fields *f)
{
	if (p + CENBYTES > slen || memcmp(strand + p, CMARK, 3) != 0) return 0;
	int v[32];
	if (p + CEN2BYTES <= slen
		&& genome_unpack_codons(v, 32, strand + p + 3, slen - p - 3, nullptr)) {
		uint8_t raw[12];
		for (int i = 0; i < 12; i++) raw[i] = (uint8_t)(v[i * 2] * 16 + v[i * 2 + 1]);
		uint32_t expected = ((uint32_t)(v[24] * 16 + v[25]) << 24)
			| ((uint32_t)(v[26] * 16 + v[27]) << 16)
			| ((uint32_t)(v[28] * 16 + v[29]) << 8)
			| (uint32_t)(v[30] * 16 + v[31]);
		if (genome_chaskey32(raw, 12) == expected) {
			f->id = raw[0];
			f->flags = raw[1];
			f->ngenes = raw[2] * 256 + raw[3];
			f->generation = raw[4] * 256 + raw[5];
			f->parity = raw[6];
			f->rawlen = ((size_t)raw[8] << 24) | ((size_t)raw[9] << 16)
				| ((size_t)raw[10] << 8) | (size_t)raw[11];
			return 2;
		}
	}
	if (genome_unpack_codons(v, 20, strand + p + 3, slen - p - 3, nullptr)) {
		uint8_t raw[6];
		for (int i = 0; i < 6; i++) raw[i] = (uint8_t)(v[i * 2] * 16 + v[i * 2 + 1]);
		uint32_t expected = ((uint32_t)(v[12] * 16 + v[13]) << 24)
			| ((uint32_t)(v[14] * 16 + v[15]) << 16)
			| ((uint32_t)(v[16] * 16 + v[17]) << 8)
			| (uint32_t)(v[18] * 16 + v[19]);
		if (genome_chaskey32(raw, 6) == expected) {
			f->id = raw[0];
			f->flags = raw[1];
			f->ngenes = raw[2] * 256 + raw[3];
			f->generation = raw[4] * 256 + raw[5];
			f->parity = 0;
			f->rawlen = 0;
			return 1;
		}
	}
	return 0;
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
	if (!err) err = &chr_err_sink;
	memset(out, 0, sizeof(*out));
	if (!strand) { *err = "expected string"; return false; }
	/* auto-detect the telomere size from the centromere marker position;
	 * the candidate is validated by the centromere tag at that offset */
	int units = (units_hint > 0) ? units_hint : 4;
	size_t pc = find_bytes(strand, slen, CMARK, 3, 0);
	cen_fields fdet;
	memset(&fdet, 0, sizeof(fdet));
	if (pc > 0 && cen_parse(strand, slen, pc, &fdet)) {
		/* a valid primer site shifts the centromere by PRIMERBYTES */
		size_t cstart = pc;
		int bc = 0;
		if (pc >= PRIMERBYTES && parse_primer(strand, slen, pc - PRIMERBYTES, &bc))
			cstart = pc - PRIMERBYTES;
		if (cstart % 3 == 0) units = (int)(cstart / 3);
	}
	size_t tb = (size_t)units * 3;
	if (slen < 2 * tb + CENBYTES) { *err = "chromosome too short"; return false; }
	int telo_ok = (memcmp(strand, TELUNIT, 3) == 0);
	for (size_t u = 1; telo_ok && u < (size_t)units; u++)
		if (memcmp(strand + u * 3, TELUNIT, 3) != 0) telo_ok = 0;
	for (size_t u = 0; telo_ok && u < (size_t)units; u++)
		if (memcmp(strand + slen - tb + u * 3, TELUNIT, 3) != 0) telo_ok = 0;
	int cen_ok = 0;
	int id = 0, flags = 0, ngenes = 0, generation = 0, parity = 0, cen_version = 0;
	size_t rawlen = 0;
	int primer = 0, primer_ok = 0, primer_bytes = 0;
	/* VIV14NB4NSH33: a valid left primer site shifts the centromere */
	size_t cenp = tb;
	if (tb + PRIMERBYTES <= slen && parse_primer(strand, slen, tb, &primer)) {
		primer_bytes = PRIMERBYTES;
		cenp = tb + PRIMERBYTES;
	}
	cen_fields f;
	memset(&f, 0, sizeof(f));
	int cv = cen_parse(strand, slen, cenp, &f);
	if (cv) {
		if (cv == 2 && slen < 2 * tb + CEN2BYTES) {
			*err = "chromosome too short";
			return false;
		}
		cen_ok = 1;
		cen_version = cv;
		id = f.id;
		flags = f.flags;
		ngenes = f.ngenes;
		generation = f.generation;
		parity = f.parity;
		rawlen = f.rawlen;
	}
	/* the reverse site must exist and carry the same barcode */
	if (primer > 0 && cenp + tb + PRIMERBYTES <= slen) {
		int rp = 0;
		if (parse_primer(strand, slen, slen - tb - PRIMERBYTES, &rp) && rp == primer)
			primer_ok = 1;
	}
	genome_gene *genes = nullptr;
	size_t gene_count = 0;
	if (cen_ok) {
		size_t cen_bytes = (cen_version == 2) ? CEN2BYTES : CENBYTES;
		size_t interior_off = cenp + cen_bytes;
		size_t interior_len = slen - interior_off - tb - (size_t)primer_bytes;
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
	out->parity = parity;
	out->rawlen = rawlen;
	out->cen_version = cen_version;
	out->primer = primer;
	out->primer_ok = primer_ok;
	out->primer_bytes = primer_bytes;
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
	if (!err) err = &chr_err_sink;
	*out = (vivi_bytes){ 0 };
	if (!rec || !rec->cen_ok) { *err = "centromere damaged"; return false; }

	if (rec->parity == 0) {
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

	/* VIV14N: data shards + parity genes over zero-padded gene shards */
	size_t total_genes = (size_t)rec->ngenes;
	size_t m = (size_t)rec->parity;
	size_t rawlen = rec->rawlen;
	if (total_genes > 255 || m > total_genes || rawlen == 0) {
		*err = "invalid parity layout";
		return false;
	}
	size_t n = total_genes - m;
	int seen[256];
	memset(seen, -1, sizeof(seen));
	for (size_t i = 0; i < rec->gene_count; i++) {
		const genome_gene *g = &rec->genes[i];
		if (g->id < 0 || (size_t)g->id >= total_genes) {
			*err = "gene id out of range";
			return false;
		}
		if (seen[g->id] >= 0) { *err = "duplicate gene id"; return false; }
		seen[g->id] = (int)i;
	}
	int all_data = 1;
	for (size_t k = 0; k < n; k++)
		if (seen[k] < 0 || !rec->genes[seen[k]].crc_ok) { all_data = 0; break; }
	uint8_t *res = vivi_alloc(rawlen);
	if (!res) { *err = "out of memory"; return false; }
	size_t pos = 0;
	if (all_data) {
		for (size_t k = 0; k < n && pos < rawlen; k++) {
			const genome_gene *g = &rec->genes[seen[k]];
			size_t take = (size_t)g->data_len;
			if (take > rawlen - pos) take = rawlen - pos;
			memcpy(res + pos, g->data, take);
			pos += take;
		}
		if (pos != rawlen) { vivi_dealloc(res); *err = "payload length mismatch"; return false; }
		out->data = res;
		out->len = rawlen;
		return true;
	}
	vivi_dealloc(res);

	/* an erasure must be covered by parity: the shard length comes from a
	 * readable parity gene */
	size_t shard_len = 0;
	for (size_t k = n; k < total_genes; k++) {
		int idx = seen[k];
		if (idx >= 0 && rec->genes[idx].crc_ok) {
			shard_len = (size_t)rec->genes[idx].rawlen;
			break;
		}
	}
	if (shard_len == 0 || rawlen > n * shard_len) { *err = "too few genes"; return false; }
	uint8_t *shardbuf[256];
	const uint8_t *refs[256];
	uint8_t present[256];
	memset(shardbuf, 0, sizeof(shardbuf));
	memset(refs, 0, sizeof(refs));
	memset(present, 0, sizeof(present));
	for (size_t k = 0; k < total_genes; k++) {
		int idx = seen[k];
		if (idx < 0 || !rec->genes[idx].crc_ok) continue;
		const genome_gene *g = &rec->genes[idx];
		shardbuf[k] = vivi_zalloc(shard_len, 1);
		if (!shardbuf[k]) goto oom;
		size_t take = (size_t)g->data_len < shard_len ? (size_t)g->data_len : shard_len;
		memcpy(shardbuf[k], g->data, take);
		refs[k] = shardbuf[k];
		present[k] = 1;
	}
	uint8_t *decoded[255];
	if (!vivi_parity_decode(decoded, refs, present, n, m, shard_len, err)) {
		for (size_t k = 0; k < total_genes; k++) vivi_dealloc(shardbuf[k]);
		return false;
	}
	res = vivi_alloc(rawlen);
	if (!res) {
		vivi_parity_release(decoded, n);
		for (size_t k = 0; k < total_genes; k++) vivi_dealloc(shardbuf[k]);
		*err = "out of memory";
		return false;
	}
	pos = 0;
	for (size_t k = 0; k < n && pos < rawlen; k++) {
		size_t take = shard_len;
		if (take > rawlen - pos) take = rawlen - pos;
		memcpy(res + pos, decoded[k], take);
		pos += take;
	}
	vivi_parity_release(decoded, n);
	for (size_t k = 0; k < total_genes; k++) vivi_dealloc(shardbuf[k]);
	if (pos != rawlen) { vivi_dealloc(res); *err = "payload length mismatch"; return false; }
	out->data = res;
	out->len = rawlen;
	return true;
oom:
	for (size_t k = 0; k < total_genes; k++) vivi_dealloc(shardbuf[k]);
	*err = "out of memory";
	return false;
}

bool chr_set_generation(vivi_bytes *out, const uint8_t *strand, size_t slen, int generation,
	const char **err)
{
	if (!err) err = &chr_err_sink;
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
	size_t pb = (size_t)rec.primer_bytes;
	uint8_t cen[CEN2BYTES];
	size_t cen_bytes;
	if (rec.cen_version == 2) {
		if (!build_centromere2(cen, rec.id, rec.flags, rec.ngenes, generation,
			rec.parity, rec.rawlen)) {
			chr_record_free(&rec);
			*err = "out of memory";
			return false;
		}
		cen_bytes = CEN2BYTES;
	} else {
		if (!build_centromere(cen, rec.id, rec.flags, rec.ngenes, generation)) {
			chr_record_free(&rec);
			*err = "out of memory";
			return false;
		}
		cen_bytes = CENBYTES;
	}
	uint8_t *ns = vivi_alloc(slen);
	if (!ns) {
		chr_record_free(&rec);
		*err = "out of memory";
		return false;
	}
	memcpy(ns, strand, tb + pb);
	memcpy(ns + tb + pb, cen, cen_bytes);
	memcpy(ns + tb + pb + cen_bytes, strand + tb + pb + cen_bytes,
		slen - tb - pb - cen_bytes);
	chr_record_free(&rec);
	out->data = ns;
	out->len = slen;
	return true;
}






