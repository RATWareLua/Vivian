/* genome.c -- genes, codon code and Chaskey-12, transliteration of genome.lua */
#include "vivi/genome.h"
#include <string.h>

static constexpr uint8_t genome_PROM[3] = { 0x1B, 0x1B, 0x1B };   /* ACGT x3 */
static constexpr uint8_t genome_TERM[3] = { 0xE4, 0xE4, 0xE4 };   /* TGAC x3 */

static uint32_t rotl32(uint32_t x, int n)
{
	return (x << n) | (x >> (32 - n));
}

/* ---- Chaskey-12 (ISO/IEC 29192-6), official reference translation ---- */

static constexpr uint32_t CHASKEY_KEY[4] = { 0x1B87'3593u, 0x9E37'79B9u, 0x85EB'CA6Bu, 0xC2B2'AE35u };
static constexpr uint32_t C87[2] = { 0x00u, 0x87u };

static void times_two(uint32_t out[4], const uint32_t in[4])
{
	uint32_t c = in[3] >> 31;
	out[0] = (in[0] << 1) ^ C87[c];
	out[1] = (in[1] << 1) | (in[0] >> 31);
	out[2] = (in[2] << 1) | (in[1] >> 31);
	out[3] = (in[3] << 1) | (in[2] >> 31);
}

static void permute12(uint32_t v[4])
{
	for (int i = 0; i < 12; i++) {
		v[0] += v[1]; v[1] = rotl32(v[1], 5); v[1] ^= v[0]; v[0] = rotl32(v[0], 16);
		v[2] += v[3]; v[3] = rotl32(v[3], 8); v[3] ^= v[2];
		v[0] += v[3]; v[3] = rotl32(v[3], 13); v[3] ^= v[0];
		v[2] += v[1]; v[1] = rotl32(v[1], 7); v[1] ^= v[2]; v[2] = rotl32(v[2], 16);
	}
}

static uint32_t load32(const uint8_t *s, size_t off)
{
	return (uint32_t)s[off] | ((uint32_t)s[off + 1] << 8)
		| ((uint32_t)s[off + 2] << 16) | ((uint32_t)s[off + 3] << 24);
}

void genome_chaskey_tag(uint32_t out[4], const uint8_t *s, size_t len, const uint32_t key[4])
{
	const uint32_t *k = key ? key : CHASKEY_KEY;
	uint32_t k1[4], k2[4];
	times_two(k1, k);
	times_two(k2, k1);
	size_t len_rem = len % 16;
	uint32_t v[4] = { k[0], k[1], k[2], k[3] };
	size_t limit;
	if (len_rem == 0) limit = (len >= 16) ? len - 16 : 0;
	else limit = len - len_rem;
	for (size_t off = 0; off + 16 <= limit; off += 16) {
		v[0] ^= load32(s, off);
		v[1] ^= load32(s, off + 4);
		v[2] ^= load32(s, off + 8);
		v[3] ^= load32(s, off + 12);
		permute12(v);
	}
	uint32_t m[4];
	const uint32_t *l;
	if (len_rem == 0 && len > 0) {
		size_t off = len - 16;
		m[0] = load32(s, off);
		m[1] = load32(s, off + 4);
		m[2] = load32(s, off + 8);
		m[3] = load32(s, off + 12);
		l = k1;
	} else {
		uint8_t p[16] = { 0 };
		size_t base = len - len_rem;
		for (size_t i = 0; i < len_rem; i++) p[i] = s[base + i];
		p[len_rem] = 0x01; /* padding bit (10*) */
		m[0] = load32(p, 0);
		m[1] = load32(p, 4);
		m[2] = load32(p, 8);
		m[3] = load32(p, 12);
		l = k2;
	}
	for (int i = 0; i < 4; i++) v[i] ^= m[i] ^ l[i];
	permute12(v);
	for (int i = 0; i < 4; i++) out[i] = v[i] ^ l[i];
}

uint32_t genome_chaskey32(const uint8_t *s, size_t len)
{
	uint32_t t[4];
	genome_chaskey_tag(t, s, len, nullptr);
	return t[0];
}

/* ---- codon packing ---- */

bool genome_pack_codons(vivi_bytes *out, const int *values, size_t nvalues, const char **err)
{
	*out = (vivi_bytes){ 0 };
	if (nvalues % 4 != 0) { *err = "codon count must be a multiple of 4"; return false; }
	vivi_buf o = { 0 };
	if (!vivi_buf_reserve(&o, nvalues / 4 + 4)) { *err = "out of memory"; return false; }
	int acc = 0, nib = 0;
	for (size_t i = 0; i < nvalues; i++) {
		int c = values[i];
		int b1 = (c >> 2) & 3;
		int b2 = c & 3;
		int w = b2 ^ 1; /* wobble != b2, so every codon ends run 1 */
		const int bases[3] = { b1, b2, w };
		for (int j = 0; j < 3; j++) {
			acc = acc * 4 + bases[j];
			nib++;
			if (nib == 4) {
				if (!vivi_buf_push(&o, (uint8_t)acc)) {
					vivi_buf_free(&o);
					*err = "out of memory";
					return false;
				}
				acc = 0;
				nib = 0;
			}
		}
	}
	if (nib != 0) {
		vivi_buf_free(&o);
		*err = "codon count must be a multiple of 4";
		return false;
	}
	out->data = o.p;
	out->len = o.len;
	vivi_buf_release(&o);
	return true;
}

bool genome_unpack_codons(int *values, size_t nvalues, const uint8_t *packed, size_t plen,
	const char **err)
{
	if (nvalues * 6 > plen * 8) { *err = "values do not fit the packed block"; return false; }
	(void)err;
	for (size_t kk = 0; kk < nvalues; kk++) {
		size_t bitpos = kk * 6;
		size_t off1 = bitpos % 8;
		size_t off2 = (bitpos + 2) % 8;
		int b1 = (int)((packed[bitpos / 8] >> (6 - off1)) & 3u);
		int b2 = (int)((packed[(bitpos + 2) / 8] >> (6 - off2)) & 3u);
		values[kk] = b1 * 4 + b2;
	}
	return true;
}

bool genome_values_from_bytes(int *values, size_t nvalues, const uint8_t *s, size_t len)
{
	if (nvalues < len * 2) return false;
	for (size_t i = 0; i < len; i++) {
		values[2 * i] = (int)((s[i] >> 4) & 15u);
		values[2 * i + 1] = (int)(s[i] & 15u);
	}
	return true;
}

bool genome_bytes_from_values(vivi_bytes *out, const int *values, size_t nvalues, size_t nbytes,
	const char **err)
{
	*out = (vivi_bytes){ 0 };
	if (nvalues < nbytes * 2) { *err = "not enough values"; return false; }
	out->data = vivi_alloc(nbytes ? nbytes : 1);
	if (!out->data) { *err = "out of memory"; return false; }
	for (size_t i = 0; i < nbytes; i++)
		out->data[i] = (uint8_t)(values[2 * i] * 16 + values[2 * i + 1]);
	out->len = nbytes;
	return true;
}

/* ---- genes ---- */

bool genome_gene_encode(vivi_bytes *out, int id, int usertype, int codon_mode, int h,
	const uint8_t *data, size_t len, const char **err)
{
	*out = (vivi_bytes){ 0 };
	if (id < 0 || id > 255) { *err = "invalid id (byte 0..255)"; return false; }
	if (usertype < 0 || usertype > 255) { *err = "invalid type (byte 0..255)"; return false; }
	if (len > 65535) { *err = "gene too large (raw > 65535)"; return false; }
	if (!codon_mode && (h < 3 || h > 12)) {
		*err = "invalid h for gene (integer 3..12)";
		return false;
	}
	vivi_bytes payload = { 0 };
	int packedlen;
	if (!codon_mode) {
		dna_opts dopts = { h, 0.05 };
		if (!dna_encode(&payload, data, len, &dopts, err)) return false;
	} else {
		size_t nv = len * 2;
		while (nv % 4 != 0) nv++;
		int *vals = vivi_alloc((nv ? nv : 1) * sizeof(int));
		if (!vals) { *err = "out of memory"; return false; }
		(void)genome_values_from_bytes(vals, nv, data, len);
		for (size_t i = len * 2; i < nv; i++) vals[i] = 0;
		if (!genome_pack_codons(&payload, vals, nv, err)) {
			vivi_dealloc(vals);
			return false;
		}
		vivi_dealloc(vals);
	}
	packedlen = (int)payload.len;
	if (packedlen > 65535) {
		vivi_bytes_free(&payload);
		*err = "gene too large (packed > 65535)";
		return false;
	}
	int typ = (usertype & 0xFE) | (codon_mode ? 1 : 0);
	uint8_t hdr[6] = { (uint8_t)id, (uint8_t)typ,
		(uint8_t)((len >> 8) & 255), (uint8_t)(len & 255),
		(uint8_t)((packedlen >> 8) & 255), (uint8_t)(packedlen & 255) };
	int hv[12];
	(void)genome_values_from_bytes(hv, 12, hdr, 6);
	vivi_bytes hdrp;
	if (!genome_pack_codons(&hdrp, hv, 12, err)) {
		vivi_bytes_free(&payload);
		return false;
	}
	uint32_t tag = genome_chaskey32(data, len);
	uint8_t tb[4] = { (uint8_t)((tag >> 24) & 255), (uint8_t)((tag >> 16) & 255),
		(uint8_t)((tag >> 8) & 255), (uint8_t)(tag & 255) };
	int tv[8];
	(void)genome_values_from_bytes(tv, 8, tb, 4);
	vivi_bytes tagp;
	if (!genome_pack_codons(&tagp, tv, 8, err)) {
		vivi_bytes_free(&payload);
		vivi_bytes_free(&hdrp);
		return false;
	}
	size_t total = 3 + 9 + payload.len + 6 + 3;
	out->data = vivi_alloc(total);
	if (!out->data) {
		vivi_bytes_free(&payload);
		vivi_bytes_free(&hdrp);
		vivi_bytes_free(&tagp);
		*err = "out of memory";
		return false;
	}
	size_t pos = 0;
	memcpy(out->data + pos, genome_PROM, 3); pos += 3;
	memcpy(out->data + pos, hdrp.data, hdrp.len); pos += hdrp.len;
	memcpy(out->data + pos, payload.data, payload.len); pos += payload.len;
	memcpy(out->data + pos, tagp.data, tagp.len); pos += tagp.len;
	memcpy(out->data + pos, genome_TERM, 3);
	out->len = total;
	vivi_bytes_free(&payload);
	vivi_bytes_free(&hdrp);
	vivi_bytes_free(&tagp);
	return true;
}

static int gene_parse_at(genome_gene *g, const uint8_t *strand, size_t slen, size_t p)
{
	memset(g, 0, sizeof(*g));
	if (p + 21 > slen) return 0;
	int hv[12];
	(void)(void)genome_unpack_codons(hv, 12, strand + p + 3, slen - p - 3, nullptr);
	int id = hv[0] * 16 + hv[1];
	int typ = hv[2] * 16 + hv[3];
	int rawlen = hv[4] * 4096 + hv[5] * 256 + hv[6] * 16 + hv[7];
	int packedlen = hv[8] * 4096 + hv[9] * 256 + hv[10] * 16 + hv[11];
	if (p + 21 + (size_t)packedlen > slen) return 0;
	size_t termi = p + 18 + (size_t)packedlen;
	if (memcmp(strand + termi, genome_TERM, 3) != 0) return 0;
	int modebit = typ & 1;
	size_t payload_off = p + 12;
	size_t payload_len = (size_t)packedlen;
	int tv[8];
	(void)genome_unpack_codons(tv, 8, strand + payload_off + payload_len, slen, nullptr);
	uint32_t expected = ((uint32_t)(tv[0] * 16 + tv[1]) << 24)
		| ((uint32_t)(tv[2] * 16 + tv[3]) << 16)
		| ((uint32_t)(tv[4] * 16 + tv[5]) << 8)
		| (uint32_t)(tv[6] * 16 + tv[7]);
	uint8_t *data = nullptr;
	size_t data_len = 0;
	int crc_ok = 0;
	if (modebit == 0) {
		vivi_bytes raw;
		if (dna_decode(&raw, strand + payload_off, payload_len, nullptr)) {
			data_len = (raw.len <= (size_t)rawlen) ? raw.len : (size_t)rawlen;
			data = vivi_alloc(data_len ? data_len : 1);
			if (data) {
				if (data_len) memcpy(data, raw.data, data_len);
				crc_ok = (data_len == (size_t)rawlen)
					&& genome_chaskey32(data, data_len) == expected;
			}
			vivi_bytes_free(&raw);
		}
	} else {
		size_t nv = payload_len * 4 / 3;
		if (payload_len % 3 == 0 && nv % 4 == 0) {
			int *vals = vivi_alloc((nv ? nv : 1) * sizeof(int));
			uint8_t *raw = vivi_alloc((nv / 2) ? nv / 2 : 1);
			if (vals && raw) {
				(void)genome_unpack_codons(vals, nv, strand + payload_off, payload_len, nullptr);
				size_t bn = nv / 2;
				for (size_t i = 0; i < bn; i++)
					raw[i] = (uint8_t)(vals[2 * i] * 16 + vals[2 * i + 1]);
				data_len = ((size_t)rawlen <= bn) ? (size_t)rawlen : bn;
				data = vivi_alloc(data_len ? data_len : 1);
				if (data) {
					memcpy(data, raw, data_len);
					crc_ok = (data_len == (size_t)rawlen)
						&& genome_chaskey32(data, data_len) == expected;
				}
			}
			vivi_dealloc(vals);
			vivi_dealloc(raw);
		}
	}
	g->id = id;
	g->type = typ;
	g->codon = modebit;
	g->rawlen = rawlen;
	g->packedlen = packedlen;
	g->offset = p;
	g->size = 21 + (size_t)packedlen;
	g->crc_ok = crc_ok;
	g->tag = expected;
	if (crc_ok) {
		g->data = data;
		g->data_len = data_len;
	} else {
		vivi_dealloc(data);
		g->data = nullptr;
		g->data_len = 0;
	}
	return 1;
}

void genome_scan_free(genome_scan_result *r)
{
	if (!r) return;
	for (size_t i = 0; i < r->count; i++) vivi_dealloc(r->genes[i].data);
	vivi_dealloc(r->genes);
	r->genes = nullptr;
	r->count = 0;
}

static size_t find_prom(const uint8_t *strand, size_t slen, size_t from)
{
	for (size_t i = from; i + 3 <= slen; i++)
		if (strand[i] == genome_PROM[0] && strand[i + 1] == genome_PROM[1]
			&& strand[i + 2] == genome_PROM[2])
			return i;
	return slen;
}

bool genome_gene_scan(genome_scan_result *out, const uint8_t *strand, size_t slen,
	const char **err)
{
	memset(out, 0, sizeof(*out));
	(void)err;
	size_t cap = 0;
	size_t pos = 0;
	while (pos + 3 <= slen) {
		size_t p = find_prom(strand, slen, pos);
		if (p >= slen) break;
		genome_gene g;
		if (gene_parse_at(&g, strand, slen, p)) {
			if (out->count == cap) {
				size_t ncap = cap ? cap * 2 : 8;
				genome_gene *ng = vivi_alloc(ncap * sizeof(genome_gene));
				if (!ng) {
					genome_scan_free(out);
					*err = "out of memory";
					return false;
				}
				if (cap) memcpy(ng, out->genes, cap * sizeof(genome_gene));
				vivi_dealloc(out->genes);
				out->genes = ng;
				cap = ncap;
			}
			out->genes[out->count++] = g;
			pos = p + g.size;
		} else {
			pos = p + 1;
		}
	}
	return true;
}

bool genome_gene_read(genome_gene *out, const uint8_t *strand, size_t slen, int id,
	const char **err)
{
	memset(out, 0, sizeof(*out));
	genome_scan_result res;
	if (!genome_gene_scan(&res, strand, slen, err)) return false;
	for (size_t i = 0; i < res.count; i++) {
		genome_gene *g = &res.genes[i];
		if (g->crc_ok && (id < 0 || g->id == id)) {
			*out = *g;
			g->data = nullptr; /* move ownership */
			genome_scan_free(&res);
			return true;
		}
	}
	genome_scan_free(&res);
	*err = "no valid gene";
	return false;
}

void genome_set_base(uint8_t *strand, size_t n, size_t idx, int digit)
{
	if (idx / 4 >= n) return;
	size_t bytei = idx / 4;
	int sh = (int)(2 * (3 - idx % 4));
	uint8_t mask = (uint8_t) ~(3u << sh);
	strand[bytei] = (uint8_t)((strand[bytei] & mask) | (((uint8_t)(digit & 3)) << sh));
}











