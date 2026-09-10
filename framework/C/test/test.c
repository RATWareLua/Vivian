/* test.c -- the C port of the framework test suites.
 *
 * Mirrors the Lua test bodies: Chaskey-12 official vectors, dna codec
 * roundtrips/constraints/errors, gene wobble tolerance, chromosome
 * roundtrips, diploid repair, organism evolution + serialization.
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef _WIN32
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>
#endif
#include "vivi/organism.h"
#include "vivi/channel.h"
#include "vivi/pool.h"
#include "vivi/parity.h"

static int g_pass = 0, g_fail = 0;

#define CHECK(cond, msg) do { \
	if (cond) g_pass++; \
	else { g_fail++; fprintf(stderr, "FAIL: %s\n", (msg)); } \
} while (0)

/* With -DVIVI_TEST_TRACK every library allocation is counted, so a
 * non-zero balance at exit is a leak even on platforms without LSan. */
#ifdef VIVI_TEST_TRACK
static long g_live;
static void *track_alloc(size_t n)
{
	void *p = malloc(n ? n : 1);
	if (p) g_live++;
	return p;
}
static void track_free(void *p)
{
	if (p) g_live--;
	free(p);
}
#endif

static uint32_t trnd_state = 0xC0FFEEu;
static uint32_t trnd(void)
{
	trnd_state ^= trnd_state << 13;
	trnd_state ^= trnd_state >> 17;
	trnd_state ^= trnd_state << 5;
	return trnd_state;
}

static void rands(uint8_t *buf, size_t n)
{
	for (size_t i = 0; i < n; i++) buf[i] = (uint8_t)(trnd() & 255);
}

static const char *ERR;

/* ---------- chaskey-12 official vectors ---------- */
static const uint32_t VK[4] = { 0x33221100u, 0x77665544u, 0xbbaa9988u, 0xffeeddccu };
static const uint8_t VV[64][8] = {
	{ 0xdd, 0x3e, 0x18, 0x49, 0xd6, 0x82, 0x45, 0x55 },
	{ 0xed, 0x1d, 0xa8, 0x9e, 0xc9, 0x31, 0x79, 0xca },
	{ 0x98, 0xfe, 0x20, 0xa3, 0x43, 0xcd, 0x66, 0x6f },
	{ 0xf6, 0xf4, 0x18, 0xac, 0xdd, 0x7d, 0x9f, 0xa1 },
	{ 0x4c, 0xf0, 0x49, 0x60, 0x09, 0x99, 0x49, 0xf3 },
	{ 0x75, 0xc8, 0x32, 0x52, 0x65, 0x3d, 0x3b, 0x57 },
	{ 0x96, 0x4b, 0x04, 0x61, 0xfb, 0xe9, 0x22, 0x73 },
	{ 0x14, 0x1f, 0xa0, 0x8b, 0xbf, 0x39, 0x96, 0x36 },
	{ 0x41, 0x2d, 0x98, 0xed, 0x93, 0x6d, 0x4a, 0xb2 },
	{ 0xfb, 0x0d, 0x98, 0xbc, 0x70, 0xe3, 0x05, 0xf9 },
	{ 0x36, 0xf8, 0x8e, 0x1f, 0xda, 0x86, 0xc8, 0xab },
	{ 0x4d, 0x1a, 0x18, 0x15, 0x86, 0x8a, 0x5a, 0xa8 },
	{ 0x7a, 0x79, 0x12, 0xc1, 0x99, 0x9e, 0xae, 0x81 },
	{ 0x9c, 0xa1, 0x11, 0x37, 0xb4, 0xa3, 0x46, 0x01 },
	{ 0x79, 0x05, 0x14, 0x2f, 0x3b, 0xe7, 0x7e, 0x67 },
	{ 0x6a, 0x3e, 0xe3, 0xd3, 0x5c, 0x04, 0x33, 0x97 },
	{ 0xd1, 0x39, 0x70, 0xd7, 0xbe, 0x9b, 0x23, 0x50 },
	{ 0x32, 0xac, 0xd9, 0x14, 0xbf, 0xda, 0x3b, 0xc8 },
	{ 0x8a, 0x58, 0xd8, 0x16, 0xcb, 0x7a, 0x14, 0x83 },
	{ 0x03, 0xf4, 0xd6, 0x66, 0x38, 0xef, 0xad, 0x8d },
	{ 0xf9, 0x93, 0x22, 0x37, 0xff, 0x05, 0xe8, 0x31 },
	{ 0xf5, 0xfe, 0xdb, 0x13, 0x48, 0x62, 0xb4, 0x71 },
	{ 0x8b, 0xb5, 0x54, 0x86, 0xf3, 0x8d, 0x57, 0xea },
	{ 0x8a, 0x3a, 0xcb, 0x94, 0xb5, 0xad, 0x59, 0x1c },
	{ 0x7c, 0xe3, 0x70, 0x87, 0x23, 0xf7, 0x49, 0x5f },
	{ 0xf4, 0x2f, 0x3d, 0x2f, 0x40, 0x57, 0x10, 0xc2 },
	{ 0xb3, 0x93, 0x3a, 0x16, 0x7e, 0x56, 0x36, 0xac },
	{ 0x89, 0x9a, 0x79, 0x45, 0x42, 0x3a, 0x5e, 0x1b },
	{ 0x65, 0xe1, 0x2d, 0xf5, 0xa6, 0x95, 0xfa, 0xc8 },
	{ 0xb8, 0x24, 0x49, 0xd8, 0xc8, 0xa0, 0x6a, 0xe9 },
	{ 0xa8, 0x50, 0xdf, 0xba, 0xde, 0xfa, 0x42, 0x29 },
	{ 0xfd, 0x42, 0xc3, 0x9d, 0x08, 0xab, 0x71, 0xa0 },
	{ 0xb4, 0x65, 0xc2, 0x41, 0x26, 0x10, 0xbf, 0x84 },
	{ 0x89, 0xc4, 0xa9, 0xdd, 0xb5, 0x3e, 0x69, 0x91 },
	{ 0x5a, 0x9a, 0xf9, 0x1e, 0xb0, 0x95, 0xd3, 0x31 },
	{ 0x8e, 0x54, 0x91, 0x4c, 0x15, 0x1e, 0x46, 0xb0 },
	{ 0xfa, 0xb8, 0xab, 0x0b, 0x5b, 0xea, 0xae, 0xc6 },
	{ 0x60, 0xad, 0x90, 0x6a, 0xcd, 0x06, 0xc8, 0x23 },
	{ 0x6b, 0x1e, 0x6b, 0xc2, 0x42, 0x6d, 0xad, 0x17 },
	{ 0x90, 0x32, 0x8f, 0xd2, 0x59, 0x88, 0x9a, 0x8f },
	{ 0xf0, 0xf7, 0x81, 0x5e, 0xe6, 0xf3, 0xd5, 0x16 },
	{ 0x97, 0xe7, 0xe2, 0xce, 0xbe, 0xa8, 0x26, 0xb8 },
	{ 0xb0, 0xfa, 0x18, 0x45, 0xf7, 0x2a, 0x76, 0xd6 },
	{ 0xa4, 0x68, 0xbd, 0xfc, 0xdf, 0x0a, 0xa9, 0xc7 },
	{ 0xda, 0x84, 0xe1, 0x13, 0x38, 0x38, 0x7d, 0xa7 },
	{ 0xb3, 0x0d, 0x5e, 0xad, 0x8e, 0x39, 0xf2, 0xbc },
	{ 0x17, 0x8a, 0x43, 0xd2, 0xa0, 0x08, 0x50, 0x3e },
	{ 0x6d, 0xfa, 0xa7, 0x05, 0xa8, 0xa0, 0x6c, 0x70 },
	{ 0xaa, 0x04, 0x7f, 0x07, 0xc5, 0xae, 0x8d, 0xb4 },
	{ 0x30, 0x5b, 0xbb, 0x42, 0x0c, 0x5d, 0x5e, 0xcc },
	{ 0x08, 0x32, 0x80, 0x31, 0x59, 0x75, 0x0f, 0x49 },
	{ 0x90, 0x80, 0x25, 0x4f, 0xb7, 0x9b, 0xab, 0x1a },
	{ 0x61, 0xc2, 0x85, 0xca, 0x24, 0x57, 0x74, 0xa4 },
	{ 0x2a, 0xae, 0x03, 0x5c, 0xfb, 0x61, 0xf9, 0x7a },
	{ 0xf5, 0x28, 0x90, 0x75, 0xc9, 0xab, 0x39, 0xe5 },
	{ 0xe6, 0x5c, 0x42, 0x37, 0x32, 0xda, 0xe7, 0x95 },
	{ 0x4b, 0x22, 0xcf, 0x0d, 0x9d, 0xa8, 0xde, 0x3d },
	{ 0x26, 0x26, 0xea, 0x2f, 0xa1, 0xf9, 0xab, 0xcf },
	{ 0xd1, 0xe1, 0x7e, 0x6e, 0xc4, 0xa8, 0x8d, 0xa6 },
	{ 0x16, 0x57, 0x44, 0x28, 0x27, 0xff, 0x64, 0x0a },
	{ 0xfd, 0x15, 0x5a, 0x40, 0xdf, 0x15, 0xf6, 0x30 },
	{ 0xff, 0xeb, 0x59, 0x6f, 0x29, 0x9f, 0x58, 0xb2 },
	{ 0xbe, 0x4e, 0xe4, 0xed, 0x39, 0x75, 0xdf, 0x87 },
	{ 0xfc, 0x7f, 0x9d, 0xf7, 0x99, 0x1b, 0x87, 0xbc },
};

static void test_chaskey(void)
{
	uint8_t m[64] = { 0 };
	int vecfail = 0;
	for (int i = 0; i < 64; i++) {
		/* official test: message = bytes 0..i-1 (length i) */
		uint32_t tag[4];
		genome_chaskey_tag(tag, m, (size_t)i, VK);
		uint8_t t8[8];
		for (int j = 0; j < 4; j++) {
			t8[j] = (uint8_t)((tag[0] >> (8 * j)) & 255u);
			t8[4 + j] = (uint8_t)((tag[1] >> (8 * j)) & 255u);
		}
		if (memcmp(t8, VV[i], 8) != 0) vecfail++;
		m[i] = (uint8_t)i; /* appended after the tag is computed */
	}
	CHECK(vecfail == 0, "chaskey-12 official vectors");
}

/* ---------- dna codec ---------- */
static const int hs[12] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

static void test_dna(void)
{
	uint8_t *zero1k = vivi_zalloc(1024, 1);
	uint8_t *ff1k = vivi_alloc(1024);
	memset(ff1k, 255, 1024);
	uint8_t *r1k = vivi_alloc(1024);
	rands(r1k, 1024);
	uint8_t *r8k = vivi_alloc(8192);
	rands(r8k, 8192);
	const uint8_t *cases[6] = { (const uint8_t *)"", (const uint8_t *)"a",
		(const uint8_t *)"Hello, DNA!", zero1k, ff1k, r8k };
	const size_t clens[6] = { 0, 1, 11, 1024, 1024, 8192 };

	for (int hi = 0; hi < 12; hi++) {
		int h = hs[hi];
		for (int ci = 0; ci < 6; ci++) {
			vivi_bytes s = { 0 }, back = { 0 };
			dna_opts o = { h, 0.05 };
			if (!dna_encode(&s, cases[ci], clens[ci], &o, &ERR)) {
				printf("FAIL: encode h=%d case=%d: %s\n", h, ci, ERR);
				g_fail++;
				continue;
			}
			if (!dna_decode(&back, s.data, s.len, &ERR)) {
				printf("FAIL: decode h=%d case=%d: %s\n", h, ci, ERR);
				g_fail++;
				vivi_bytes_free(&s);
				vivi_bytes_free(&back);
				continue;
			}
			CHECK(back.len == clens[ci] && (clens[ci] == 0 || memcmp(back.data, cases[ci], clens[ci]) == 0),
				"dna roundtrip");
			CHECK((int)dna_max_homopolymer(s.data, s.len) <= h, "dna homopolymer");
			double g = dna_gc_content(s.data, s.len);
			CHECK(g >= 0.45 && g <= 0.55, "dna gc band");
			CHECK(dna_validate(s.data, s.len, h, 0.05, &ERR), "dna validate");
			vivi_bytes_free(&s);
			vivi_bytes_free(&back);
		}
	}
	vivi_dealloc(zero1k);
	vivi_dealloc(ff1k);
	vivi_dealloc(r1k);
	vivi_dealloc(r8k);

	vivi_bytes r;
	CHECK(!dna_decode(&r, nullptr, 0, &ERR), "decode empty rejected");
	CHECK(!dna_decode(&r, (const uint8_t *)"\0\0", 2, &ERR), "invalid prefix rejected");
	CHECK(!dna_decode(&r, (const uint8_t *)"\0\0\0\1AAA", 7, &ERR), "truncated rejected");
	vivi_bytes full;
	(void)dna_encode(&full, (const uint8_t *)"truncate me please, payload", 27, nullptr, &ERR);
	CHECK(!dna_decode(&r, full.data, full.len / 2, &ERR), "truncated strand rejected");
	vivi_bytes_free(&full);

	vivi_bytes c1 = { 0 }, c2 = { 0 }, c3 = { 0 };
	(void)dna_encode(&c1, (const uint8_t *)"xy", 2, nullptr, &ERR);
	(void)dna_complement(&c2, c1.data, c1.len, &ERR);
	(void)dna_complement(&c3, c2.data, c2.len, &ERR);
	CHECK(c3.len == c1.len && memcmp(c3.data, c1.data, c1.len) == 0, "complement involution");
	vivi_bytes_free(&c2);
	(void)dna_reverse_complement(&c2, c1.data, c1.len, &ERR);
	CHECK(c2.len == c1.len && c2.data[0] == (uint8_t)(c1.data[c1.len - 1] ^ 0xFF),
		"reverse complement");
	vivi_bytes_free(&c2);
	(void)dna_to_ascii(&c2, c1.data, c1.len, -1, &ERR);
	vivi_bytes_free(&c3);
	(void)dna_from_ascii(&c3, (const char *)c2.data, c2.len, &ERR);
	CHECK(c3.len == c1.len && memcmp(c3.data, c1.data, c1.len) == 0, "ascii pack/unpack");
	CHECK(!dna_from_ascii(&r, "ACGX", 4, &ERR), "ascii rejects bad base");
	vivi_bytes_free(&c1);
	vivi_bytes_free(&c2);
	vivi_bytes_free(&c3);

	/* trailing partial groups are left-aligned (Lua reference bytes) */
	(void)dna_from_ascii(&c1, "ACG", 3, &ERR);
	CHECK(c1.len == 1 && c1.data[0] == 0x18, "from_ascii partial group ACG");
	vivi_bytes_free(&c1);
	(void)dna_from_ascii(&c1, "AC", 2, &ERR);
	CHECK(c1.len == 1 && c1.data[0] == 0x10, "from_ascii partial group AC");
	vivi_bytes_free(&c1);
}

/* ---------- genes ---------- */
static void test_genome(void)
{
	for (int mode = 0; mode < 2; mode++) {
		const size_t lens[5] = { 0, 1, 3, 100, 1000 };
		for (int ci = 0; ci < 5; ci++) {
			size_t len = lens[ci];
			uint8_t *data = vivi_alloc(len ? len : 1);
			rands(data, len);
			vivi_bytes g;
			if (!genome_gene_encode(&g, 7, 0, mode, 3, data, len, &ERR)) {
				printf("FAIL: gene encode mode=%d len=%zu: %s\n", mode, len, ERR);
				g_fail++;
				vivi_dealloc(data);
				continue;
			}
			genome_gene rd;
			CHECK(genome_gene_read(&rd, g.data, g.len, -1, &ERR), "gene roundtrip");
			CHECK(rd.id == 7 && rd.crc_ok && rd.data_len == len
				&& memcmp(rd.data, data, len) == 0, "gene data intact");
			vivi_dealloc(rd.data);
			genome_scan_result sc;
			(void)genome_gene_scan(&sc, g.data, g.len, &ERR);
			CHECK(sc.count == 1, "scan single gene");
			genome_scan_free(&sc);
			vivi_bytes_free(&g);
			vivi_dealloc(data);
		}
	}
	/* multi-gene scan */
	vivi_bytes g1, g2, g3;
	uint8_t big[300];
	rands(big, 300);
	(void)genome_gene_encode(&g1, 1, 0, 0, 3, (const uint8_t *)"alpha-data", 10, &ERR);
	(void)genome_gene_encode(&g2, 2, 0, 1, 3, (const uint8_t *)"beta-data", 9, &ERR);
	(void)genome_gene_encode(&g3, 3, 0, 0, 3, big, 300, &ERR);
	size_t clen = g1.len + g2.len + g3.len;
	uint8_t *chrom = vivi_alloc(clen);
	size_t pos = 0;
	memcpy(chrom + pos, g1.data, g1.len); pos += g1.len;
	memcpy(chrom + pos, g2.data, g2.len); pos += g2.len;
	memcpy(chrom + pos, g3.data, g3.len);
	genome_scan_result sc;
	(void)genome_gene_scan(&sc, chrom, clen, &ERR);
	CHECK(sc.count == 3, "scan finds 3 genes");
	if (sc.count == 3) {
		CHECK(sc.genes[0].id == 1 && sc.genes[1].id == 2 && sc.genes[2].id == 3,
			"scan order/ids");
		CHECK(sc.genes[0].data_len == 10 && memcmp(sc.genes[0].data, "alpha-data", 10) == 0,
			"gene 1 data");
		CHECK(sc.genes[1].data_len == 9 && memcmp(sc.genes[1].data, "beta-data", 9) == 0,
			"gene 2 data");
	}
	genome_scan_free(&sc);
	vivi_dealloc(chrom);
	vivi_bytes_free(&g1);
	vivi_bytes_free(&g2);
	vivi_bytes_free(&g3);
}

/* ---------- chromosome ---------- */
static void test_chromosome(void)
{
	const char *err;
	uint8_t data[500];
	rands(data, 500);
	vivi_bytes chr;
	chr_opts co = { 1024, 0, 3, 4, 0, 0 };
	CHECK(chr_encode(&chr, 1, data, 500, &co, &err), "chr encode");
	chr_record rec;
	CHECK(chr_parse(&rec, chr.data, chr.len, -1, &err), "chr parse");
	CHECK(rec.id == 1 && rec.ngenes == 1 && rec.generation == 0, "chr header");
	CHECK(rec.telomere_ok && rec.cen_ok, "chr structures");
	vivi_bytes back;
	CHECK(chr_read(&back, &rec, &err), "chr read");
	CHECK(back.len == 500 && memcmp(back.data, data, 500) == 0, "chr roundtrip");
	CHECK((int)dna_max_homopolymer(chr.data, chr.len) <= 3, "chr homopolymer");
	chr_record_free(&rec);
	vivi_bytes_free(&back);
	vivi_bytes_free(&chr);

	/* multi-gene */
	uint8_t big[5000];
	rands(big, 5000);
	chr_opts co2 = { 1024, 0, 3, 4, 0, 0 };
	(void)chr_encode(&chr, 2, big, 5000, &co2, &err);
	(void)chr_parse(&rec, chr.data, chr.len, -1, &err);
	CHECK(rec.ngenes == 5, "5 genes for 5000B @1024");
	(void)chr_read(&back, &rec, &err);
	CHECK(back.len == 5000 && memcmp(back.data, big, 5000) == 0, "multi-gene roundtrip");
	chr_record_free(&rec);
	vivi_bytes_free(&back);
	vivi_bytes_free(&chr);

	/* empty */
	(void)chr_encode(&chr, 3, (const uint8_t *)"", 0, &co, &err);
	(void)chr_parse(&rec, chr.data, chr.len, -1, &err);
	(void)chr_read(&back, &rec, &err);
	CHECK(back.len == 0, "empty chromosome");
	chr_record_free(&rec);
	vivi_bytes_free(&back);
	vivi_bytes_free(&chr);

	/* generation + damage detection */
	{
		uint8_t dd[300];
		rands(dd, 300);
		(void)chr_encode(&chr, 1, dd, 300, &co, &err);
		(void)chr_set_generation(&back, chr.data, chr.len, 7, &err);
		CHECK(back.len == chr.len, "set_generation preserves length");
		(void)chr_parse(&rec, back.data, back.len, -1, &err);
		CHECK(rec.generation == 7 && rec.id == 1, "generation rewritten");
		vivi_bytes_free(&back);
		(void)chr_read(&back, &rec, &err);
		CHECK(back.len == 300 && memcmp(back.data, dd, 300) == 0, "read after set_generation");
		chr_record_free(&rec);
		vivi_bytes_free(&back);
		vivi_bytes_free(&chr);
	}
}

/* ---------- prng ---------- */
static void test_prng(void)
{
	vivi_prng a, b;
	vivi_prng_init(&a, 42);
	vivi_prng_init(&b, 42);
	uint64_t x = vivi_prng_next64(&a);
	uint64_t y = vivi_prng_next64(&b);
	CHECK(x == y, "prng deterministic");
	vivi_prng_init(&b, 43);
	CHECK(vivi_prng_next64(&b) != x, "prng seed changes stream");
	CHECK(!vivi_prng_chance(&a, 0.0), "prng chance p=0");
	CHECK(vivi_prng_chance(&a, 1.0), "prng chance p=1");
	vivi_prng_init(&a, 0);
	uint64_t z1 = vivi_prng_next64(&a);
	uint64_t z2 = vivi_prng_next64(&a);
	CHECK(z1 != 0 && z1 != z2, "prng works from seed 0");
}

/* ---------- channel ---------- */
static void test_channel(void)
{
	uint8_t src[8];
	rands(src, sizeof(src));
	const char *err;
	vivi_read rd;
	vivi_channel_opts z = { 0, 0, 0, 0, 7 };
	CHECK(vivi_channel_read(&rd, src, sizeof(src), &z, &err), "channel identity");
	CHECK(!rd.dropped && rd.bases == 32 && rd.strand.len == 8
		&& memcmp(rd.strand.data, src, 8) == 0, "channel identity bytes");
	vivi_bytes_free(&rd.strand);

	vivi_channel_opts drop = { 0, 0, 0, 1.0, 7 };
	CHECK(vivi_channel_read(&rd, src, sizeof(src), &drop, &err), "channel drop");
	CHECK(rd.dropped && rd.strand.data == nullptr && rd.bases == 0, "channel dropped read");

	vivi_channel_opts sub = { 1.0, 0, 0, 0, 123 };
	CHECK(vivi_channel_read(&rd, src, sizeof(src), &sub, &err), "channel substitution");
	CHECK(!rd.dropped && rd.bases == 32 && rd.strand.len == 8, "substitution keeps shape");
	int changed = 0;
	for (int i = 0; i < 32; i++) {
		int a = (src[i / 4] >> (2 * (3 - i % 4))) & 3;
		int b = (rd.strand.data[i / 4] >> (2 * (3 - i % 4))) & 3;
		if (a != b) changed++;
	}
	CHECK(changed == 32, "every base substituted");
	vivi_read rd2;
	CHECK(vivi_channel_read(&rd2, src, sizeof(src), &sub, &err), "channel substitution 2");
	CHECK(rd2.strand.len == rd.strand.len
		&& memcmp(rd2.strand.data, rd.strand.data, rd.strand.len) == 0,
		"channel deterministic per seed");
	vivi_bytes_free(&rd2.strand);
	vivi_bytes_free(&rd.strand);

	vivi_channel_opts ins = { 0, 1.0, 0, 0, 4 };
	CHECK(vivi_channel_read(&rd, src, sizeof(src), &ins, &err), "channel insertion");
	CHECK(!rd.dropped && rd.bases == 64 && rd.strand.len == 16, "insertion doubles bases");
	vivi_bytes_free(&rd.strand);

	vivi_channel_opts del = { 0, 0, 1.0, 0, 4 };
	CHECK(vivi_channel_read(&rd, src, sizeof(src), &del, &err), "channel deletion");
	CHECK(!rd.dropped && rd.bases == 0 && rd.strand.data == nullptr, "deletion empties read");

	vivi_channel_opts bad = { 1.5, 0, 0, 0, 1 };
	CHECK(!vivi_channel_read(&rd, src, sizeof(src), &bad, &err), "channel rejects bad rate");
	(void)err;

	/* an identity channel is a lossless read of a real strand */
	uint8_t pl[64];
	rands(pl, sizeof(pl));
	vivi_bytes enc = { 0 }, dec = { 0 };
	dna_opts dopts = { 3, 0.05 };
	CHECK(dna_encode(&enc, pl, sizeof(pl), &dopts, &err), "channel encode");
	CHECK(vivi_channel_read(&rd, enc.data, enc.len, &z, &err), "channel identity on strand");
	CHECK(dna_decode(&dec, rd.strand.data, rd.strand.len, &err), "channel identity decode");
	CHECK(dec.len == sizeof(pl) && memcmp(dec.data, pl, sizeof(pl)) == 0,
		"channel identity roundtrip");
	vivi_bytes_free(&rd.strand);
	vivi_bytes_free(&enc);
	vivi_bytes_free(&dec);
}

/* ---------- pool / random access ---------- */
static void test_pool(void)
{
	const char *err;
	uint8_t payload[600];
	rands(payload, sizeof(payload));
	chr_opts co = { 64, 0, 3, 4, 0, 0 };
	vivi_bytes chr = { 0 };
	CHECK(chr_encode(&chr, 0, payload, sizeof(payload), &co, &err), "pool chromosome");
	vivi_pool pool = { 0 };
	CHECK(vivi_pool_add_chromosome(&pool, chr.data, chr.len, &err), "pool add chromosome");
	CHECK(pool.count == 10, "pool holds every gene");

	vivi_amp_opts ao = { 0.0, 0.0, 7, { 0.0, 0.0, 0.0, 0.0, 0 } };
	vivi_amp_result a1, a2;
	CHECK(vivi_pool_amplify(&a1, &pool, 3, &ao, &err), "pool amplify");
	CHECK(vivi_pool_amplify(&a2, &pool, 3, &ao, &err), "pool amplify again");
	CHECK(a1.id == 3 && !a1.read.dropped && !a2.read.dropped, "pool product on target");
	CHECK(a1.read.strand.len == a2.read.strand.len
		&& memcmp(a1.read.strand.data, a2.read.strand.data, a1.read.strand.len) == 0,
		"pool deterministic per seed");
	genome_gene g;
	CHECK(genome_gene_read(&g, a1.read.strand.data, a1.read.strand.len, 3, &err),
		"pool amplicon reads");
	CHECK(g.data_len == 64 && memcmp(g.data, payload + 3 * 64, 64) == 0, "pool on-target data");
	vivi_dealloc(g.data);
	vivi_bytes_free(&a1.read.strand);
	vivi_bytes_free(&a2.read.strand);

	vivi_amp_opts fail = { 1.0, 0.0, 7, { 0.0, 0.0, 0.0, 0.0, 0 } };
	CHECK(vivi_pool_amplify(&a1, &pool, 3, &fail, &err), "pool primer failure");
	CHECK(a1.read.dropped && a1.id == -1, "pool dropped product");

	vivi_amp_opts xtalk = { 0.0, 1.0, 7, { 0.0, 0.0, 0.0, 0.0, 0 } };
	CHECK(vivi_pool_amplify(&a1, &pool, 3, &xtalk, &err), "pool cross-talk");
	CHECK(a1.id >= 0 && a1.id != 3, "pool off-target id");
	vivi_bytes_free(&a1.read.strand);

	CHECK(!vivi_pool_amplify(&a1, &pool, 999, &ao, &err), "pool unknown target");

	vivi_pool_free(&pool);
	vivi_bytes_free(&chr);
}

/* ---------- parity / outer code ---------- */
static void test_parity(void)
{
	const char *err;
	enum { N = 4, M = 2, LEN = 32 };
	uint8_t databuf[N][LEN];
	uint8_t *data[N], *shards[N + M], *out[N];
	uint8_t present[N + M];
	for (int i = 0; i < N; i++) {
		rands(databuf[i], LEN);
		data[i] = databuf[i];
	}
	CHECK(vivi_parity_encode(shards, (const uint8_t *const *)data, N, M, LEN, &err),
		"parity encode");
	CHECK(memcmp(shards[0], data[0], LEN) == 0
		&& memcmp(shards[N - 1], data[N - 1], LEN) == 0, "parity is systematic");

	memset(present, 1, sizeof(present));
	CHECK(vivi_parity_decode(out, (const uint8_t *const *)shards, present, N, M, LEN, &err),
		"parity decode all present");
	int same = 1;
	for (int i = 0; i < N; i++) if (memcmp(out[i], data[i], LEN) != 0) same = 0;
	CHECK(same, "parity roundtrip");
	vivi_parity_release(out, N);

	int okall = 1;
	for (int e = 0; e < N + M; e++) {
		present[e] = 0;
		if (!vivi_parity_decode(out, (const uint8_t *const *)shards, present, N, M, LEN, &err))
			okall = 0;
		else {
			for (int i = 0; i < N; i++) if (memcmp(out[i], data[i], LEN) != 0) okall = 0;
			vivi_parity_release(out, N);
		}
		present[e] = 1;
	}
	CHECK(okall, "parity tolerates any single erasure");

	okall = 1;
	for (int a = 0; a < N + M; a++) {
		for (int b = a + 1; b < N + M; b++) {
			present[a] = 0;
			present[b] = 0;
			if (!vivi_parity_decode(out, (const uint8_t *const *)shards, present, N, M, LEN, &err))
				okall = 0;
			else {
				for (int i = 0; i < N; i++) if (memcmp(out[i], data[i], LEN) != 0) okall = 0;
				vivi_parity_release(out, N);
			}
			present[a] = 1;
			present[b] = 1;
		}
	}
	CHECK(okall, "parity tolerates any two erasures");

	present[0] = present[1] = present[2] = 0;
	CHECK(!vivi_parity_decode(out, (const uint8_t *const *)shards, present, N, M, LEN, &err),
		"parity rejects three erasures");
	memset(present, 1, sizeof(present));
	vivi_parity_release(shards, N + M);

	/* n = 1 is a repetition code: parity equals the single data shard */
	uint8_t one[8];
	rands(one, sizeof(one));
	uint8_t *d1[1] = { one };
	uint8_t *s2[2];
	CHECK(vivi_parity_encode(s2, (const uint8_t *const *)d1, 1, 1, sizeof(one), &err),
		"parity n=1 encode");
	CHECK(memcmp(s2[0], one, sizeof(one)) == 0 && memcmp(s2[1], one, sizeof(one)) == 0,
		"parity n=1 repeats data");
	vivi_parity_release(s2, 2);

	uint8_t *dummy[1] = { one };
	uint8_t *big[256];
	CHECK(!vivi_parity_encode(big, (const uint8_t *const *)dummy, 255, 1, 1, &err),
		"parity rejects bad geometry");
}

/* ---------- cells ---------- */
static void test_cells(void)
{
	uint8_t cdata[3000];
	rands(cdata, 3000);
	chr_opts co = { 1024, 0, 3, 4, 0, 0 };
	const char *err;
	vivi_cell *c;
	(void)cell_new(&c, 1, cdata, 3000, &co, &ERR);
	vivi_bytes dmg;
	(void)cell_damage_strand(&dmg, c->hom[0], c->hlen[0], 200, 42, &err);
	vivi_dealloc(c->hom[0]);
	c->hom[0] = dmg.data;
	c->hlen[0] = dmg.len;
	cell_report rep = cell_checkpoint(c);
	CHECK(rep.repaired + rep.structural >= 1, "repair happened");
	vivi_bytes data;
	CHECK(cell_read(&data, c, &rep, &err), "cell read after repair");
	CHECK(data.len == 3000 && memcmp(data.data, cdata, 3000) == 0, "repaired data intact");
	vivi_bytes_free(&data);

	/* damage transliteration pin (must match the Lua reference bytes) */
	{
		const uint8_t raw[4] = { 0x00, 0x11, 0x22, 0x33 };
		const uint8_t want[4] = { 0x20, 0x01, 0x22, 0x33 };
		vivi_bytes dd = { 0 };
		(void)cell_damage_strand(&dd, raw, 4, 2, 1, &err);
		CHECK(dd.len == 4 && memcmp(dd.data, want, 4) == 0, "damage matches Lua reference");
		vivi_bytes_free(&dd);
	}


	/* same-locus damage in both -> dead -> stem rescue */
	vivi_cell *stem;
	(void)cell_stem(&stem, 1, cdata, 3000, &co, &err);
	vivi_cell *c2;
	(void)cell_new(&c2, 1, cdata, 3000, &co, &err);
	(void)cell_damage_strand(&dmg, c2->hom[0], c2->hlen[0], 500, 777, &err);
	vivi_dealloc(c2->hom[0]);
	c2->hom[0] = dmg.data;
	c2->hlen[0] = dmg.len;
	(void)cell_damage_strand(&dmg, c2->hom[1], c2->hlen[1], 500, 777, &err);
	vivi_dealloc(c2->hom[1]);
	c2->hom[1] = dmg.data;
	c2->hlen[1] = dmg.len;
	cell_checkpoint(c2);
	CHECK(!cell_read(&data, c2, &rep, &err), "dead cell unreadable without stem");
	cell_attach_stem(c2, stem);
	CHECK(cell_read(&data, c2, &rep, &err), "stem rescue");
	CHECK(data.len == 3000 && memcmp(data.data, cdata, 3000) == 0, "rescued data intact");
	vivi_bytes_free(&data);
	CHECK(c2->generation == 0, "renewed generation 0");

	/* mitosis + Hayflick */
	vivi_cell *par, *dau;
	(void)cell_new(&par, 1, cdata, 3000, &co, &err);
	CHECK(cell_mitosis(&dau, par, &err), "mitosis");
	CHECK(dau->generation == 1 && par->generation == 1, "mitosis generation");
	chr_record rec;
	(void)chr_parse(&rec, dau->hom[0], dau->hlen[0], -1, &err);
	CHECK(rec.generation == 1, "centromere generation synced");
	chr_record_free(&rec);
	cell_free(dau);
	vivi_cell *aged;
	(void)cell_new(&aged, 1, cdata, 3000, &co, &err);
	aged->max_gen = 2;
	CHECK(cell_mitosis(&dau, aged, &err), "mitosis 1 ok");
	cell_free(dau);
	CHECK(cell_mitosis(&dau, aged, &err), "mitosis 2 ok");
	cell_free(dau);
	dau = nullptr;
	{
		const char *e2 = nullptr;
		int sen = !cell_mitosis(&dau, aged, &e2) && e2 != nullptr
			&& strcmp(e2, "senescent") == 0;
		CHECK(sen, "Hayflick limit");
		if (dau) cell_free(dau);
	}
	cell_free(aged);
	cell_free(par);
	cell_free(c);
	cell_free(c2);
	cell_free(stem);
}

/* ---------- organisms ---------- */
static void test_organisms(void)
{
	chr_opts o1 = { 1024, 0, 3, 4, 0, 0 };
	int ids[2] = { 0, 1 };
	chr_opts opts[2] = { o1, o1 };
	uint8_t d0[300], d1[200];
	rands(d0, 300);
	rands(d1, 200);
	const uint8_t *datas[2] = { d0, d1 };
	const size_t lens[2] = { 300, 200 };
	vivi_organism *org;
	CHECK(organism_new(&org, ids, opts, datas, lens, 2, 60, &ERR), "organism new");
	{
		int dup_ids[2] = { 0, 0 };
		vivi_organism *dup = nullptr;
		CHECK(!organism_new(&dup, dup_ids, opts, datas, lens, 2, 60, &ERR),
			"duplicate chr_id rejected");
	}
	vivi_bytes *rd = nullptr;
	cell_report rep;
	CHECK(organism_read(&rd, org, &rep, &ERR), "organism read");
	CHECK(rd[0].len == 300 && memcmp(rd[0].data, d0, 300) == 0, "chr 0 data");
	CHECK(rd[1].len == 200 && memcmp(rd[1].data, d1, 200) == 0, "chr 1 data");
	vivi_bytes_free_n(rd, 2);

	/* damage one chromosome in one homolog -> checkpoint heals */
	{
		vivi_bytes dmg;
		(void)cell_damage_strand(&dmg, org->hom[0][1], org->hlen[0][1], 150, 5, &ERR);
		vivi_dealloc(org->hom[0][1]);
		org->hom[0][1] = dmg.data;
		org->hlen[0][1] = dmg.len;
	}
	rep = organism_checkpoint(org);
	CHECK(rep.repaired + rep.structural >= 1, "organism checkpoint repaired");
	CHECK(organism_read(&rd, org, &rep, &ERR), "read after repair");
	CHECK(rd[1].len == 200 && memcmp(rd[1].data, d1, 200) == 0, "chr 1 repaired");
	vivi_bytes_free_n(rd, 2);

	/* stem + same-locus lethal damage + auto-renewal */
	vivi_organism *niche;
	(void)organism_stem(&niche, ids, opts, datas, lens, 2, 60, &ERR);
	vivi_organism *org2;
	(void)organism_new(&org2, ids, opts, datas, lens, 2, 60, &ERR);
	organism_attach_stem(org2, niche);
	{
		vivi_bytes dmg;
		(void)cell_damage_strand(&dmg, org2->hom[0][1], org2->hlen[0][1], 500, 777, &ERR);
		vivi_dealloc(org2->hom[0][1]);
		org2->hom[0][1] = dmg.data;
		org2->hlen[0][1] = dmg.len;
		(void)cell_damage_strand(&dmg, org2->hom[1][1], org2->hlen[1][1], 500, 777, &ERR);
		vivi_dealloc(org2->hom[1][1]);
		org2->hom[1][1] = dmg.data;
		org2->hlen[1][1] = dmg.len;
	}
	CHECK(organism_read(&rd, org2, &rep, &ERR), "read auto-renews from stem");
	CHECK(rd[1].len == 200 && memcmp(rd[1].data, d1, 200) == 0, "renewed data intact");
	CHECK(rep.renewed, "renewal reported");
	vivi_bytes_free_n(rd, 2);

	/* renewal from a stem with a different chromosome count must resize */
	{
		int one[1] = { 0 };
		chr_opts oneopts[1] = { o1 };
		const uint8_t *onedata[1] = { d0 };
		const size_t onelen[1] = { 300 };
		vivi_organism *small = nullptr;
		(void)organism_new(&small, one, oneopts, onedata, onelen, 1, 60, &ERR);
		organism_attach_stem(small, niche);
		organism_kill(small);
		vivi_bytes *rr = nullptr;
		int ok = organism_read(&rr, small, &rep, &ERR);
		CHECK(ok, "renew from differently sized stem");
		CHECK(ok && small->nchr == 2 && rr[1].len == 200
			&& memcmp(rr[1].data, d1, 200) == 0, "renewed karyotype matches stem");
		if (rr) vivi_bytes_free_n(rr, small->nchr);
		organism_free(small);
	}


	/* mutate */
	org_mut_result mu = { 0 };
	CHECK(organism_mutate(&mu, org, 10, 33, &ERR), "mutate");
	CHECK(mu.chr >= 0, "mutate returns chr");
	CHECK(organism_read(&rd, org, &rep, &ERR), "read after mutate");
	CHECK(rd[mu.chr].len == (mu.chr == 0 ? 300 : 200), "mutate keeps length");
	vivi_bytes_free_n(rd, 2);
	vivi_bytes_free(&mu.data);

	/* cross */
	vivi_organism *pb;
	uint8_t e0[300], e1[200];
	rands(e0, 300);
	rands(e1, 200);
	const uint8_t *edatas[2] = { e0, e1 };
	const size_t elens[2] = { 300, 200 };
	(void)organism_new(&pb, ids, opts, edatas, elens, 2, 60, &ERR);
	vivi_organism *child;
	int crok = organism_cross(&child, org, pb, 77, &ERR);
	CHECK(crok, "cross");
	if (!crok) { fprintf(stderr, "cross error: %s\n", ERR); }
	CHECK(organism_read(&rd, child, &rep, &ERR), "cross child readable");
	/* every child gene must come from one of the parents */
	int allin = 1;
	for (int i = 0; i < 2; i++) {
		chr_record rc;
		(void)chr_parse(&rc, child->hom[0][i], child->hlen[0][i], -1, &ERR);
		for (size_t k = 0; k < rc.gene_count; k++) {
			genome_gene *g = &rc.genes[k];
			int found = 0;
			for (int p = 0; p < 2 && !found; p++) {
				vivi_organism *parent = p ? pb : org;
				chr_record rp;
				if (!chr_parse(&rp, parent->hom[0][i], parent->hlen[0][i], -1, &ERR))
					continue;
				for (size_t q = 0; q < rp.gene_count && !found; q++)
					if (rp.genes[q].crc_ok && rp.genes[q].data_len == g->data_len
						&& memcmp(rp.genes[q].data, g->data, g->data_len) == 0)
						found = 1;
				chr_record_free(&rp);
			}
			if (!found) allin = 0;
		}
		chr_record_free(&rc);
	}
	CHECK(allin, "child genes are per-gene mosaics of the parents");
	vivi_bytes_free_n(rd, 2);

	/* serialize roundtrip */
	vivi_bytes ser;
	CHECK(organism_serialize(&ser, child, &ERR), "serialize");
	CHECK(ser.len >= 4 && memcmp(ser.data, "VIV1", 4) == 0, "VIV1 container magic");
	vivi_organism *loaded = nullptr;
	CHECK(!organism_deserialize(&loaded, ser.data, ser.len - 1, &ERR),
		"truncated container rejected");
	CHECK(organism_deserialize(&loaded, ser.data, ser.len, &ERR), "deserialize");
	CHECK(organism_read(&rd, loaded, &rep, &ERR), "deserialized reads");
	CHECK(rd[0].len == 300 && rd[1].len == 200, "serialize roundtrip data");
	vivi_bytes_free_n(rd, 2);
	organism_free(loaded);
	loaded = nullptr;
	CHECK(!organism_deserialize(&loaded, (const uint8_t *)"garbage!", 8, &ERR),
		"bad magic rejected");
	if (loaded) organism_free(loaded);
	vivi_bytes_free(&ser);
	organism_free(child);
	organism_free(pb);
	organism_free(niche);
	organism_free(org2);
	organism_free(org);
}

/* ---------- VIV14 container ---------- */
static void test_viv14(void)
{
	const char *err;
	chr_opts o1 = { 1024, 0, 3, 4, 0, 0 };
	int ids[1] = { 0 };
	uint8_t d0[300];
	rands(d0, 300);
	const uint8_t *datas[1] = { d0 };
	const size_t lens[1] = { 300 };
	vivi_organism *org = nullptr;
	CHECK(organism_new(&org, ids, &o1, datas, lens, 1, 7, &err), "viv14 organism");

	vivi_bytes v1 = { 0 }, v14 = { 0 };
	CHECK(organism_serialize(&v1, org, &err), "viv14 serialize v1");
	CHECK(organism_container_version(v1.data, v1.len) == 1, "viv14 sniffs v1");
	CHECK(organism_serialize14(&v14, org, &err), "viv14 serialize");
	CHECK(v14.len > 5 && memcmp(v14.data, "VIV14", 5) == 0, "viv14 magic");
	CHECK(organism_container_version(v14.data, v14.len) == 14, "viv14 sniffs v14");
	CHECK(organism_container_version((const uint8_t *)"junk", 4) == 0, "viv14 rejects junk");

	/* damage homolog 0 only: only VIV14 must keep both homologs */
	vivi_bytes dmg = { 0 };
	(void)cell_damage_strand(&dmg, org->hom[0][0], org->hlen[0][0], 40, 11, &err);
	vivi_dealloc(org->hom[0][0]);
	org->hom[0][0] = dmg.data;
	org->hlen[0][0] = dmg.len;
	size_t dl0 = org->hlen[0][0], dl1 = org->hlen[1][0];
	uint8_t *stored0 = vivi_alloc(dl0);
	uint8_t *stored1 = vivi_alloc(dl1);
	memcpy(stored0, org->hom[0][0], dl0);
	memcpy(stored1, org->hom[1][0], dl1);

	vivi_bytes s14 = { 0 };
	CHECK(organism_serialize14(&s14, org, &err), "viv14 serialize damaged pair");
	vivi_organism *loaded = nullptr;
	CHECK(organism_deserialize(&loaded, s14.data, s14.len, &err), "viv14 load pair");
	CHECK(loaded->nchr == 1 && loaded->max_gen == 7 && loaded->generation == 0,
		"viv14 metadata");
	CHECK(loaded->hlen[0][0] == dl0 && memcmp(loaded->hom[0][0], stored0, dl0) == 0
		&& loaded->hlen[1][0] == dl1 && memcmp(loaded->hom[1][0], stored1, dl1) == 0,
		"viv14 preserves both homologs");
	CHECK(memcmp(loaded->hom[0][0], loaded->hom[1][0], dl0) != 0,
		"viv14 keeps the divergence");
	vivi_bytes *rd = nullptr;
	cell_report rep;
	CHECK(organism_read(&rd, loaded, &rep, &err), "viv14 pair repairs");
	CHECK(rd[0].len == 300 && memcmp(rd[0].data, d0, 300) == 0, "viv14 repaired payload");
	vivi_bytes_free_n(rd, loaded->nchr);
	organism_free(loaded);

	/* VIV1 stores only homolog 0, and its parser rejects a strand with a
	 * damaged centromere -- exactly what VIV14 is there to preserve */
	vivi_bytes v1d = { 0 };
	CHECK(organism_serialize(&v1d, org, &err), "viv14 serialize damaged as v1");
	vivi_organism *old = nullptr;
	CHECK(!organism_deserialize(&old, v1d.data, v1d.len, &err),
		"v1 cannot load the damaged homolog");
	if (old) organism_free(old);

	vivi_organism *bad = nullptr;
	CHECK(!organism_deserialize(&bad, s14.data, s14.len - 1, &err), "viv14 truncated rejected");
	if (bad) organism_free(bad);
	uint8_t *buf = vivi_alloc(s14.len);
	memcpy(buf, s14.data, s14.len);
	buf[5] = 2;
	CHECK(!organism_deserialize(&bad, buf, s14.len, &err), "viv14 bad flags rejected");
	if (bad) organism_free(bad);
	vivi_dealloc(buf);

	CHECK(organism_deserialize(&bad, v1.data, v1.len, &err), "viv14 dispatcher reads v1");
	organism_free(bad);

	vivi_dealloc(stored0);
	vivi_dealloc(stored1);
	vivi_bytes_free(&s14);
	vivi_bytes_free(&v1d);
	vivi_bytes_free(&v14);
	vivi_bytes_free(&v1);
	organism_free(org);
}

/* ---------- VIV14N parity genes ---------- */
static void test_viv14n(void)
{
	const char *err;
	uint8_t payload[700];
	rands(payload, sizeof(payload));
	chr_opts co = { 128, 0, 3, 4, 0, 2 };   /* 6 data genes + 2 parity genes */
	vivi_bytes chr = { 0 };
	CHECK(chr_encode(&chr, 1, payload, sizeof(payload), &co, &err), "viv14n encode");
	chr_record rec;
	CHECK(chr_parse(&rec, chr.data, chr.len, -1, &err), "viv14n parse");
	CHECK(rec.cen_ok && rec.telomere_ok && rec.cen_version == 2, "viv14n centromere v2");
	CHECK(rec.ngenes == 8 && rec.parity == 2 && rec.rawlen == sizeof(payload),
		"viv14n header");
	CHECK(rec.gene_count == 8, "viv14n genes");
	int data_genes = 0, par_genes = 0;
	for (size_t k = 0; k < rec.gene_count; k++) {
		if (rec.genes[k].type & 2) par_genes++;
		else data_genes++;
	}
	CHECK(data_genes == 6 && par_genes == 2, "viv14n gene types");
	vivi_bytes back = { 0 };
	CHECK(chr_read(&back, &rec, &err), "viv14n read");
	CHECK(back.len == sizeof(payload) && memcmp(back.data, payload, sizeof(payload)) == 0,
		"viv14n roundtrip");
	vivi_bytes_free(&back);

	/* break two data genes (tag fails) and recover through parity */
	for (int k = 0; k < 2; k++)
		chr.data[rec.genes[k].offset + 12] ^= 0x01;   /* first payload base */
	chr_record rec2;
	CHECK(chr_parse(&rec2, chr.data, chr.len, -1, &err), "viv14n parse damaged");
	CHECK(rec2.gene_count == 8, "viv14n framing survives damage");
	int damaged = 0;
	for (size_t k = 0; k < rec2.gene_count; k++)
		if (!rec2.genes[k].crc_ok) damaged++;
	CHECK(damaged == 2, "viv14n two genes damaged");
	CHECK(chr_read(&back, &rec2, &err), "viv14n parity read");
	CHECK(back.len == sizeof(payload) && memcmp(back.data, payload, sizeof(payload)) == 0,
		"viv14n parity recovery");
	vivi_bytes_free(&back);

	/* the same corruption without parity is unrecoverable */
	chr_opts plain = { 128, 0, 3, 4, 0, 0 };
	vivi_bytes chr0 = { 0 };
	CHECK(chr_encode(&chr0, 1, payload, sizeof(payload), &plain, &err), "viv14n plain encode");
	chr_record rec0;
	CHECK(chr_parse(&rec0, chr0.data, chr0.len, -1, &err), "viv14n plain parse");
	for (int k = 0; k < 2; k++)
		chr0.data[rec0.genes[k].offset + 12] ^= 0x01;
	chr_record rec0b;
	CHECK(chr_parse(&rec0b, chr0.data, chr0.len, -1, &err), "viv14n plain reparse");
	CHECK(!chr_read(&back, &rec0b, &err), "viv14n plain cannot recover");
	chr_record_free(&rec0);
	chr_record_free(&rec0b);
	vivi_bytes_free(&chr0);

	/* generation rewrite keeps the VIV14N centromere */
	vivi_bytes gen = { 0 };
	CHECK(chr_set_generation(&gen, chr.data, chr.len, 9, &err), "viv14n set generation");
	chr_record recg;
	CHECK(chr_parse(&recg, gen.data, gen.len, -1, &err), "viv14n reparsed generation");
	CHECK(recg.cen_version == 2 && recg.generation == 9 && recg.parity == 2
		&& recg.rawlen == sizeof(payload), "viv14n generation preserved");
	chr_record_free(&recg);
	vivi_bytes_free(&gen);

	/* container VIV14N roundtrip keeps the parity option */
	int ids[1] = { 1 };
	const uint8_t *datas[1] = { payload };
	const size_t lens[1] = { sizeof(payload) };
	vivi_organism *org = nullptr;
	CHECK(organism_new(&org, ids, &co, datas, lens, 1, 60, &err), "viv14n organism");
	vivi_bytes ser = { 0 };
	CHECK(organism_serialize14n(&ser, org, &err), "viv14n serialize");
	CHECK(organism_container_version(ser.data, ser.len) == 141, "viv14n sniff");
	vivi_organism *loaded = nullptr;
	CHECK(organism_deserialize(&loaded, ser.data, ser.len, &err), "viv14n deserialize");
	CHECK(loaded->chr_opts[0].parity == 2, "viv14n option preserved");
	vivi_bytes *rd = nullptr;
	cell_report rep;
	CHECK(organism_read(&rd, loaded, &rep, &err), "viv14n organism read");
	CHECK(rd[0].len == sizeof(payload) && memcmp(rd[0].data, payload, sizeof(payload)) == 0,
		"viv14n organism payload");
	vivi_bytes_free_n(rd, loaded->nchr);
	organism_free(loaded);
	vivi_bytes_free(&ser);
	organism_free(org);

	chr_record_free(&rec);
	chr_record_free(&rec2);
	vivi_bytes_free(&chr);
}

int main(void)
{
#ifdef _WIN32
	_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
	_CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDERR);
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
#ifdef VIVI_TEST_TRACK
	vivi_set_allocator(track_alloc, track_free);
#endif
	test_chaskey();
	test_dna();
	test_genome();
	test_chromosome();
	test_prng();
	test_channel();
	test_pool();
	test_parity();
	test_cells();
	test_organisms();
	test_viv14();
	test_viv14n();
	dna_free_caches();
#ifdef VIVI_TEST_TRACK
	CHECK(g_live == 0, "no leaked allocations");
#endif
	printf("\npass=%d fail=%d\n", g_pass, g_fail);
	return g_fail ? 1 : 0;
}












