/* density.c -- packing density: data bits per nucleotide.
 *
 * Alphabet ceiling: 4 bases = 2 bits -> 2.0 bits/nt.
 * Everything above 2.0 is impossible; everything below it is overhead:
 *   - codec expansion (homopolymer-breaking fallback + GC fill)
 *   - gene framing (PROM/header/tag/TERM = 84 nt per gene)
 *   - centromere + telomeres (per chromosome)
 * Run it: clang -std=c23 -O2 -Iinclude -o density.exe test/density.c src\*.c
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "vivi/vivi.h"
#include "vivi/dna.h"
#include "vivi/chromosome.h"

static uint32_t st = 0x2463u;
static uint32_t rnd(void)
{
	st ^= st << 13;
	st ^= st >> 17;
	st ^= st << 5;
	return st;
}

typedef void (*fill_fn)(uint8_t *, size_t);

static void fill_random(uint8_t *b, size_t n)
{
	for (size_t i = 0; i < n; i++) b[i] = (uint8_t)(rnd() >> 24);
}

static void fill_text(uint8_t *b, size_t n)
{
	static const char T[] = "METABOLISM:FAST;BRAIN:CURIOUS;SPEED:12;COLOR:RED; 0123456789 ";
	for (size_t i = 0; i < n; i++) b[i] = (uint8_t)T[i % (sizeof(T) - 1)];
}

static void fill_zeros(uint8_t *b, size_t n)
{
	memset(b, 0, n);
}

static const struct { fill_fn fill; const char *name; } PATS[] = {
	{ fill_random, "random" },
	{ fill_text, "text  " },
	{ fill_zeros, "zeros " },
};

static void codec_table(void)
{
	const size_t sizes[] = { 64, 256, 1024, 4096, 16384, 65536, 0 };
	printf("== codec only (dna_encode, h=3, gc_eps=0.05): bits/nt ==\n");
	printf("%-8s", "bytes");
	for (size_t p = 0; p < 3; p++) printf("%10s", PATS[p].name);
	printf("\n");
	for (size_t si = 0; si < sizeof(sizes) / sizeof(sizes[0]); si++) {
		size_t n = sizes[si];
		uint8_t *data = vivi_alloc(n);
		printf("%-8zu", n);
		for (size_t p = 0; p < 3; p++) {
			PATS[p].fill(data, n);
			vivi_bytes strand = { 0 };
			dna_opts o = { 3, 0.05 };
			const char *err = nullptr;
			if (!dna_encode(&strand, data, n, &o, &err)) {
				printf("      FAIL");
				continue;
			}
			printf("%10.3f", (double)(8 * n) / ((double)strand.len * 4.0));
			vivi_bytes_free(&strand);
		}
		printf("\n");
		vivi_dealloc(data);
	}
	printf("  (alphabet ceiling = 2.000; header is 8 bytes inside)\n\n");
}

static void chromosome_row(size_t n, fill_fn fill, const char *pname)
{
	uint8_t *data = vivi_alloc(n);
	fill(data, n);
	chr_opts co = { 1024, 0, 3, 4, 0, 0, 0 };
	const char *err = nullptr;
	vivi_bytes strand = { 0 };
	if (!chr_encode(&strand, 1, data, n, &co, &err)) {
		printf("%s %zuB: chr_encode FAIL: %s\n", pname, n, err);
		vivi_dealloc(data);
		return;
	}
	chr_record rec;
	if (!chr_parse(&rec, strand.data, strand.len, -1, &err)) {
		printf("%s %zuB: chr_parse FAIL: %s\n", pname, n, err);
		vivi_bytes_free(&strand);
		vivi_dealloc(data);
		return;
	}
	size_t nt = strand.len * 4;
	size_t telo_nt = 2 * rec.telomere_bytes * 4;
	size_t genes_nt = 0, payload_nt = 0;
	for (size_t i = 0; i < rec.gene_count; i++) {
		genes_nt += rec.genes[i].size * 4;
		payload_nt += (size_t)rec.genes[i].packedlen * 4;
	}
	size_t cen_nt = nt - telo_nt - genes_nt;
	size_t frame_nt = genes_nt - payload_nt;
	long expand_nt = (long)payload_nt - (long)(4 * n);
	double gc = dna_gc_content(strand.data, strand.len);
	size_t run = dna_max_homopolymer(strand.data, strand.len);
	printf("%s %7zuB nt=%-8zu eff=%5.3f gc=%.3f run=%zu | telo=%zu cen=%-4zu frame=%-6zu expand=%+ld nt (%+.1f%%)\n",
		pname, n, nt, (double)(8 * n) / (double)nt, gc, run,
		telo_nt, cen_nt, frame_nt, expand_nt,
		4.0 * (double)expand_nt / (4.0 * (double)n));
	chr_record_free(&rec);
	vivi_bytes_free(&strand);
	vivi_dealloc(data);
}

static void chromosome_table(void)
{
	const size_t sizes[] = { 64, 256, 1024, 4096, 16384, 65536, 0 };
	printf("== full chromosome (gene_raw=1024, dense, h=3, units=4) ==\n");
	for (size_t si = 0; si < sizeof(sizes) / sizeof(sizes[0]); si++) {
		size_t n = sizes[si];
		for (size_t p = 0; p < 3; p++)
			chromosome_row(n, PATS[p].fill, PATS[p].name);
	}
	printf("  gene frame = 84 nt/gene (PROM 12 + header 36 + tag 24 + TERM 12)\n");
	printf("  'expand' = payload re-encoding cost from homopolymer/GC constraints\n");
}

int main(void)
{
	codec_table();
	chromosome_table();
	return 0;
}
