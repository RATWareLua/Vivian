/* gene.c -- libFuzzer harness for gene scanning and reading */
#include <stdlib.h>
#include "vivi/organism.h"

static vivi_bytes g_seed_dense, g_seed_codon;
static int g_seed_done;

static void seed_cleanup(void)
{
	vivi_bytes_free(&g_seed_dense);
	vivi_bytes_free(&g_seed_codon);
}

static void seed_once(void)
{
	if (g_seed_done) return;
	g_seed_done = 1;
	static const uint8_t payload[] = "gene fuzz payload 0123456789";
	(void)genome_gene_encode(&g_seed_dense, 7, 0, 0, 3, payload,
		sizeof(payload) - 1, nullptr);
	(void)genome_gene_encode(&g_seed_codon, 9, 0, 1, 3, payload,
		sizeof(payload) - 1, nullptr);
	atexit(seed_cleanup);
}

static void read_one(const uint8_t *data, size_t size)
{
	genome_gene g;
	if (genome_gene_read(&g, data, size, -1, nullptr))
		vivi_dealloc(g.data);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
	seed_once();
	if (g_seed_dense.data) read_one(g_seed_dense.data, g_seed_dense.len);
	if (g_seed_codon.data) read_one(g_seed_codon.data, g_seed_codon.len);
	genome_scan_result res;
	if (genome_gene_scan(&res, data, size, nullptr)) genome_scan_free(&res);
	read_one(data, size);
	return 0;
}
