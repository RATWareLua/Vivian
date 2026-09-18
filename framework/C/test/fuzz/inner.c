#include <stdlib.h>
#include "vivi/genome.h"
#include "vivi/chromosome.h"

static vivi_bytes g_seed;
static int g_seed_done;

static void seed_cleanup(void)
{
	vivi_bytes_free(&g_seed);
}

static void seed_once(void)
{
	if (g_seed_done) return;
	g_seed_done = 1;
	static const uint8_t payload[] = "inner fuzz seed payload";
	chr_opts co = { 64, 0, 3, 4, 0, 0, 0, 16 };
	(void)chr_encode(&g_seed, 2, payload, sizeof(payload) - 1, &co, nullptr);
	atexit(seed_cleanup);
}

static void scan_one(const uint8_t *data, size_t size)
{
	genome_scan_result res;
	if (genome_gene_scan(&res, data, size, nullptr)) genome_scan_free(&res);
	chr_record rec;
	if (!chr_parse(&rec, data, size, -1, nullptr)) return;
	vivi_bytes out = { 0 };
	if (chr_read(&out, &rec, nullptr)) vivi_bytes_free(&out);
	chr_record_free(&rec);
}

static void inner_one(const uint8_t *data, size_t size)
{
	if (size < 1) return;
	int id = (int)data[0];
	size_t len = size - 1;
	if (len > 200) len = 200;
	int inner_m = 2 + (int)(data[size - 1] % 63);
	if (len + (size_t)inner_m > 255) return;
	vivi_bytes gene = { 0 };
	if (!genome_gene_encode_inner(&gene, id, 0, 0, 3, inner_m, data + 1, len, nullptr))
		return;
	genome_gene g;
	if (genome_gene_read(&g, gene.data, gene.len, -1, nullptr))
		vivi_dealloc(g.data);
	vivi_bytes_free(&gene);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
	seed_once();
	if (g_seed.data) scan_one(g_seed.data, g_seed.len);
	scan_one(data, size);
	inner_one(data, size);
	return 0;
}
