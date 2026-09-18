/* dna.c -- libFuzzer harness for the strand codec, ASCII view and channel */
#include <stdlib.h>
#include "vivi/organism.h"
#include "vivi/channel.h"

static vivi_bytes g_seed;
static int g_seed_done;

static void seed_cleanup(void) { vivi_bytes_free(&g_seed); }

static void seed_once(void)
{
	if (g_seed_done) return;
	g_seed_done = 1;
	static const uint8_t payload[] = "Vivian fuzz seed payload";
	dna_opts o = { 3, 0.05 };
	(void)dna_encode(&g_seed, payload, sizeof(payload) - 1, &o, nullptr);
	atexit(seed_cleanup);
}

static double frac(const uint8_t *data, size_t size, size_t idx, double scale)
{
	if (idx >= size) return 0.0;
	return ((double)data[idx] / 255.0) * scale;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
	seed_once();
	vivi_bytes out = { 0 };
	if (g_seed.data && dna_decode(&out, g_seed.data, g_seed.len, nullptr))
		vivi_bytes_free(&out);
	out = (vivi_bytes){ 0 };
	if (dna_decode(&out, data, size, nullptr)) vivi_bytes_free(&out);
	vivi_bytes ascii = { 0 };
	(void)dna_from_ascii(&ascii, (const char *)data, size, nullptr);
	vivi_bytes_free(&ascii);
	vivi_bytes comp = { 0 };
	if (dna_complement(&comp, data, size, nullptr)) vivi_bytes_free(&comp);
	comp = (vivi_bytes){ 0 };
	if (dna_reverse_complement(&comp, data, size, nullptr)) vivi_bytes_free(&comp);
	vivi_channel_opts ch = { 0.01, 0.001, 0.001, 0.01, (uint32_t)size + 1u,
		frac(data, size, 0, 0.05), frac(data, size, 1, 0.05),
		frac(data, size, 2, 0.02), frac(data, size, 3, 0.02),
		(uint32_t)(size > 4 ? data[4] % 8u : 0u),
		frac(data, size, 5, 0.02) };
	vivi_read rd;
	if (vivi_channel_read(&rd, data, size, &ch, nullptr))
		vivi_read_free(&rd);
	vivi_consensus_opts cc = { ch, 1u + (uint32_t)(size % 32), (int)(size & 1) };
	vivi_read cr;
	if (vivi_consensus_read(&cr, data, size, &cc, nullptr))
		vivi_read_free(&cr);
	return 0;
}
