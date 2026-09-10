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
	vivi_channel_opts ch = { 0.01, 0.001, 0.001, 0.01, (uint32_t)size + 1u };
	vivi_read rd;
	if (vivi_channel_read(&rd, data, size, &ch, nullptr))
		vivi_bytes_free(&rd.strand);
	return 0;
}
