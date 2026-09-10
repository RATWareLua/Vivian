/* sim.c -- channel experiments: encode -> damage -> repair -> CSV (hosted).
 *
 * Runs a diploid, single-chromosome organism through the error channel
 * `--trials` times and reports how often organism_read() recovers the
 * payload exactly. One CSV row per invocation; loop over parameters to
 * build sweeps. See docs/research.md.
 */
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vivi/organism.h"
#include "vivi/channel.h"

static void usage(void)
{
	printf(
		"vivi_sim -- synthesis/sequencing channel experiment\n"
		"usage: vivi_sim (--in FILE | --size N) [options]\n"
		"  --in FILE      payload file (max 16 MiB)\n"
		"  --size N       synthetic payload of N random bytes\n"
		"  --out FILE     write the CSV row to FILE instead of stdout\n"
		"  --trials N     number of trials (default 1000)\n"
		"  --p-sub P      substitution probability per base (default 0)\n"
		"  --p-ins P      insertion probability per base (default 0)\n"
		"  --p-del P      deletion probability per base (default 0)\n"
		"  --p-drop P     whole-read loss probability (default 0)\n"
		"  --seed S       base seed (default 1)\n"
		"  --h H          dense homopolymer limit (default 3)\n"
		"  --gene-raw R   raw bytes per gene (default 1024)\n"
		"  --units U      telomere repeat units (default 4)\n"
		"  --header       print the CSV header first\n"
		"columns: p_sub,p_ins,p_del,p_drop,trials,success,wrong,failed,dropped,\n"
		"         rate,avg_repaired,avg_structural,avg_dead,avg_anomaly\n");
}

static int parse_u32(const char *s, uint32_t *out)
{
	char *end = nullptr;
	if (!s || !*s) return 0;
	unsigned long long v = strtoull(s, &end, 0);
	if (*end || v > 0xFFFF'FFFFull) return 0;
	*out = (uint32_t)v;
	return 1;
}

static int parse_double(const char *s, double *out)
{
	char *end = nullptr;
	if (!s || !*s) return 0;
	double v = strtod(s, &end);
	if (*end || !(v >= 0.0) || v > 1.0) return 0;
	*out = v;
	return 1;
}

int main(int argc, char **argv)
{
	const char *in_path = nullptr, *out_path = nullptr;
	uint32_t size = 0, trials = 1000, seed = 1;
	int have_size = 0, header = 0;
	double p_sub = 0.0, p_ins = 0.0, p_del = 0.0, p_drop = 0.0;
	chr_opts co = { 1024, 0, 3, 4, 0 };

	for (int i = 1; i < argc; i++) {
		const char *a = argv[i];
		const char *v = (i + 1 < argc) ? argv[i + 1] : nullptr;
		if (!strcmp(a, "--in") && v) { in_path = v; i++; }
		else if (!strcmp(a, "--out") && v) { out_path = v; i++; }
		else if (!strcmp(a, "--size") && parse_u32(v, &size)) { have_size = 1; i++; }
		else if (!strcmp(a, "--trials") && parse_u32(v, &trials)) { i++; }
		else if (!strcmp(a, "--seed") && parse_u32(v, &seed)) { i++; }
		else if (!strcmp(a, "--p-sub") && parse_double(v, &p_sub)) { i++; }
		else if (!strcmp(a, "--p-ins") && parse_double(v, &p_ins)) { i++; }
		else if (!strcmp(a, "--p-del") && parse_double(v, &p_del)) { i++; }
		else if (!strcmp(a, "--p-drop") && parse_double(v, &p_drop)) { i++; }
		else if (!strcmp(a, "--h")) {
			uint32_t h;
			if (!parse_u32(v, &h) || h < 3 || h > 12) { usage(); return 2; }
			co.h = (int)h; i++;
		}
		else if (!strcmp(a, "--gene-raw")) {
			uint32_t g;
			if (!parse_u32(v, &g) || g < 16 || g > 65535) { usage(); return 2; }
			co.gene_raw = (int)g; i++;
		}
		else if (!strcmp(a, "--units")) {
			uint32_t u;
			if (!parse_u32(v, &u) || u < 1) { usage(); return 2; }
			co.units = (int)u; i++;
		}
		else if (!strcmp(a, "--header")) { header = 1; }
		else { usage(); return 2; }
	}
	if ((!in_path && !have_size) || (in_path && have_size)) { usage(); return 2; }
	if (trials == 0) { usage(); return 2; }

	uint8_t *payload = nullptr;
	size_t plen = 0;
	if (in_path) {
		FILE *f = fopen(in_path, "rb");
		if (!f) { fprintf(stderr, "vivi_sim: cannot open %s\n", in_path); return 1; }
		fseek(f, 0, SEEK_END);
		long sz = ftell(f);
		fseek(f, 0, SEEK_SET);
		if (sz < 0 || (unsigned long)sz > 16u * 1024u * 1024u) {
			fprintf(stderr, "vivi_sim: file too large (max 16 MiB)\n");
			fclose(f);
			return 1;
		}
		plen = (size_t)sz;
		payload = vivi_alloc(plen ? plen : 1);
		if (!payload || (plen && fread(payload, 1, plen, f) != plen)) {
			fprintf(stderr, "vivi_sim: read error\n");
			vivi_dealloc(payload);
			fclose(f);
			return 1;
		}
		fclose(f);
	} else {
		plen = size;
		payload = vivi_alloc(plen ? plen : 1);
		if (!payload) { fprintf(stderr, "vivi_sim: out of memory\n"); return 1; }
		vivi_rng r;
		vivi_rng_init(&r, seed);
		for (size_t i = 0; i < plen; i++) payload[i] = (uint8_t)(vivi_rng_next(&r) & 255u);
	}

	const char *err = nullptr;
	vivi_organism *org = nullptr;
	int ids[1] = { 0 };
	chr_opts opts[1] = { co };
	const uint8_t *datas[1] = { payload };
	const size_t lens[1] = { plen };
	if (!organism_new(&org, ids, opts, datas, lens, 1, 60, &err)) {
		fprintf(stderr, "vivi_sim: organism_new: %s\n", err ? err : "?");
		vivi_dealloc(payload);
		return 1;
	}

	uint8_t *pristine[2];
	size_t prlen[2];
	for (int h = 0; h < 2; h++) {
		prlen[h] = org->hlen[h][0];
		pristine[h] = vivi_alloc(prlen[h] ? prlen[h] : 1);
		memcpy(pristine[h], org->hom[h][0], prlen[h]);
	}

	unsigned long long success = 0, wrong = 0, failed = 0, dropped = 0;
	unsigned long long sum_repaired = 0, sum_struct = 0, sum_dead = 0, sum_anom = 0;
	for (uint32_t t = 0; t < trials; t++) {
		for (int h = 0; h < 2; h++) {
			vivi_channel_opts ch = { p_sub, p_ins, p_del, p_drop,
				seed + t * 2u + (uint32_t)h + 1u };
			vivi_read rd;
			if (!vivi_channel_read(&rd, pristine[h], prlen[h], &ch, &err)) {
				fprintf(stderr, "vivi_sim: channel: %s\n", err ? err : "?");
				return 1;
			}
			vivi_dealloc(org->hom[h][0]);
			org->hom[h][0] = rd.strand.data;
			org->hlen[h][0] = rd.strand.len;
			if (rd.dropped) dropped++;
		}
		vivi_bytes *data = nullptr;
		cell_report rep;
		if (organism_read(&data, org, &rep, &err)) {
			if (data[0].len == plen
				&& (plen == 0 || memcmp(data[0].data, payload, plen) == 0))
				success++;
			else
				wrong++;
			vivi_bytes_free_n(data, org->nchr);
		} else {
			failed++;
		}
		sum_repaired += (unsigned long long)rep.repaired;
		sum_struct += (unsigned long long)rep.structural;
		sum_dead += (unsigned long long)rep.dead;
		sum_anom += (unsigned long long)rep.anomaly;
	}

	FILE *out = stdout;
	if (out_path && !(out = fopen(out_path, "wb"))) {
		fprintf(stderr, "vivi_sim: cannot write %s\n", out_path);
		return 1;
	}
	if (header)
		fprintf(out, "p_sub,p_ins,p_del,p_drop,trials,success,wrong,failed,dropped,"
			"rate,avg_repaired,avg_structural,avg_dead,avg_anomaly\n");
	double dtrials = (double)trials;
	fprintf(out, "%g,%g,%g,%g,%u,%llu,%llu,%llu,%llu,%.6f,%.4f,%.4f,%.4f,%.4f\n",
		p_sub, p_ins, p_del, p_drop, trials,
		success, wrong, failed, dropped, (double)success / dtrials,
		(double)sum_repaired / dtrials, (double)sum_struct / dtrials,
		(double)sum_dead / dtrials, (double)sum_anom / dtrials);
	if (out != stdout) fclose(out);

	for (int h = 0; h < 2; h++) vivi_dealloc(pristine[h]);
	organism_free(org);
	vivi_dealloc(payload);
	dna_free_caches();
	return 0;
}
