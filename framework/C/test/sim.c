/* sim.c -- channel experiments: encode -> damage -> repair -> CSV (hosted).
 *
 * Three modes share one simulator:
 *   whole   -- read the whole diploid chromosome through the repair machinery
 *   access  -- amplify a single gene id from an oligo pool ("PCR primer")
 *   library -- n data genes + m parity genes + K replicas each (outer code)
 * One CSV row per invocation; loop over parameters to build sweeps.
 * See docs/research.md.
 */
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vivi/organism.h"
#include "vivi/channel.h"
#include "vivi/pool.h"
#include "vivi/parity.h"

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
		"  --access       random-access mode: amplify one gene id per trial\n"
		"  --p-access P   primer failure probability (default 0)\n"
		"  --p-cross P    off-target amplification probability (default 0)\n"
		"  --p-primer P   on-strand primer-site dropout (Banshee, default 0)\n"
		"  --library      library mode: gene molecules + outer code\n"
		"  --replicas K   copies of every molecule (default 1)\n"
		"  --parity M     parity molecules over the data genes (default 0)\n"
		"  --header       print the CSV header first\n"
		"whole:   mode,p_sub,p_ins,p_del,p_drop,trials,success,wrong,failed,dropped,\n"
		"         rate,avg_repaired,avg_structural,avg_dead,avg_anomaly\n"
		"access:  mode,p_sub,p_ins,p_del,p_drop,p_access,p_cross,p_primer,trials,success,\n"
		"         wrong,failed,dropped,cross,rate\n"
		"library: mode,p_sub,p_ins,p_del,p_drop,replicas,parity,genes,trials,\n"
		"         success,wrong,failed,dropped,rate,avg_present\n");
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

/* library mode: n data genes, m parity molecules, K replicas each; a gene
 * counts as present when any replica survives the channel intact */
static int run_library(const uint8_t *payload, size_t plen, const chr_opts *co,
	uint32_t trials, uint32_t seed, uint32_t replicas, uint32_t parity_m,
	double p_sub, double p_ins, double p_del, double p_drop,
	int header, const char *out_path)
{
	const char *err = nullptr;
	size_t shard_len = (size_t)co->gene_raw;
	size_t n = (plen + shard_len - 1) / shard_len;
	size_t m = parity_m;
	if (n + m > 255) {
		fprintf(stderr, "vivi_sim: too many shards (n+m must be <= 255)\n");
		return 1;
	}
	uint8_t **data = vivi_alloc(n * sizeof(uint8_t *));
	uint8_t **shards = vivi_alloc((n + m) * sizeof(uint8_t *));
	uint8_t **mol = vivi_alloc((n + m) * sizeof(uint8_t *));
	size_t *mlen = vivi_alloc((n + m) * sizeof(size_t));
	uint8_t **recv = vivi_zalloc(n + m, sizeof(uint8_t *));
	uint8_t *present = vivi_zalloc(n + m, 1);
	uint8_t **out = vivi_alloc(n * sizeof(uint8_t *));
	if (!data || !shards || !mol || !mlen || !recv || !present || !out) {
		fprintf(stderr, "vivi_sim: out of memory\n");
		return 1;
	}
	for (size_t i = 0; i < n; i++) {
		data[i] = vivi_zalloc(shard_len, 1);
		size_t off = i * shard_len;
		size_t take = (plen - off < shard_len) ? plen - off : shard_len;
		memcpy(data[i], payload + off, take);
	}
	if (!vivi_parity_encode(shards, (const uint8_t *const *)data, n, m, shard_len, &err)) {
		fprintf(stderr, "vivi_sim: parity encode: %s\n", err ? err : "?");
		return 1;
	}
	for (size_t s = 0; s < n + m; s++) {
		size_t rawlen = shard_len;
		if (s < n) {
			size_t off = s * shard_len;
			rawlen = (plen - off < shard_len) ? plen - off : shard_len;
		}
		vivi_bytes g = { 0 };
		if (!genome_gene_encode(&g, (int)s, 0, 0, co->h ? co->h : 3,
			shards[s], rawlen, &err)) {
			fprintf(stderr, "vivi_sim: gene encode: %s\n", err ? err : "?");
			return 1;
		}
		mol[s] = g.data;
		mlen[s] = g.len;
	}
	vivi_parity_release(shards, n + m);
	vivi_dealloc(shards);
	for (size_t i = 0; i < n; i++) vivi_dealloc(data[i]);
	vivi_dealloc(data);

	long long success = 0, wrong = 0, failed = 0, dropped = 0;
	unsigned long long sum_present = 0;
	for (uint32_t t = 0; t < trials; t++) {
		for (size_t s = 0; s < n + m; s++) present[s] = 0;
		size_t present_count = 0;
		for (size_t s = 0; s < n + m; s++) {
			for (uint32_t r = 0; r < replicas; r++) {
				vivi_channel_opts ch = { p_sub, p_ins, p_del, p_drop,
					seed + t * 977u + (uint32_t)s * 31u + r + 1u };
				vivi_read rd;
				if (!vivi_channel_read(&rd, mol[s], mlen[s], &ch, &err)) {
					fprintf(stderr, "vivi_sim: channel: %s\n", err ? err : "?");
					return 1;
				}
				if (rd.dropped) {
					dropped++;
					continue;
				}
				genome_gene g;
				if (genome_gene_read(&g, rd.strand.data, rd.strand.len, (int)s, &err)) {
					if (!present[s]) {
						recv[s] = vivi_zalloc(shard_len, 1);
						memcpy(recv[s], g.data, g.data_len);
						present[s] = 1;
						present_count++;
					}
					vivi_dealloc(g.data);
				}
				vivi_bytes_free(&rd.strand);
				if (present[s]) break;
			}
		}
		sum_present += (unsigned long long)present_count;
		if (present_count < n) {
			failed++;
		} else if (!vivi_parity_decode(out, (const uint8_t *const *)recv, present,
			n, m, shard_len, &err)) {
			fprintf(stderr, "vivi_sim: parity decode: %s\n", err ? err : "?");
			return 1;
		} else {
			int match = 1;
			size_t off = 0;
			for (size_t i = 0; i < n; i++) {
				size_t take = (plen - off < shard_len) ? plen - off : shard_len;
				if (memcmp(out[i], payload + off, take) != 0) match = 0;
				off += take;
			}
			if (match) success++;
			else wrong++;
			vivi_parity_release(out, n);
		}
		for (size_t s = 0; s < n + m; s++) {
			vivi_dealloc(recv[s]);
			recv[s] = nullptr;
		}
	}

	FILE *fout = stdout;
	if (out_path && !(fout = fopen(out_path, "wb"))) {
		fprintf(stderr, "vivi_sim: cannot write %s\n", out_path);
		return 1;
	}
	if (header)
		fprintf(fout, "mode,p_sub,p_ins,p_del,p_drop,replicas,parity,genes,trials,"
			"success,wrong,failed,dropped,rate,avg_present\n");
	double dt = (double)trials;
	fprintf(fout, "library,%g,%g,%g,%g,%u,%u,%zu,%u,%lld,%lld,%lld,%lld,%.6f,%.4f\n",
		p_sub, p_ins, p_del, p_drop, replicas, parity_m, n, trials,
		success, wrong, failed, dropped, (double)success / dt,
		(double)sum_present / dt);
	if (fout != stdout) fclose(fout);

	for (size_t s = 0; s < n + m; s++) vivi_dealloc(mol[s]);
	vivi_dealloc(mol);
	vivi_dealloc(mlen);
	vivi_dealloc(recv);
	vivi_dealloc(present);
	vivi_dealloc(out);
	return 0;
}

int main(int argc, char **argv)
{
	const char *in_path = nullptr, *out_path = nullptr;
	uint32_t size = 0, trials = 1000, seed = 1;
	uint32_t replicas = 1, parity_m = 0;
	int have_size = 0, header = 0, access_mode = 0, library_mode = 0;
	double p_sub = 0.0, p_ins = 0.0, p_del = 0.0, p_drop = 0.0;
	double p_access = 0.0, p_cross = 0.0, p_primer = 0.0;
	chr_opts co = { 1024, 0, 3, 4, 0, 0, 0 };

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
		else if (!strcmp(a, "--p-access") && parse_double(v, &p_access)) { i++; }
		else if (!strcmp(a, "--p-cross") && parse_double(v, &p_cross)) { i++; }
		else if (!strcmp(a, "--p-primer") && parse_double(v, &p_primer)) { i++; }
		else if (!strcmp(a, "--access")) { access_mode = 1; }
		else if (!strcmp(a, "--library")) { library_mode = 1; }
		else if (!strcmp(a, "--replicas") && parse_u32(v, &replicas) && replicas >= 1) { i++; }
		else if (!strcmp(a, "--parity") && parse_u32(v, &parity_m)) { i++; }
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
	if (access_mode && library_mode) { usage(); return 2; }

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
		vivi_prng r;
		vivi_prng_init(&r, seed);
		for (size_t i = 0; i < plen; i++) payload[i] = (uint8_t)(vivi_prng_next(&r) & 255u);
	}
	if ((access_mode || library_mode) && plen == 0) {
		fprintf(stderr, "vivi_sim: --access/--library need a non-empty payload\n");
		vivi_dealloc(payload);
		return 1;
	}
	if (library_mode) {
		int rc = run_library(payload, plen, &co, trials, seed, replicas, parity_m,
			p_sub, p_ins, p_del, p_drop, header, out_path);
		vivi_dealloc(payload);
		return rc;
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

	long long success = 0, wrong = 0, failed = 0, dropped = 0, cross = 0;
	long long sum_repaired = 0, sum_struct = 0, sum_dead = 0, sum_anom = 0;

	if (access_mode) {
		vivi_pool pool = { 0 };
		if (!vivi_pool_add_chromosome(&pool, pristine[0], prlen[0], &err)) {
			fprintf(stderr, "vivi_sim: pool: %s\n", err ? err : "?");
			return 1;
		}
		for (uint32_t t = 0; t < trials; t++) {
			vivi_prng pick_rng;
			vivi_prng_init(&pick_rng, (uint64_t)seed + t + 1u);
			int gid = pool.ids[vivi_prng_next(&pick_rng) % (uint32_t)pool.count];
			vivi_amp_opts ao = { p_access, p_cross, seed + t + 1u,
				{ p_sub, p_ins, p_del, 0.0, 0 }, p_primer };
			vivi_amp_result ar;
			if (!vivi_pool_amplify(&ar, &pool, gid, &ao, &err)) {
				fprintf(stderr, "vivi_sim: amplify: %s\n", err ? err : "?");
				return 1;
			}
			if (ar.read.dropped) {
				dropped++;
				continue;
			}
			if (ar.id != gid) cross++;
			genome_gene g;
			if (genome_gene_read(&g, ar.read.strand.data, ar.read.strand.len, gid, &err)) {
				size_t off = (size_t)gid * (size_t)co.gene_raw;
				size_t exp = (off < plen)
					? ((plen - off < (size_t)co.gene_raw) ? plen - off : (size_t)co.gene_raw)
					: 0;
				if (g.data_len == exp && (exp == 0 || memcmp(g.data, payload + off, exp) == 0))
					success++;
				else
					wrong++;
				vivi_dealloc(g.data);
			} else {
				failed++;
			}
			vivi_bytes_free(&ar.read.strand);
		}
		FILE *out = stdout;
		if (out_path && !(out = fopen(out_path, "wb"))) {
			fprintf(stderr, "vivi_sim: cannot write %s\n", out_path);
			return 1;
		}
		if (header)
			fprintf(out, "mode,p_sub,p_ins,p_del,p_drop,p_access,p_cross,p_primer,trials,"
				"success,wrong,failed,dropped,cross,rate\n");
		fprintf(out, "access,%g,%g,%g,%g,%g,%g,%g,%u,%lld,%lld,%lld,%lld,%lld,%.6f\n",
			p_sub, p_ins, p_del, p_drop, p_access, p_cross, p_primer, trials,
			success, wrong, failed, dropped, cross, (double)success / (double)trials);
		if (out != stdout) fclose(out);
		vivi_pool_free(&pool);
	} else {
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
			sum_repaired += (long long)rep.repaired;
			sum_struct += (long long)rep.structural;
			sum_dead += (long long)rep.dead;
			sum_anom += (long long)rep.anomaly;
		}
		FILE *out = stdout;
		if (out_path && !(out = fopen(out_path, "wb"))) {
			fprintf(stderr, "vivi_sim: cannot write %s\n", out_path);
			return 1;
		}
		if (header)
			fprintf(out, "mode,p_sub,p_ins,p_del,p_drop,trials,success,wrong,failed,dropped,"
				"rate,avg_repaired,avg_structural,avg_dead,avg_anomaly\n");
		double dt = (double)trials;
		fprintf(out, "whole,%g,%g,%g,%g,%u,%lld,%lld,%lld,%lld,%.6f,%.4f,%.4f,%.4f,%.4f\n",
			p_sub, p_ins, p_del, p_drop, trials,
			success, wrong, failed, dropped, (double)success / dt,
			(double)sum_repaired / dt, (double)sum_struct / dt,
			(double)sum_dead / dt, (double)sum_anom / dt);
		if (out != stdout) fclose(out);
	}

	for (int h = 0; h < 2; h++) vivi_dealloc(pristine[h]);
	organism_free(org);
	vivi_dealloc(payload);
	dna_free_caches();
	return 0;
}
