#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vivi/vivi.h"
#include "vivi/genome.h"
#include "vivi/channel.h"
#include "vivi/parity.h"

enum {
	READ_ERR = -1,
	READ_NONE = 0,
	READ_OK = 1,
	READ_DROPPED = 2
};

typedef struct {
	const uint8_t *payload;
	size_t plen;
	int gene_raw;
	int h;
	uint32_t trials;
	uint32_t seed;
	double p_sub;
	double p_ins;
	double p_del;
	size_t n;
} bench_ctx;

typedef struct {
	long long success;
	long long wrong;
	long long failed;
	long long dropped;
	double storage_x;
	uint32_t trials_run;
	int skipped;
} bench_result;

static void usage(void)
{
	fprintf(stderr,
		"vivi_bench -- recovery strategies at equal storage redundancy\n"
		"usage: vivi_bench (--size N | --in FILE) [options]\n"
		"  --size N       synthetic payload of N random bytes\n"
		"  --in FILE      payload file (max 16 MiB)\n"
		"  --trials T     trials per strategy (default 1000)\n"
		"  --p-sub P      substitution probability per base (default 0)\n"
		"  --p-ins P      insertion probability per base (default 0)\n"
		"  --p-del P      deletion probability per base (default 0)\n"
		"  --seed S       base seed (default 1)\n"
		"  --gene-raw R   raw bytes per gene (default 64)\n"
		"  --header       print the CSV header first\n"
		"  --out FILE     write the CSV to FILE instead of stdout\n"
		"strategies: plain, replica 2/3, parity 2/4/8, inner 2/4/8, coverage 3/5\n"
		"columns: strategy,param,gene_raw,genes,storage_x,trials,success,wrong,failed,dropped,rate\n");
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

static size_t gene_offset(const bench_ctx *c, size_t i)
{
	return i * (size_t)c->gene_raw;
}

static size_t gene_take(const bench_ctx *c, size_t i)
{
	size_t off = gene_offset(c, i);
	if (off >= c->plen) return 0;
	size_t left = c->plen - off;
	return left < (size_t)c->gene_raw ? left : (size_t)c->gene_raw;
}

static size_t max_take(const bench_ctx *c)
{
	return c->plen < (size_t)c->gene_raw ? c->plen : (size_t)c->gene_raw;
}

static int read_gene(genome_gene *g, const vivi_bytes *mol, size_t id,
	const bench_ctx *c, uint32_t t, uint32_t r, uint32_t coverage, const char **err)
{
	vivi_channel_opts ch = { c->p_sub, c->p_ins, c->p_del, 0.0,
		c->seed + t * 977u + (uint32_t)id * 31u + r + 1u,
		0.0, 0.0, 0.0, 0.0, 0u, 0.0 };
	vivi_read rd;
	if (coverage > 1) {
		vivi_consensus_opts cn = { ch, coverage, 0 };
		if (!vivi_consensus_read(&rd, mol->data, mol->len, &cn, err)) return READ_ERR;
	} else {
		if (!vivi_channel_read(&rd, mol->data, mol->len, &ch, err)) return READ_ERR;
	}
	if (rd.dropped) {
		vivi_read_free(&rd);
		return READ_DROPPED;
	}
	int ok = genome_gene_read(g, rd.strand.data, rd.strand.len, (int)id, err)
		? READ_OK : READ_NONE;
	vivi_read_free(&rd);
	return ok;
}

static bool build_plain(vivi_bytes *mol, const bench_ctx *c, const char **err)
{
	for (size_t i = 0; i < c->n; i++) {
		size_t off = gene_offset(c, i);
		if (!genome_gene_encode(&mol[i], (int)i, 0, 0, c->h,
			c->payload + off, gene_take(c, i), err))
			return false;
	}
	return true;
}

static bool build_inner(vivi_bytes *mol, const bench_ctx *c, int inner_m, const char **err)
{
	for (size_t i = 0; i < c->n; i++) {
		size_t off = gene_offset(c, i);
		if (!genome_gene_encode_inner(&mol[i], (int)i, 0, 0, c->h, inner_m,
			c->payload + off, gene_take(c, i), err))
			return false;
	}
	return true;
}

static void free_gene_data(uint8_t **datas, size_t n)
{
	if (!datas) return;
	for (size_t i = 0; i < n; i++) vivi_dealloc(datas[i]);
	vivi_dealloc(datas);
}

static void free_molecules(vivi_bytes *mol, size_t n)
{
	if (!mol) return;
	for (size_t i = 0; i < n; i++) vivi_dealloc(mol[i].data);
	vivi_dealloc(mol);
}

static bool run_gene_strategy(const bench_ctx *c, const vivi_bytes *mol,
	uint32_t copies, uint32_t coverage, bench_result *res, const char **err)
{
	memset(res, 0, sizeof(*res));
	double total = 0.0;
	for (size_t i = 0; i < c->n; i++)
		total += (double)mol[i].len;
	res->storage_x = total * (double)copies / (double)c->plen;
	res->trials_run = c->trials;

	uint8_t **datas = vivi_zalloc(c->n, sizeof(uint8_t *));
	size_t *dlens = vivi_alloc((c->n ? c->n : 1) * sizeof(size_t));
	if (!datas || !dlens) {
		free_gene_data(datas, c->n);
		vivi_dealloc(dlens);
		*err = "out of memory";
		return false;
	}

	const uint8_t *const base = c->payload;
	for (uint32_t t = 0; t < c->trials; t++) {
		int all = 1;
		for (size_t i = 0; i < c->n; i++) {
			datas[i] = nullptr;
			dlens[i] = 0;
		}
		for (size_t i = 0; i < c->n; i++) {
			int have = 0;
			genome_gene g;
			if (coverage > 1) {
				int rc = read_gene(&g, &mol[i], i, c, t, 0, coverage, err);
				if (rc == READ_ERR) {
					free_gene_data(datas, c->n);
					vivi_dealloc(dlens);
					return false;
				}
				if (rc == READ_DROPPED) {
					res->dropped++;
				} else if (rc == READ_OK) {
					datas[i] = g.data;
					dlens[i] = g.data_len;
					have = 1;
				}
			} else {
				for (uint32_t r = 0; r < copies; r++) {
					int rc = read_gene(&g, &mol[i], i, c, t, r, 1, err);
					if (rc == READ_ERR) {
						free_gene_data(datas, c->n);
						vivi_dealloc(dlens);
						return false;
					}
					if (rc == READ_DROPPED) {
						res->dropped++;
						continue;
					}
					if (rc == READ_OK) {
						datas[i] = g.data;
						dlens[i] = g.data_len;
						have = 1;
						break;
					}
				}
			}
			if (!have) all = 0;
		}
		if (!all) {
			res->failed++;
		} else {
			int match = 1;
			for (size_t i = 0; i < c->n; i++) {
				size_t take = gene_take(c, i);
				if (dlens[i] != take
					|| (take && memcmp(datas[i], base + gene_offset(c, i), take) != 0))
					match = 0;
			}
			if (match) res->success++;
			else res->wrong++;
		}
		for (size_t i = 0; i < c->n; i++) {
			vivi_dealloc(datas[i]);
			datas[i] = nullptr;
		}
	}
	free_gene_data(datas, c->n);
	vivi_dealloc(dlens);
	return true;
}

static bool run_parity(const bench_ctx *c, uint32_t m, bench_result *res, const char **err)
{
	memset(res, 0, sizeof(*res));
	size_t n = c->n;
	size_t shard_len = (size_t)c->gene_raw;
	if (n + (size_t)m > 255) {
		res->skipped = 1;
		return true;
	}

	uint8_t **data = vivi_zalloc(n, sizeof(uint8_t *));
	uint8_t **shards = vivi_zalloc(n + m, sizeof(uint8_t *));
	vivi_bytes *mol = vivi_zalloc(n + m, sizeof(vivi_bytes));
	size_t *mlen = vivi_alloc((n + m) * sizeof(size_t));
	uint8_t **recv = vivi_zalloc(n + m, sizeof(uint8_t *));
	uint8_t *present = vivi_zalloc(n + m, 1);
	uint8_t **out = vivi_zalloc(n, sizeof(uint8_t *));
	if (!data || !shards || !mol || !mlen || !recv || !present || !out) {
		*err = "out of memory";
		goto fail;
	}
	for (size_t i = 0; i < n; i++) {
		data[i] = vivi_zalloc(shard_len, 1);
		if (!data[i]) {
			*err = "out of memory";
			goto fail;
		}
		memcpy(data[i], c->payload + gene_offset(c, i), gene_take(c, i));
	}
	if (!vivi_parity_encode(shards, (const uint8_t *const *)data, n, m, shard_len, err))
		goto fail;
	for (size_t s = 0; s < n + m; s++) {
		size_t rawlen = (s < n) ? gene_take(c, s) : shard_len;
		if (!genome_gene_encode(&mol[s], (int)s, 0, 0, c->h, shards[s], rawlen, err))
			goto fail;
		mlen[s] = mol[s].len;
	}
	vivi_parity_release(shards, n + m);
	vivi_dealloc(shards);
	shards = nullptr;
	for (size_t i = 0; i < n; i++) vivi_dealloc(data[i]);
	vivi_dealloc(data);
	data = nullptr;

	double total = 0.0;
	for (size_t s = 0; s < n + m; s++) total += (double)mlen[s];
	res->storage_x = total / (double)c->plen;
	res->trials_run = c->trials;

	for (uint32_t t = 0; t < c->trials; t++) {
		for (size_t s = 0; s < n + m; s++) present[s] = 0;
		size_t present_count = 0;
		for (size_t s = 0; s < n + m; s++) {
			genome_gene g;
			int rc = read_gene(&g, &mol[s], s, c, t, 0, 1, err);
			if (rc == READ_ERR) goto fail;
			if (rc == READ_DROPPED) {
				res->dropped++;
				continue;
			}
			if (rc == READ_OK) {
				recv[s] = vivi_zalloc(shard_len, 1);
				if (!recv[s]) {
					vivi_dealloc(g.data);
					*err = "out of memory";
					goto fail;
				}
				memcpy(recv[s], g.data, g.data_len);
				present[s] = 1;
				present_count++;
				vivi_dealloc(g.data);
			}
		}
		if (present_count < n) {
			res->failed++;
		} else if (!vivi_parity_decode(out, (const uint8_t *const *)recv, present,
			n, m, shard_len, err)) {
			res->failed++;
		} else {
			int match = 1;
			size_t off = 0;
			for (size_t i = 0; i < n; i++) {
				size_t take = gene_take(c, i);
				if (memcmp(out[i], c->payload + off, take) != 0) match = 0;
				off += take;
			}
			if (match) res->success++;
			else res->wrong++;
			vivi_parity_release(out, n);
		}
		for (size_t s = 0; s < n + m; s++) {
			vivi_dealloc(recv[s]);
			recv[s] = nullptr;
		}
	}

	free_molecules(mol, n + m);
	vivi_dealloc(mlen);
	vivi_dealloc(recv);
	vivi_dealloc(present);
	vivi_dealloc(out);
	return true;

fail:
	if (data) {
		for (size_t i = 0; i < n; i++) vivi_dealloc(data[i]);
		vivi_dealloc(data);
	}
	if (shards) {
		vivi_parity_release(shards, n + m);
		vivi_dealloc(shards);
	}
	free_molecules(mol, n + m);
	if (recv) {
		for (size_t s = 0; s < n + m; s++) vivi_dealloc(recv[s]);
		vivi_dealloc(recv);
	}
	if (out) {
		vivi_parity_release(out, n);
		vivi_dealloc(out);
	}
	vivi_dealloc(mlen);
	vivi_dealloc(present);
	return false;
}

static void emit_row(FILE *f, const char *strategy, uint32_t param, const bench_ctx *c,
	const bench_result *r)
{
	double rate = r->trials_run ? (double)r->success / (double)r->trials_run : 0.0;
	fprintf(f, "%s,%u,%d,%zu,%.3f,%u,%lld,%lld,%lld,%lld,%.6f\n",
		strategy, param, c->gene_raw, c->n, r->storage_x, r->trials_run,
		r->success, r->wrong, r->failed, r->dropped, rate);
	if (r->skipped)
		fprintf(stderr, "bench: %s param=%u skipped (not applicable)\n", strategy, param);
	else
		fprintf(stderr, "bench: %s param=%u storage_x=%.3f success=%lld rate=%.6f\n",
			strategy, param, r->storage_x, r->success, rate);
}

static const uint32_t REPLICA_KS[2] = { 2, 3 };
static const uint32_t PARITY_MS[3] = { 2, 4, 8 };
static const uint32_t INNER_MS[3] = { 2, 4, 8 };
static const uint32_t COVERAGE_KS[2] = { 3, 5 };

int main(int argc, char **argv)
{
	const char *in_path = nullptr, *out_path = nullptr;
	uint32_t size = 0, trials = 1000, seed = 1;
	int have_size = 0, header = 0;
	double p_sub = 0.0, p_ins = 0.0, p_del = 0.0;
	bench_ctx ctx = { 0 };
	ctx.gene_raw = 64;
	ctx.h = 3;
	uint8_t *payload = nullptr;
	FILE *fout = stdout;
	vivi_bytes *plain = nullptr;
	bench_result res, cov[2];
	const char *err = nullptr;
	double plain_storage = 0.0;
	int rc = 0;

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
		else if (!strcmp(a, "--gene-raw")) {
			uint32_t g;
			if (!parse_u32(v, &g) || g < 16 || g > 65535) { usage(); return 2; }
			ctx.gene_raw = (int)g;
			i++;
		}
		else if (!strcmp(a, "--header")) { header = 1; }
		else { usage(); return 2; }
	}
	if ((!in_path && !have_size) || (in_path && have_size)) { usage(); return 2; }
	if (trials == 0) { usage(); return 2; }

	size_t plen = 0;
	if (in_path) {
		FILE *f = fopen(in_path, "rb");
		if (!f) {
			fprintf(stderr, "vivi_bench: cannot open %s\n", in_path);
			return 1;
		}
		fseek(f, 0, SEEK_END);
		long sz = ftell(f);
		fseek(f, 0, SEEK_SET);
		if (sz < 0 || (unsigned long)sz > 16u * 1024u * 1024u) {
			fprintf(stderr, "vivi_bench: file too large (max 16 MiB)\n");
			fclose(f);
			return 1;
		}
		plen = (size_t)sz;
		payload = vivi_alloc(plen ? plen : 1);
		if (!payload || (plen && fread(payload, 1, plen, f) != plen)) {
			fprintf(stderr, "vivi_bench: read error\n");
			vivi_dealloc(payload);
			fclose(f);
			return 1;
		}
		fclose(f);
	} else {
		plen = size;
		payload = vivi_alloc(plen ? plen : 1);
		if (!payload) {
			fprintf(stderr, "vivi_bench: out of memory\n");
			return 1;
		}
		vivi_prng r;
		vivi_prng_init(&r, seed);
		for (size_t i = 0; i < plen; i++) payload[i] = (uint8_t)(vivi_prng_next(&r) & 255u);
	}
	if (plen == 0) {
		fprintf(stderr, "vivi_bench: empty payload\n");
		vivi_dealloc(payload);
		return 1;
	}

	ctx.payload = payload;
	ctx.plen = plen;
	ctx.trials = trials;
	ctx.seed = seed;
	ctx.p_sub = p_sub;
	ctx.p_ins = p_ins;
	ctx.p_del = p_del;
	ctx.n = (plen + (size_t)ctx.gene_raw - 1) / (size_t)ctx.gene_raw;
	if (ctx.n > 255) {
		fprintf(stderr, "vivi_bench: too many genes (%zu > 255); use a larger --gene-raw\n", ctx.n);
		vivi_dealloc(payload);
		return 1;
	}

	if (out_path && !(fout = fopen(out_path, "wb"))) {
		fprintf(stderr, "vivi_bench: cannot write %s\n", out_path);
		vivi_dealloc(payload);
		return 1;
	}
	if (header)
		fprintf(fout, "strategy,param,gene_raw,genes,storage_x,trials,success,wrong,failed,dropped,rate\n");

	plain = vivi_zalloc(ctx.n, sizeof(vivi_bytes));
	if (!plain || !build_plain(plain, &ctx, &err)) {
		fprintf(stderr, "vivi_bench: encode: %s\n", err ? err : "?");
		rc = 1;
		goto cleanup;
	}
	if (!run_gene_strategy(&ctx, plain, 1, 1, &res, &err)) {
		fprintf(stderr, "vivi_bench: plain: %s\n", err ? err : "?");
		rc = 1;
		goto cleanup;
	}
	emit_row(fout, "plain", 0, &ctx, &res);
	plain_storage = res.storage_x;

	for (int i = 0; i < 2; i++) {
		if (!run_gene_strategy(&ctx, plain, REPLICA_KS[i], 1, &res, &err)) {
			fprintf(stderr, "vivi_bench: replica: %s\n", err ? err : "?");
			rc = 1;
			goto cleanup;
		}
		emit_row(fout, "replica", REPLICA_KS[i], &ctx, &res);
	}

	for (int i = 0; i < 2; i++) {
		if (p_ins != 0.0 || p_del != 0.0) {
			memset(&cov[i], 0, sizeof(cov[i]));
			cov[i].skipped = 1;
			cov[i].storage_x = plain_storage;
		} else if (!run_gene_strategy(&ctx, plain, 1, COVERAGE_KS[i], &cov[i], &err)) {
			fprintf(stderr, "vivi_bench: coverage: %s\n", err ? err : "?");
			rc = 1;
			goto cleanup;
		}
	}

	free_molecules(plain, ctx.n);
	plain = nullptr;

	for (int i = 0; i < 3; i++) {
		if (!run_parity(&ctx, PARITY_MS[i], &res, &err)) {
			fprintf(stderr, "vivi_bench: parity: %s\n", err ? err : "?");
			rc = 1;
			goto cleanup;
		}
		emit_row(fout, "parity", PARITY_MS[i], &ctx, &res);
	}

	for (int k = 0; k < 3; k++) {
		if (max_take(&ctx) + INNER_MS[k] > 255) {
			memset(&res, 0, sizeof(res));
			res.skipped = 1;
			emit_row(fout, "inner", INNER_MS[k], &ctx, &res);
			continue;
		}
		vivi_bytes *inner = vivi_zalloc(ctx.n, sizeof(vivi_bytes));
		if (!inner || !build_inner(inner, &ctx, (int)INNER_MS[k], &err)) {
			fprintf(stderr, "vivi_bench: encode: %s\n", err ? err : "?");
			free_molecules(inner, ctx.n);
			rc = 1;
			goto cleanup;
		}
		if (!run_gene_strategy(&ctx, inner, 1, 1, &res, &err)) {
			fprintf(stderr, "vivi_bench: inner: %s\n", err ? err : "?");
			free_molecules(inner, ctx.n);
			rc = 1;
			goto cleanup;
		}
		emit_row(fout, "inner", INNER_MS[k], &ctx, &res);
		free_molecules(inner, ctx.n);
	}

	for (int i = 0; i < 2; i++)
		emit_row(fout, "coverage", COVERAGE_KS[i], &ctx, &cov[i]);

cleanup:
	free_molecules(plain, ctx.n);
	if (fout != stdout) fclose(fout);
	vivi_dealloc(payload);
	dna_free_caches();
	return rc;
}
