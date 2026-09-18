#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vivi/organism.h"

typedef struct {
	int gene_raw, parity, primer, codon, h, units, max_gen;
	int damage;
	uint32_t seed;
} cli_opts;

static const char *prog = "vivi_cli";

static void usage(void)
{
	printf(
		"vivi_cli -- file front-end for the Vivian C framework\n"
		"usage: vivi_cli <command> <input> [output] [options]\n"
		"  encode  <in> <out.vivi>  wrap a file as a VIV14NB4NSH33 organism\n"
		"  decode  <in.vivi> <out>  read the payload (repairing with the homolog)\n"
		"  repair  <in.vivi> <out>  checkpoint and re-serialize the container\n"
		"  inspect <in.vivi>        print container and chromosome structure\n"
		"  verify  <in.vivi>        check integrity without touching the data\n"
		"options:\n"
		"  --gene-raw N   raw bytes per gene (16..65535, default 1024)\n"
		"  --parity M     parity genes per chromosome (0..16, default 0)\n"
		"  --primer B     on-strand barcode (0..65535, default 0 = off)\n"
		"  --codon        codon payload instead of dense\n"
		"  --h H          dense homopolymer limit (3..12, default 3)\n"
		"  --units U      telomere repeat units (>= 1, default 4)\n"
		"  --max-gen G    Hayflick limit (default 60)\n"
		"  --damage N     decode: flip N bases before reading\n"
		"  --seed S       damage seed (default 1)\n");
}

static int parse_int(const char *s, long lo, long hi, long *out)
{
	char *end = nullptr;
	if (!s || !*s) return 0;
	long v = strtol(s, &end, 0);
	if (*end || v < lo || v > hi) return 0;
	*out = v;
	return 1;
}

static uint8_t *read_file(const char *path, size_t *out_len, const char **err)
{
	FILE *f = fopen(path, "rb");
	if (!f) { *err = "cannot open input file"; return nullptr; }
	if (fseek(f, 0, SEEK_END) != 0) { fclose(f); *err = "cannot seek input file"; return nullptr; }
	long n = ftell(f);
	if (n < 0) { fclose(f); *err = "cannot size input file"; return nullptr; }
	rewind(f);
	size_t len = (size_t)n;
	uint8_t *buf = malloc(len ? len : 1);
	if (!buf) { fclose(f); *err = "out of memory"; return nullptr; }
	if (len && fread(buf, 1, len, f) != len) {
		free(buf);
		fclose(f);
		*err = "cannot read input file";
		return nullptr;
	}
	fclose(f);
	*out_len = len;
	return buf;
}

static int write_file(const char *path, const uint8_t *data, size_t len, const char **err)
{
	FILE *f = fopen(path, "wb");
	if (!f) { *err = "cannot open output file"; return 0; }
	if (len && fwrite(data, 1, len, f) != len) {
		fclose(f);
		*err = "cannot write output file";
		return 0;
	}
	if (fclose(f) != 0) { *err = "cannot close output file"; return 0; }
	return 1;
}

static int build_organism(vivi_organism **out, const uint8_t *payload, size_t plen,
	const cli_opts *co, const char **err)
{
	if (plen == 0) { *err = "empty input has no genes"; return 0; }
	int genes_per = 256;
	if (co->parity > 0) genes_per = 255 - co->parity;
	size_t cap = (size_t)genes_per * (size_t)co->gene_raw;
	size_t nchr = (plen + cap - 1) / cap;
	if (nchr == 0) nchr = 1;
	if (nchr > 255) { *err = "payload too large for this layout"; return 0; }
	int *ids = malloc(nchr * sizeof(int));
	const uint8_t **datas = malloc(nchr * sizeof(uint8_t *));
	size_t *lens = malloc(nchr * sizeof(size_t));
	chr_opts *opts = malloc(nchr * sizeof(chr_opts));
	if (!ids || !datas || !lens || !opts) {
		free(ids); free(datas); free(lens); free(opts);
		*err = "out of memory";
		return 0;
	}
	for (size_t i = 0; i < nchr; i++) {
		size_t off = i * cap;
		size_t take = (plen - off < cap) ? (plen - off) : cap;
		ids[i] = (int)i;
		datas[i] = payload + off;
		lens[i] = take;
		memset(&opts[i], 0, sizeof(chr_opts));
		opts[i].gene_raw = co->gene_raw;
		opts[i].codon = co->codon;
		opts[i].h = co->h;
		opts[i].units = co->units;
		opts[i].parity = co->parity;
		opts[i].primer = co->primer;
	}
	int ok = organism_new(out, ids, opts, datas, lens, nchr, co->max_gen, err);
	free(ids); free(datas); free(lens); free(opts);
	return ok;
}

static void print_report(const cell_report *rep)
{
	printf("repaired=%d structural=%d dead=%d anomaly=%d renewed=%d failed=%d\n",
		rep->repaired, rep->structural, rep->dead, rep->anomaly,
		rep->renewed, rep->failed);
}

static size_t max_data_genes(const chr_record *a, const chr_record *b)
{
	size_t ga = a->gene_count > (size_t)a->parity ? a->gene_count - (size_t)a->parity : 0;
	size_t gb = b->gene_count > (size_t)b->parity ? b->gene_count - (size_t)b->parity : 0;
	return ga > gb ? ga : gb;
}

static int gene_intact(const chr_record *r, int id)
{
	for (size_t i = 0; i < r->gene_count; i++)
		if (r->genes[i].id == id && r->genes[i].crc_ok) return 1;
	return 0;
}

static int do_encode(const uint8_t *in, size_t in_len, const char *out_path,
	const cli_opts *co)
{
	const char *err = nullptr;
	vivi_organism *o = nullptr;
	if (!build_organism(&o, in, in_len, co, &err)) {
		fprintf(stderr, "%s: encode: %s\n", prog, err ? err : "?");
		return 1;
	}
	vivi_bytes ser = { 0 };
	if (!organism_serialize14nb(&ser, o, &err)) {
		fprintf(stderr, "%s: serialize: %s\n", prog, err ? err : "?");
		organism_free(o);
		return 1;
	}
	int ok = write_file(out_path, ser.data, ser.len, &err);
	if (!ok) fprintf(stderr, "%s: %s\n", prog, err ? err : "?");
	else printf("encoded %zu bytes -> %zu bytes across %zu chromosome(s)\n",
		in_len, ser.len, o->nchr);
	vivi_bytes_free(&ser);
	organism_free(o);
	return ok ? 0 : 1;
}

static int do_decode(const uint8_t *in, size_t in_len, const char *out_path,
	const cli_opts *co)
{
	const char *err = nullptr;
	vivi_organism *o = nullptr;
	if (!organism_deserialize(&o, in, in_len, &err)) {
		fprintf(stderr, "%s: load: %s\n", prog, err ? err : "?");
		return 1;
	}
	if (co->damage > 0 && !organism_damage(o, co->damage, co->seed, &err)) {
		fprintf(stderr, "%s: damage: %s\n", prog, err ? err : "?");
		organism_free(o);
		return 1;
	}
	if (co->damage > 0) printf("inflicted %d base change(s), seed %u\n", co->damage, co->seed);
	vivi_bytes *data = nullptr;
	cell_report rep;
	if (!organism_read(&data, o, &rep, &err)) {
		fprintf(stderr, "%s: read: %s\n", prog, err ? err : "?");
		organism_free(o);
		return 1;
	}
	print_report(&rep);
	size_t total = 0;
	for (size_t i = 0; i < o->nchr; i++) total += data[i].len;
	uint8_t *flat = malloc(total ? total : 1);
	if (!flat) {
		vivi_bytes_free_n(data, o->nchr);
		organism_free(o);
		fprintf(stderr, "%s: out of memory\n", prog);
		return 1;
	}
	size_t pos = 0;
	for (size_t i = 0; i < o->nchr; i++) {
		memcpy(flat + pos, data[i].data, data[i].len);
		pos += data[i].len;
	}
	int ok = write_file(out_path, flat, total, &err);
	if (!ok) fprintf(stderr, "%s: %s\n", prog, err ? err : "?");
	else printf("decoded %zu bytes from %zu chromosome(s)\n", total, o->nchr);
	free(flat);
	vivi_bytes_free_n(data, o->nchr);
	organism_free(o);
	return ok ? 0 : 1;
}

static int do_repair(const uint8_t *in, size_t in_len, const char *out_path)
{
	const char *err = nullptr;
	vivi_organism *o = nullptr;
	if (!organism_deserialize(&o, in, in_len, &err)) {
		fprintf(stderr, "%s: load: %s\n", prog, err ? err : "?");
		return 1;
	}
	cell_report rep;
	if (!organism_checkpoint_checked(o, &rep, &err)) {
		fprintf(stderr, "%s: checkpoint: %s\n", prog, err ? err : "?");
		organism_free(o);
		return 1;
	}
	print_report(&rep);
	vivi_bytes ser = { 0 };
	if (!organism_serialize14nb(&ser, o, &err)) {
		fprintf(stderr, "%s: serialize: %s\n", prog, err ? err : "?");
		organism_free(o);
		return 1;
	}
	int ok = write_file(out_path, ser.data, ser.len, &err);
	if (!ok) fprintf(stderr, "%s: %s\n", prog, err ? err : "?");
	else printf("wrote repaired container (%zu bytes)\n", ser.len);
	vivi_bytes_free(&ser);
	organism_free(o);
	return ok ? 0 : 1;
}

static int do_inspect(const uint8_t *in, size_t in_len)
{
	const char *err = nullptr;
	vivi_organism *o = nullptr;
	if (!organism_deserialize(&o, in, in_len, &err)) {
		fprintf(stderr, "%s: load: %s\n", prog, err ? err : "?");
		return 1;
	}
	int rev = organism_container_version(in, in_len);
	printf("container: %s (revision %d)\n",
		rev == 143 ? "VIV14NB4NSH33" : rev == 141 ? "VIV14N" :
		rev == 14 ? "VIV14" : rev == 1 ? "VIV1" : "unknown", rev);
	printf("chromosomes: %zu   generation: %d   max_gen: %d   dead: %d\n",
		o->nchr, o->generation, o->max_gen, o->dead);
	for (size_t i = 0; i < o->nchr; i++) {
		chr_record rec;
		const char *e = nullptr;
		if (chr_parse(&rec, o->hom[0][i], o->hlen[0][i], -1, &e)) {
			int gene_raw = o->chr_opts[i].gene_raw ? o->chr_opts[i].gene_raw : 1024;
			char rawbuf[32];
			if (rec.cen_version == 2) snprintf(rawbuf, sizeof(rawbuf), "%zu", rec.rawlen);
			else snprintf(rawbuf, sizeof(rawbuf), "n/a");
			printf("  chr %d: genes=%zu parity=%d rawlen=%s gen=%d cen=%d tel=%d primer=%d gene_raw=%d\n",
				o->chr_ids[i], rec.gene_count, rec.parity, rawbuf, rec.generation,
				rec.cen_ok, rec.telomere_ok, rec.primer, gene_raw);
			chr_record_free(&rec);
		} else {
			printf("  chr %d: unreadable (%s)\n", o->chr_ids[i], e ? e : "?");
		}
	}
	organism_free(o);
	return 0;
}

static int do_verify(const uint8_t *in, size_t in_len)
{
	const char *err = nullptr;
	vivi_organism *o = nullptr;
	if (!organism_deserialize(&o, in, in_len, &err)) {
		fprintf(stderr, "%s: load: %s\n", prog, err ? err : "?");
		return 2;
	}
	int bad = 0;
	for (size_t i = 0; i < o->nchr; i++) {
		chr_record a, b;
		const char *ea = nullptr, *eb = nullptr;
		int pa = chr_parse(&a, o->hom[0][i], o->hlen[0][i], -1, &ea);
		int pb = chr_parse(&b, o->hom[1][i], o->hlen[1][i], -1, &eb);
		if (!pa && !pb) {
			printf("  chr %d: BROKEN (both homologs unreadable)\n", o->chr_ids[i]);
			bad = 1;
			continue;
		}
		size_t total = 0;
		if (pa && pb) total = max_data_genes(&a, &b);
		else total = pa ? max_data_genes(&a, &a) : max_data_genes(&b, &b);
		size_t clean = 0, recoverable = 0;
		for (size_t g = 0; g < total; g++) {
			int ia = pa && gene_intact(&a, (int)g);
			int ib = pb && gene_intact(&b, (int)g);
			if (ia && ib) clean++;
			if (ia || ib) recoverable++;
		}
		int cen = (pa && a.cen_ok) || (pb && b.cen_ok);
		if (recoverable == total && cen && clean == total) {
			printf("  chr %d: intact (%zu loci, centromere ok)\n", o->chr_ids[i], total);
		} else if (recoverable == total && cen) {
			printf("  chr %d: damaged but recoverable (clean %zu/%zu, centromere ok)\n",
				o->chr_ids[i], clean, total);
			bad = 1;
		} else {
			printf("  chr %d: BROKEN (recoverable %zu/%zu, centromere %s)\n",
				o->chr_ids[i], recoverable, total, cen ? "ok" : "BROKEN");
			bad = 1;
		}
		if (pa) chr_record_free(&a);
		if (pb) chr_record_free(&b);
	}
	organism_free(o);
	printf(bad ? "verify: DAMAGE DETECTED\n" : "verify: OK\n");
	return bad ? 1 : 0;
}

int main(int argc, char **argv)
{
	cli_opts co = { 1024, 0, 0, 0, 3, 4, 60, 0, 1 };
	if (argc < 3) { usage(); return 2; }
	const char *cmd = argv[1];
	const char *in_path = argv[2];
	const char *out_path = nullptr;
	int argi = 3;
	if (strcmp(cmd, "inspect") != 0 && strcmp(cmd, "verify") != 0) {
		if (argc < 4) { usage(); return 2; }
		out_path = argv[3];
		argi = 4;
	}
	for (int i = argi; i < argc; i++) {
		const char *a = argv[i];
		long v;
		if (strcmp(a, "--codon") == 0) { co.codon = 1; continue; }
		if (i + 1 >= argc) { fprintf(stderr, "%s: missing value for %s\n", prog, a); return 2; }
		const char *val = argv[++i];
		if (strcmp(a, "--gene-raw") == 0) {
			if (!parse_int(val, 16, 65535, &v)) { fprintf(stderr, "%s: bad --gene-raw\n", prog); return 2; }
			co.gene_raw = (int)v;
		} else if (strcmp(a, "--parity") == 0) {
			if (!parse_int(val, 0, 16, &v)) { fprintf(stderr, "%s: bad --parity\n", prog); return 2; }
			co.parity = (int)v;
		} else if (strcmp(a, "--primer") == 0) {
			if (!parse_int(val, 0, 65535, &v)) { fprintf(stderr, "%s: bad --primer\n", prog); return 2; }
			co.primer = (int)v;
		} else if (strcmp(a, "--h") == 0) {
			if (!parse_int(val, 3, 12, &v)) { fprintf(stderr, "%s: bad --h\n", prog); return 2; }
			co.h = (int)v;
		} else if (strcmp(a, "--units") == 0) {
			if (!parse_int(val, 1, 1000000, &v)) { fprintf(stderr, "%s: bad --units\n", prog); return 2; }
			co.units = (int)v;
		} else if (strcmp(a, "--max-gen") == 0) {
			if (!parse_int(val, 0, 65535, &v)) { fprintf(stderr, "%s: bad --max-gen\n", prog); return 2; }
			co.max_gen = (int)v;
		} else if (strcmp(a, "--damage") == 0) {
			if (!parse_int(val, 0, 1 << 20, &v)) { fprintf(stderr, "%s: bad --damage\n", prog); return 2; }
			co.damage = (int)v;
		} else if (strcmp(a, "--seed") == 0) {
			if (!parse_int(val, 0, 2147483647L, &v)) { fprintf(stderr, "%s: bad --seed\n", prog); return 2; }
			co.seed = (uint32_t)v;
		} else {
			fprintf(stderr, "%s: unknown option %s\n", prog, a);
			return 2;
		}
	}

	const char *err = nullptr;
	size_t in_len = 0;
	uint8_t *in = read_file(in_path, &in_len, &err);
	if (!in) { fprintf(stderr, "%s: %s\n", prog, err ? err : "?"); return 2; }

	int rc;
	if (strcmp(cmd, "encode") == 0) rc = do_encode(in, in_len, out_path, &co);
	else if (strcmp(cmd, "decode") == 0) rc = do_decode(in, in_len, out_path, &co);
	else if (strcmp(cmd, "repair") == 0) rc = do_repair(in, in_len, out_path);
	else if (strcmp(cmd, "inspect") == 0) rc = do_inspect(in, in_len);
	else if (strcmp(cmd, "verify") == 0) rc = do_verify(in, in_len);
	else { usage(); rc = 2; }
	free(in);
	return rc;
}
