#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vivi/channel.h"
#include "vivi/dna.h"

typedef struct {
	double q_mean, q_sd, q_slope, q_min, q_max;
	double q_bias;
	uint32_t seed;
} calib_model;

static double u01(vivi_prng *p)
{
	return (double)vivi_prng_next64(p) / 18446744073709551616.0;
}

static double gauss(vivi_prng *p)
{
	double u1 = u01(p), u2 = u01(p);
	if (u1 < 1e-12) u1 = 1e-12;
	return sqrt(-2.0 * log(u1)) * cos(6.28318530717958647692 * u2);
}

static void set_digit(uint8_t *strand, size_t i, int d)
{
	size_t bytei = i / 4;
	int sh = 2 * (3 - (int)(i % 4));
	strand[bytei] = (uint8_t)((strand[bytei] & ~(3u << sh)) | ((uint32_t)d << sh));
}

static double q_true_at(const calib_model *m, size_t i, size_t n, vivi_prng *r)
{
	double frac = n ? (double)i / (double)n : 0.0;
	double q = m->q_mean - m->q_slope * frac + m->q_sd * gauss(r);
	if (q < m->q_min) q = m->q_min;
	if (q > m->q_max) q = m->q_max;
	return q;
}

static int calib_read(vivi_read *out, const uint8_t *strand, size_t slen, size_t bases,
	const calib_model *m, vivi_prng *r, const char **err)
{
	if (!vivi_read_from_bytes(out, strand, slen, bases, err)) return 0;
	if (!vivi_read_alloc_quality(out, 0, err)) { vivi_read_free(out); return 0; }
	for (size_t i = 0; i < bases; i++) {
		double qt = q_true_at(m, i, bases, r);
		double perr = pow(10.0, -qt / 10.0);
		if (u01(r) < perr) {
			int d = vivi_read_base(out, i);
			d = (d + 1 + (int)(vivi_prng_next(r) % 3u)) & 3;
			set_digit(out->strand.data, i, d);
		}
		double qrep = qt + m->q_bias;
		if (qrep < 0.0) qrep = 0.0;
		if (qrep > 255.0) qrep = 255.0;
		(void)vivi_read_set_quality(out, i, (uint8_t)(qrep + 0.5));
	}
	return 1;
}

static void usage(void)
{
	printf(
		"vivi_calib -- reliability calibration for a nominal short-read model\n"
		"usage: vivi_calib [options]\n"
		"  --size N       payload bytes (default 64)\n"
		"  --trials T     number of trials (default 1000)\n"
		"  --coverage C   reads voted per trial (default 5)\n"
		"  --q-mean M     mean Phred quality (default 36)\n"
		"  --q-sd S       quality std-dev (default 3)\n"
		"  --q-slope K    quality drop across the read (default 8)\n"
		"  --q-bias B     reported-quality offset, i.e. miscalibration (default 0)\n"
		"  --seed S       rng seed (default 1)\n"
		"  --header       print the CSV header\n");
}

int main(int argc, char **argv)
{
	calib_model m = { 36.0, 3.0, 8.0, 2.0, 45.0, 0.0, 1 };
	uint32_t size = 64, trials = 1000, coverage = 5;
	int header = 0;
	for (int i = 1; i < argc; i++) {
		const char *a = argv[i];
		const char *v = (i + 1 < argc) ? argv[i + 1] : nullptr;
		if (!strcmp(a, "--header")) { header = 1; continue; }
		if (!v) { usage(); return 2; }
		if (!strcmp(a, "--size")) { size = (uint32_t)strtoul(v, nullptr, 0); i++; }
		else if (!strcmp(a, "--trials")) { trials = (uint32_t)strtoul(v, nullptr, 0); i++; }
		else if (!strcmp(a, "--coverage")) { coverage = (uint32_t)strtoul(v, nullptr, 0); i++; }
		else if (!strcmp(a, "--q-mean")) { m.q_mean = strtod(v, nullptr); i++; }
		else if (!strcmp(a, "--q-sd")) { m.q_sd = strtod(v, nullptr); i++; }
		else if (!strcmp(a, "--q-slope")) { m.q_slope = strtod(v, nullptr); i++; }
		else if (!strcmp(a, "--q-bias")) { m.q_bias = strtod(v, nullptr); i++; }
		else if (!strcmp(a, "--seed")) { m.seed = (uint32_t)strtoul(v, nullptr, 0); i++; }
		else { usage(); return 2; }
	}
	if (size == 0 || trials == 0 || coverage == 0) { usage(); return 2; }

	const char *err = nullptr;
	uint8_t *payload = vivi_alloc(size);
	if (!payload) { fprintf(stderr, "vivi_calib: out of memory\n"); return 1; }
	vivi_prng r;
	vivi_prng_init(&r, m.seed);
	for (uint32_t i = 0; i < size; i++) payload[i] = (uint8_t)(vivi_prng_next(&r) & 255u);
	vivi_bytes strand = { 0 };
	dna_opts d = { 6, 0.05 };
	if (!dna_encode(&strand, payload, size, &d, &err)) {
		fprintf(stderr, "vivi_calib: encode: %s\n", err ? err : "?");
		return 1;
	}
	size_t bases = strand.len * 4;

	/* reliability of the reported quality, before and after a bias fit */
	long *n = calloc(256, sizeof(long));
	long *e = calloc(256, sizeof(long));
	long *nc = calloc(256, sizeof(long));
	long *ec = calloc(256, sizeof(long));
	if (!n || !e || !nc || !ec) { fprintf(stderr, "vivi_calib: out of memory\n"); return 1; }

	long long hard_bad = 0, soft_bad = 0;
	vivi_read *reads = calloc(coverage, sizeof(vivi_read));
	if (!reads) { fprintf(stderr, "vivi_calib: out of memory\n"); return 1; }
	const vivi_read **rp = calloc(coverage, sizeof(vivi_read *));
	if (!rp) { fprintf(stderr, "vivi_calib: out of memory\n"); return 1; }

	for (uint32_t t = 0; t < trials; t++) {
		for (uint32_t c = 0; c < coverage; c++) {
			vivi_prng rr;
			vivi_prng_init(&rr, ((uint64_t)m.seed << 32) ^ ((uint64_t)t * 1000003u + c + 1u));
			if (!calib_read(&reads[c], strand.data, strand.len, bases, &m, &rr, &err)) {
				fprintf(stderr, "vivi_calib: read: %s\n", err ? err : "?");
				return 1;
			}
			rp[c] = &reads[c];
			for (size_t i = 0; i < bases; i++) {
				int q = reads[c].qual[i];
				int was = (int)((strand.data[i / 4] >> (2 * (3 - i % 4))) & 3u);
				int got = vivi_read_base(&reads[c], i);
				n[q]++;
				if (got != was) e[q]++;
				double qc = q - m.q_bias;
				if (qc < 0.0) qc = 0.0;
				int qci = (int)(qc + 0.5);
				nc[qci]++;
				if (got != was) ec[qci]++;
			}
		}
		vivi_read hard, soft;
		if (!vivi_consensus_vote(&hard, rp, coverage, 0, &err)
			|| !vivi_consensus_vote(&soft, rp, coverage, 1, &err)) {
			fprintf(stderr, "vivi_calib: vote: %s\n", err ? err : "?");
			return 1;
		}
		for (size_t i = 0; i < bases; i++) {
			int was = (int)((strand.data[i / 4] >> (2 * (3 - i % 4))) & 3u);
			if (vivi_read_base(&hard, i) != was) hard_bad++;
			if (vivi_read_base(&soft, i) != was) soft_bad++;
		}
		vivi_read_free(&hard);
		vivi_read_free(&soft);
		for (uint32_t c = 0; c < coverage; c++) vivi_read_free(&reads[c]);
	}

	long N = 0, Nc = 0;
	double ece = 0.0, ece_c = 0.0;
	for (int q = 0; q < 256; q++) { N += n[q]; Nc += nc[q]; }
	for (int q = 1; q < 256; q++) {
		if (n[q]) {
			double pred = pow(10.0, -q / 10.0);
			double emp = (double)e[q] / (double)n[q];
			if (fabs(emp - pred) > 0.0) ece += (double)n[q] / (double)N * fabs(emp - pred);
		}
		if (nc[q]) {
			double pred = pow(10.0, -q / 10.0);
			double emp = (double)ec[q] / (double)nc[q];
			ece_c += (double)nc[q] / (double)Nc * fabs(emp - pred);
		}
	}

	if (header)
		printf("kind,key,value\n");
	printf("calib,params,size=%u trials=%u coverage=%u q_mean=%g q_sd=%g q_slope=%g q_bias=%g seed=%u\n",
		size, trials, coverage, m.q_mean, m.q_sd, m.q_slope, m.q_bias, m.seed);
	for (int q = 1; q < 256; q++) {
		if (!n[q]) continue;
		double pred = pow(10.0, -q / 10.0);
		double emp = (double)e[q] / (double)n[q];
		printf("reliability,%d,%.6g,%.6g,%ld\n", q, pred, emp, n[q]);
	}
	printf("ece,reported,%.6f,%ld\n", ece, N);
	printf("ece,corrected,%.6f,%ld\n", ece_c, Nc);
	printf("consensus,hard_bad,%lld,%zu\n", hard_bad, bases);
	printf("consensus,soft_bad,%lld,%zu\n", soft_bad, bases);

	free(n); free(e); free(nc); free(ec);
	free(reads); free(rp);
	vivi_bytes_free(&strand);
	vivi_dealloc(payload);
	return 0;
}
