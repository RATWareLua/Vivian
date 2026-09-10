/* xcheck.c -- deterministic C side of the C<->Lua byte-parity check.
 *
 * Prints one labeled hex line per scenario; framework/luau/xcheck.lua
 * prints the same lines and CI diffs the two outputs. Keep the scenario
 * list, payload seeds and all option values in sync with the Lua file. */
#include "vivi/vivi.h"
#include "vivi/genome.h"
#include "vivi/chromosome.h"
#include "vivi/cell.h"
#include "vivi/organism.h"
#include <stdio.h>
#include <string.h>

static void hexline(const char *label, const uint8_t *p, size_t n)
{
	printf("%s ", label);
	for (size_t i = 0; i < n; i++) printf("%02x", p[i]);
	printf("\n");
}

static void payload_seed(uint8_t *p, size_t n, uint32_t seed)
{
	uint32_t x = seed;
	for (size_t i = 0; i < n; i++) {
		x = x * 1664525u + 1013904223u;
		p[i] = (uint8_t)(x >> 24);
	}
}

static void scenario_chromosomes(void)
{
	uint8_t data[200];
	payload_seed(data, sizeof(data), 0xC0FFEE01u);
	hexline("payload200", data, sizeof(data));

	chr_opts o1 = { 64, 0, 3, 3, 9, 0 };   /* gene_raw, codon, h, units, flags, parity */
	vivi_bytes c1 = { 0 };
	if (!chr_encode(&c1, 7, data, sizeof(data), &o1, nullptr)) {
		printf("chr_dense fail\n");
		return;
	}
	hexline("chr_dense", c1.data, c1.len);

	chr_opts o2 = { 64, 1, 3, 3, 9, 0 };
	vivi_bytes c2 = { 0 };
	if (!chr_encode(&c2, 7, data, sizeof(data), &o2, nullptr)) {
		printf("chr_codon fail\n");
		return;
	}
	hexline("chr_codon", c2.data, c2.len);

	vivi_bytes c3 = { 0 };
	if (!chr_set_generation(&c3, c1.data, c1.len, 5, nullptr)) {
		printf("chr_gen5 fail\n");
		return;
	}
	hexline("chr_gen5", c3.data, c3.len);

	chr_opts o4 = { 64, 0, 3, 3, 9, 2 };
	vivi_bytes c4 = { 0 };
	if (!chr_encode(&c4, 7, data, sizeof(data), &o4, nullptr)) {
		printf("chr_par2 fail\n");
		return;
	}
	hexline("chr_par2", c4.data, c4.len);

	chr_opts o5 = { 48, 1, 3, 3, 9, 3 };
	vivi_bytes c5 = { 0 };
	if (!chr_encode(&c5, 8, data, sizeof(data), &o5, nullptr)) {
		printf("chr_par3_codon fail\n");
		return;
	}
	hexline("chr_par3_codon", c5.data, c5.len);

	chr_record r4;
	if (!chr_parse(&r4, c4.data, c4.len, -1, nullptr)) {
		printf("chr_par2_fields fail\n");
		return;
	}
	printf("chr_par2_fields %d %d %zu %d %d\n", r4.cen_version, r4.parity,
		r4.rawlen, r4.ngenes, (int)r4.gene_count);
	chr_record_free(&r4);

	/* damage two data genes in the parity chromosome, then recover */
	uint8_t *broken = vivi_alloc(c4.len);
	if (broken) {
		memcpy(broken, c4.data, c4.len);
		chr_record rb;
		(void)chr_parse(&rb, broken, c4.len, -1, nullptr);
		for (int k = 0; k < 2; k++) broken[rb.genes[k].offset + 12] ^= 0x01;
		chr_record_free(&rb);
		chr_record rb2;
		(void)chr_parse(&rb2, broken, c4.len, -1, nullptr);
		int damaged = 0;
		for (size_t k = 0; k < rb2.gene_count; k++)
			if (!rb2.genes[k].crc_ok) damaged++;
		printf("chr_par2_damaged %d %d\n", damaged, rb2.ngenes);
		vivi_bytes rec = { 0 };
		if (chr_read(&rec, &rb2, nullptr)) hexline("chr_par2_recover", rec.data, rec.len);
		else printf("chr_par2_recover fail\n");
		vivi_bytes_free(&rec);
		chr_record_free(&rb2);
		vivi_dealloc(broken);
	}

	/* the same corruption without parity cannot be recovered */
	uint8_t *plain = vivi_alloc(c1.len);
	if (plain) {
		memcpy(plain, c1.data, c1.len);
		chr_record rp;
		(void)chr_parse(&rp, plain, c1.len, -1, nullptr);
		for (int k = 0; k < 2; k++) plain[rp.genes[k].offset + 12] ^= 0x01;
		chr_record_free(&rp);
		chr_record rp2;
		(void)chr_parse(&rp2, plain, c1.len, -1, nullptr);
		vivi_bytes rec2 = { 0 };
		printf("chr_plain_recover %s\n", chr_read(&rec2, &rp2, nullptr) ? "ok" : "fail");
		vivi_bytes_free(&rec2);
		chr_record_free(&rp2);
		vivi_dealloc(plain);
	}

	/* generation rewrite keeps the CEN2 centromere */
	vivi_bytes c6 = { 0 };
	if (!chr_set_generation(&c6, c4.data, c4.len, 7, nullptr)) {
		printf("chr_par2_gen7 fail\n");
		return;
	}
	hexline("chr_par2_gen7", c6.data, c6.len);

	/* parity-only reconstruction: erase every data gene, keep the parity */
	chr_opts o7 = { 64, 0, 3, 3, 9, 4 };
	vivi_bytes c7 = { 0 };
	if (!chr_encode(&c7, 9, data, sizeof(data), &o7, nullptr)) {
		printf("chr_par4 fail\n");
		return;
	}
	uint8_t *broken2 = vivi_alloc(c7.len);
	if (broken2) {
		memcpy(broken2, c7.data, c7.len);
		chr_record rb3;
		(void)chr_parse(&rb3, broken2, c7.len, -1, nullptr);
		for (size_t k = 0; k < 4; k++) broken2[rb3.genes[k].offset + 12] ^= 0x01;
		chr_record_free(&rb3);
		chr_record rb4;
		(void)chr_parse(&rb4, broken2, c7.len, -1, nullptr);
		int alive_par = 0;
		for (size_t k = 0; k < rb4.gene_count; k++)
			if (rb4.genes[k].crc_ok && (rb4.genes[k].type & 2)) alive_par++;
		printf("chr_par4_alive %d %d\n", alive_par, rb4.ngenes);
		vivi_bytes rec3 = { 0 };
		if (chr_read(&rec3, &rb4, nullptr)) hexline("chr_par4_recover", rec3.data, rec3.len);
		else printf("chr_par4_recover fail\n");
		vivi_bytes_free(&rec3);
		chr_record_free(&rb4);
		vivi_dealloc(broken2);
	}

	vivi_bytes_free(&c1);
	vivi_bytes_free(&c2);
	vivi_bytes_free(&c3);
	vivi_bytes_free(&c4);
	vivi_bytes_free(&c5);
	vivi_bytes_free(&c6);
	vivi_bytes_free(&c7);
}

static void scenario_organisms(void)
{
	uint8_t d0[300], d1[180], e0[300], e1[180];
	payload_seed(d0, sizeof(d0), 1);
	payload_seed(d1, sizeof(d1), 2);
	payload_seed(e0, sizeof(e0), 3);
	payload_seed(e1, sizeof(e1), 4);

	int ids[2] = { 0, 1 };
	chr_opts opts[2] = { { 64, 0, 3, 3, 0, 0 }, { 64, 0, 3, 3, 5, 2 } };
	const uint8_t *datas[2] = { d0, d1 };
	size_t lens[2] = { sizeof(d0), sizeof(d1) };
	vivi_organism *org = nullptr;
	if (!organism_new(&org, ids, opts, datas, lens, 2, 40, nullptr)) {
		printf("org fail\n");
		return;
	}

	vivi_bytes s1 = { 0 }, s14 = { 0 }, s14n = { 0 };
	(void)organism_serialize(&s1, org, nullptr);
	(void)organism_serialize14(&s14, org, nullptr);
	(void)organism_serialize14n(&s14n, org, nullptr);
	hexline("org1", s1.data, s1.len);
	hexline("org14", s14.data, s14.len);
	hexline("org14n", s14n.data, s14n.len);
	printf("versions %d %d %d\n",
		organism_container_version(s14.data, s14.len),
		organism_container_version(s14n.data, s14n.len),
		organism_container_version((const uint8_t *)"VIV14NB4NSH33xx", 15));

	/* VIV14N roundtrip keeps the parity option */
	vivi_organism *loaded = nullptr;
	if (organism_deserialize(&loaded, s14n.data, s14n.len, nullptr)) {
		vivi_bytes rt = { 0 };
		(void)organism_serialize14n(&rt, loaded, nullptr);
		hexline("org14n_rt", rt.data, rt.len);
		printf("org14n_parity %d\n", loaded->chr_opts[1].parity);
		vivi_bytes_free(&rt);
		organism_free(loaded);
	} else {
		printf("org14n_rt fail\n");
	}

	/* VIV14 container carrying CEN2 strands infers the parity count */
	vivi_organism *inf = nullptr;
	if (organism_deserialize(&inf, s14.data, s14.len, nullptr)) {
		printf("org14_infer_parity %d\n", inf->chr_opts[1].parity);
		organism_free(inf);
	} else {
		printf("org14_infer_parity fail\n");
	}

	/* cross before the damage/replication scenarios mutate org */
	const uint8_t *datasB[2] = { e0, e1 };
	vivi_organism *orgB = nullptr;
	(void)organism_new(&orgB, ids, opts, datasB, lens, 2, 40, nullptr);
	vivi_organism *child = nullptr;
	if (organism_cross(&child, org, orgB, 3, nullptr)) {
		vivi_bytes sc = { 0 };
		(void)organism_serialize14n(&sc, child, nullptr);
		hexline("org14n_cross", sc.data, sc.len);
		vivi_bytes_free(&sc);
		organism_free(child);
	} else {
		printf("org14n_cross fail\n");
	}
	organism_free(orgB);

	(void)organism_damage(org, 3, 7, nullptr);
	vivi_bytes sd = { 0 };
	(void)organism_serialize14n(&sd, org, nullptr);
	hexline("org14n_damage", sd.data, sd.len);

	(void)organism_replicate(org, nullptr);
	vivi_bytes sr = { 0 };
	(void)organism_serialize14n(&sr, org, nullptr);
	hexline("org14n_replicate", sr.data, sr.len);

	vivi_organism *m = nullptr;
	if (organism_deserialize(&m, s14n.data, s14n.len, nullptr)) {
		org_mut_result mr = { 0 };
		if (organism_mutate(&mr, m, 2, 5, nullptr)) {
			printf("org_mutate_chr %d\n", mr.chr);
			vivi_bytes sm = { 0 };
			(void)organism_serialize14n(&sm, m, nullptr);
			hexline("org14n_mutate", sm.data, sm.len);
			vivi_bytes_free(&sm);
			vivi_bytes_free(&mr.data);
		} else {
			printf("org_mutate fail\n");
		}
		organism_free(m);
	} else {
		printf("org_mutate_load fail\n");
	}

	vivi_bytes_free(&s1);
	vivi_bytes_free(&s14);
	vivi_bytes_free(&s14n);
	vivi_bytes_free(&sd);
	vivi_bytes_free(&sr);
	organism_free(org);
}

int main(void)
{
	scenario_chromosomes();
	scenario_organisms();
	return 0;
}
