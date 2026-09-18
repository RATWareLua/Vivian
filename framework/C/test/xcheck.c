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
#include "vivi/rs.h"
#include "vivi/channel.h"
#include "vivi/pool.h"
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

	chr_opts o1 = { 64, 0, 3, 3, 9, 0, 0, 0 };   /* gene_raw, codon, h, units, flags, parity */
	vivi_bytes c1 = { 0 };
	if (!chr_encode(&c1, 7, data, sizeof(data), &o1, nullptr)) {
		printf("chr_dense fail\n");
		return;
	}
	hexline("chr_dense", c1.data, c1.len);

	chr_opts o2 = { 64, 1, 3, 3, 9, 0, 0, 0 };
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

	chr_opts o4 = { 64, 0, 3, 3, 9, 2, 0, 0 };
	vivi_bytes c4 = { 0 };
	if (!chr_encode(&c4, 7, data, sizeof(data), &o4, nullptr)) {
		printf("chr_par2 fail\n");
		return;
	}
	hexline("chr_par2", c4.data, c4.len);

	chr_opts o5 = { 48, 1, 3, 3, 9, 3, 0, 0 };
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
	chr_opts o7 = { 64, 0, 3, 3, 9, 4, 0, 0 };
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

	/* VIV14NB4NSH33: on-strand primer sites */
	chr_opts o8 = { 64, 0, 3, 3, 9, 2, 1234, 0 };
	vivi_bytes c8 = { 0 };
	if (!chr_encode(&c8, 11, data, sizeof(data), &o8, nullptr)) {
		printf("chr_banshee fail\n");
		return;
	}
	hexline("chr_banshee", c8.data, c8.len);
	chr_record r8;
	if (!chr_parse(&r8, c8.data, c8.len, -1, nullptr)) {
		printf("chr_banshee_fields fail\n");
		return;
	}
	printf("chr_banshee_fields %d %d %d %d %d %d\n", r8.cen_version, r8.parity,
		r8.primer, r8.primer_ok, r8.ngenes, (int)r8.gene_count);
	chr_record_free(&r8);

	/* destroy the reverse site: data stays readable, access is lost */
	uint8_t *pb = vivi_alloc(c8.len);
	if (pb) {
		memcpy(pb, c8.data, c8.len);
		pb[c8.len - 24] ^= 0x40;   /* first marker byte of the reverse site */
		chr_record r9;
		(void)chr_parse(&r9, pb, c8.len, -1, nullptr);
		printf("chr_banshee_dark %d %d %d\n", r9.primer, r9.primer_ok,
			chr_amplifiable(&r9) ? 1 : 0);
		vivi_bytes rr = { 0 };
		if (chr_read(&rr, &r9, nullptr)) hexline("chr_banshee_dark_read", rr.data, rr.len);
		else printf("chr_banshee_dark_read fail\n");
		vivi_bytes_free(&rr);
		chr_record_free(&r9);
		vivi_dealloc(pb);
	}

	vivi_bytes c9 = { 0 };
	if (!chr_set_generation(&c9, c8.data, c8.len, 11, nullptr)) {
		printf("chr_banshee_gen11 fail\n");
		return;
	}
	hexline("chr_banshee_gen11", c9.data, c9.len);

	vivi_bytes_free(&c8);
	vivi_bytes_free(&c9);

	/* inner Reed-Solomon genes */
	chr_opts oi = { 64, 0, 3, 3, 9, 0, 0, 16 };
	vivi_bytes ci = { 0 };
	if (!chr_encode(&ci, 13, data, sizeof(data), &oi, nullptr)) {
		printf("chr_inner fail\n");
		return;
	}
	hexline("chr_inner", ci.data, ci.len);

	chr_record ric;
	if (chr_parse(&ric, ci.data, ci.len, -1, nullptr)) {
		vivi_bytes ird = { 0 };
		if (chr_read(&ird, &ric, nullptr)) hexline("chr_inner_read", ird.data, ird.len);
		else printf("chr_inner_read fail\n");
		vivi_bytes_free(&ird);
		chr_record_free(&ric);
	} else {
		printf("chr_inner_read fail\n");
	}

	chr_record ri;
	if (!chr_parse(&ri, ci.data, ci.len, -1, nullptr)) {
		printf("chr_inner_fields fail\n");
		return;
	}
	printf("chr_inner_fields %d %d %zu\n", ri.cen_version, ri.inner, ri.gene_count);
	chr_record_free(&ri);

	uint8_t *ib = vivi_alloc(ci.len);
	if (ib) {
		memcpy(ib, ci.data, ci.len);
		chr_record rib;
		(void)chr_parse(&rib, ib, ci.len, -1, nullptr);
		ib[rib.genes[0].offset + 12] ^= 0x01;
		chr_record_free(&rib);
		chr_record rib2;
		(void)chr_parse(&rib2, ib, ci.len, -1, nullptr);
		int idamaged = 0, ifixed = 0;
		for (size_t k = 0; k < rib2.gene_count; k++) {
			if (!rib2.genes[k].crc_ok) idamaged++;
			ifixed += rib2.genes[k].inner_fixed;
		}
		printf("chr_inner_damaged %d %d %d\n", idamaged, ifixed, rib2.ngenes);
		vivi_bytes irec = { 0 };
		if (chr_read(&irec, &rib2, nullptr)) hexline("chr_inner_recover", irec.data, irec.len);
		else printf("chr_inner_recover fail\n");
		vivi_bytes_free(&irec);
		chr_record_free(&rib2);
		vivi_dealloc(ib);
	}

	chr_opts oip = { 64, 0, 3, 3, 9, 2, 0, 16 };
	vivi_bytes cip = { 0 };
	if (!chr_encode(&cip, 14, data, sizeof(data), &oip, nullptr)) {
		printf("chr_inner_parity fail\n");
		return;
	}
	hexline("chr_inner_parity", cip.data, cip.len);
	vivi_bytes_free(&ci);
	vivi_bytes_free(&cip);

	/* inner RS codec directly: encode, correct two bytes, reject three */
	uint8_t kdata[16];
	payload_seed(kdata, sizeof(kdata), 0x5EEDu);
	vivi_bytes rscw = { 0 };
	if (!rs_encode(&rscw, kdata, sizeof(kdata), 4, nullptr)) {
		printf("rs_encode fail\n");
		return;
	}
	hexline("rs_codeword", rscw.data, rscw.len);
	uint8_t rsfix[32];
	memcpy(rsfix, rscw.data, rscw.len);
	rsfix[2] ^= 0x11;
	rsfix[9] ^= 0x22;
	size_t rscor = 0;
	if (!rs_decode(rsfix, rscw.len, sizeof(kdata), &rscor, nullptr)) {
		printf("rs_correct fail\n");
	} else {
		printf("rs_correct %zu\n", rscor);
		hexline("rs_corrected", rsfix, sizeof(kdata));
	}
	memcpy(rsfix, rscw.data, rscw.len);
	rsfix[0] ^= 0x07;
	rsfix[5] ^= 0x08;
	rsfix[12] ^= 0x09;
	size_t rscor3 = 0;
	int rsok3 = rs_decode(rsfix, rscw.len, sizeof(kdata), &rscor3, nullptr);
	printf("rs_toomany %s\n", rsok3 ? "ok" : "fail");
	if (rsok3) hexline("rs_toomany_data", rsfix, sizeof(kdata));
	vivi_bytes_free(&rscw);

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
	chr_opts opts[2] = { { 64, 0, 3, 3, 0, 0, 0, 0 }, { 64, 0, 3, 3, 5, 2, 0, 0 } };
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
		printf("org14n_parity %d\n", organism_chr_opts_at(loaded, 1)->parity);
		vivi_bytes_free(&rt);
		organism_free(loaded);
	} else {
		printf("org14n_rt fail\n");
	}

	/* VIV14 container carrying CEN2 strands infers the parity count */
	vivi_organism *inf = nullptr;
	if (organism_deserialize(&inf, s14.data, s14.len, nullptr)) {
		printf("org14_infer_parity %d\n", organism_chr_opts_at(inf, 1)->parity);
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

	/* VIV14NB4NSH33 container */
	chr_opts bopts[2] = { { 64, 0, 3, 3, 0, 0, 0, 0 }, { 64, 0, 3, 3, 5, 2, 1234, 0 } };
	vivi_organism *bo = nullptr;
	if (organism_new(&bo, ids, bopts, datas, lens, 2, 40, nullptr)) {
		vivi_bytes bser = { 0 };
		if (organism_serialize14nb(&bser, bo, nullptr)) {
			hexline("org14nb", bser.data, bser.len);
			vivi_organism *bl = nullptr;
			if (organism_deserialize(&bl, bser.data, bser.len, nullptr)) {
				vivi_bytes brt = { 0 };
				(void)organism_serialize14nb(&brt, bl, nullptr);
				hexline("org14nb_rt", brt.data, brt.len);
				printf("org14nb_opts %d %d\n", organism_chr_opts_at(bl, 1)->parity,
					organism_chr_opts_at(bl, 1)->primer);
				vivi_bytes_free(&brt);
				organism_free(bl);
			} else {
				printf("org14nb_load fail\n");
			}
			vivi_bytes_free(&bser);
		} else {
			printf("org14nb fail\n");
		}
		organism_free(bo);
	} else {
		printf("org14nb_org fail\n");
	}

	/* inner Reed-Solomon inferred from the genes */
	chr_opts iopts[2] = { { 64, 0, 3, 3, 0, 0, 0, 0 },
		{ 64, 0, 3, 3, 5, 2, 0, 16 } };
	vivi_organism *io = nullptr;
	if (organism_new(&io, ids, iopts, datas, lens, 2, 40, nullptr)) {
		vivi_bytes iser = { 0 };
		if (organism_serialize14nb(&iser, io, nullptr)) {
			hexline("org14nb_inner", iser.data, iser.len);
			vivi_organism *il = nullptr;
			if (organism_deserialize(&il, iser.data, iser.len, nullptr)) {
				printf("org14nb_inner_infer %d\n", organism_chr_opts_at(il, 1)->inner);
				organism_free(il);
			} else {
				printf("org14nb_inner_load fail\n");
			}
			vivi_bytes_free(&iser);
		} else {
			printf("org14nb_inner fail\n");
		}
		organism_free(io);
	} else {
		printf("org14nb_inner_org fail\n");
	}

	organism_free(org);
}

static void chan_case(const char *label, const uint8_t *data, size_t n,
	const vivi_channel_opts *opts, int hexout, int numout)
{
	const char *err = nullptr;
	vivi_read rd;
	if (!vivi_channel_read(&rd, data, n, opts, &err)) {
		printf("%s fail\n", label);
		return;
	}
	if (hexout) hexline(label, rd.strand.data, rd.strand.len);
	if (numout) printf("%s_n %zu %d\n", label, rd.bases, rd.dropped);
	vivi_bytes_free(&rd.strand);
}

static void consensus_case(const char *label, const uint8_t *data, size_t n,
	const vivi_consensus_opts *opts, int hexout, int numout)
{
	const char *err = nullptr;
	vivi_read rd;
	if (!vivi_consensus_read(&rd, data, n, opts, &err)) {
		printf("%s fail\n", label);
		return;
	}
	if (hexout) hexline(label, rd.strand.data, rd.strand.len);
	if (numout) printf("%s_n %zu %d\n", label, rd.bases, rd.dropped);
	vivi_bytes_free(&rd.strand);
}

static void scenario_research(void)
{
	const char *err = nullptr;
	uint8_t chdata[32];
	payload_seed(chdata, sizeof(chdata), 0x51AB1Eu);
	hexline("ch_payload", chdata, sizeof(chdata));

	chan_case("ch_identity", chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 7 }, 1, 1);
	chan_case("ch_drop", chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 7, .p_drop = 0.5 }, 1, 1);
	chan_case("ch_sub", chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 7, .p_sub = 0.1 }, 1, 1);
	chan_case("ch_ins", chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 7, .p_ins = 0.1 }, 1, 1);
	chan_case("ch_del", chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 7, .p_del = 0.1 }, 1, 1);
	chan_case("ch_trunc", chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 7, .p_trunc = 1.0 }, 0, 1);
	chan_case("ch_burst", chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 7, .p_burst = 1.0, .burst_len = 4 }, 1, 0);
	chan_case("ch_gc", chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 7, .p_sub_gc = 1.0 }, 1, 0);
	chan_case("ch_hp", chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 7, .p_sub_hp = 1.0 }, 1, 0);
	chan_case("ch_combo", chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 11, .p_sub = 0.03, .p_ins = 0.02,
			.p_del = 0.02, .p_sub_gc = 0.05, .p_sub_hp = 0.05,
			.p_burst = 0.04, .burst_len = 3, .p_trunc = 0.1 }, 1, 0);

	consensus_case("ch_cons0", chdata, sizeof(chdata),
		&(vivi_consensus_opts){ .ch = { .seed = 7 }, .coverage = 0 }, 1, 0);
	consensus_case("ch_cons1", chdata, sizeof(chdata),
		&(vivi_consensus_opts){ .ch = { .seed = 7 }, .coverage = 1 }, 1, 0);

	vivi_read p0, c0, c1;
	(void)vivi_channel_read(&p0, chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 7 }, &err);
	(void)vivi_consensus_read(&c0, chdata, sizeof(chdata),
		&(vivi_consensus_opts){ .ch = { .seed = 7 }, .coverage = 0 }, &err);
	(void)vivi_consensus_read(&c1, chdata, sizeof(chdata),
		&(vivi_consensus_opts){ .ch = { .seed = 7 }, .coverage = 1 }, &err);
	int eq0 = !p0.dropped && p0.bases == sizeof(chdata) * 4
		&& p0.strand.len == sizeof(chdata)
		&& memcmp(p0.strand.data, chdata, sizeof(chdata)) == 0;
	int eq1 = !c0.dropped && c0.bases == sizeof(chdata) * 4
		&& c0.strand.len == sizeof(chdata)
		&& memcmp(c0.strand.data, chdata, sizeof(chdata)) == 0;
	int eq2 = !c1.dropped && c1.bases == sizeof(chdata) * 4
		&& c1.strand.len == sizeof(chdata)
		&& memcmp(c1.strand.data, chdata, sizeof(chdata)) == 0;
	printf("ch_cons_eq %d %d %d\n", eq0, eq1, eq2);
	vivi_bytes_free(&p0.strand);
	vivi_bytes_free(&c0.strand);
	vivi_bytes_free(&c1.strand);

	consensus_case("ch_cons5", chdata, sizeof(chdata),
		&(vivi_consensus_opts){ .ch = { .seed = 7, .p_sub = 0.05 }, .coverage = 5 }, 1, 0);

	vivi_read ir;
	if (vivi_consensus_read(&ir, chdata, sizeof(chdata),
			&(vivi_consensus_opts){ .ch = { .seed = 7, .p_ins = 1.0 }, .coverage = 5 }, &err)) {
		printf("ch_cons_indel ok\n");
		vivi_bytes_free(&ir.strand);
	} else {
		printf("ch_cons_indel %s\n", err);
	}

	uint8_t pdata[200];
	payload_seed(pdata, sizeof(pdata), 0xC0FFEE07u);
	vivi_bytes pchr = { 0 };
	chr_opts po = { .gene_raw = 64, .h = 3, .units = 3, .flags = 9 };
	if (!chr_encode(&pchr, 7, pdata, sizeof(pdata), &po, &err)) {
		printf("ch_pool_chr fail\n");
		return;
	}
	vivi_pool pool = { 0 };
	if (!vivi_pool_add_chromosome(&pool, pchr.data, pchr.len, &err)) {
		printf("ch_pool_add fail\n");
		vivi_bytes_free(&pchr);
		return;
	}
	vivi_amp_opts ao = { .p_access = 0.0, .p_cross = 0.01, .seed = 0x1234,
		.p_primer = 0.05, .coverage = 3 };
	ao.ch.p_sub = 0.01;
	vivi_amp_result ar;
	if (!vivi_pool_amplify(&ar, &pool, 1, &ao, &err)) {
		printf("ch_pool_amp fail\n");
	} else {
		hexline("ch_pool", ar.read.strand.data, ar.read.strand.len);
		printf("ch_pool_n %d %d\n", ar.id, ar.read.dropped);
		vivi_bytes_free(&ar.read.strand);
	}

	vivi_read srd;
	if (vivi_channel_read_soft(&srd, chdata, sizeof(chdata),
			&(vivi_channel_opts){ .seed = 7 }, &err)) {
		hexline("ch_soft_clean", srd.strand.data, srd.strand.len);
		int hi = 0;
		for (size_t i = 0; i < srd.bases; i++)
			if (srd.qual[i] == VIVI_QUAL_HI) hi++;
		printf("ch_soft_clean_n %zu %d %d\n", srd.bases, srd.dropped, hi);
		vivi_read_free(&srd);
	} else {
		printf("ch_soft_clean fail\n");
	}

	vivi_read srd2;
	if (vivi_channel_read_soft(&srd2, chdata, sizeof(chdata),
			&(vivi_channel_opts){ .seed = 7, .p_sub = 1.0 }, &err)) {
		int lo = 0;
		for (size_t i = 0; i < srd2.bases; i++)
			if (srd2.qual[i] == VIVI_QUAL_LO) lo++;
		printf("ch_soft_sub_n %d\n", lo);
		vivi_read_free(&srd2);
	} else {
		printf("ch_soft_sub fail\n");
	}

	vivi_read hardc, softc;
	int okh = vivi_consensus_read(&hardc, chdata, sizeof(chdata),
		&(vivi_consensus_opts){ .ch = { .seed = 7, .p_sub = 0.3 },
			.coverage = 3, .soft = 0 }, &err);
	int oks = vivi_consensus_read(&softc, chdata, sizeof(chdata),
		&(vivi_consensus_opts){ .ch = { .seed = 7, .p_sub = 0.3 },
			.coverage = 3, .soft = 1 }, &err);
	if (okh && oks) {
		int hardm = 0, softm = 0;
		for (size_t i = 0; i < sizeof(chdata); i++) {
			if (hardc.strand.data[i] != chdata[i]) hardm++;
			if (softc.strand.data[i] != chdata[i]) softm++;
		}
		printf("ch_soft_cons %d %d\n", hardm, softm);
	} else {
		printf("ch_soft_cons fail\n");
	}
	if (okh) vivi_read_free(&hardc);
	if (oks) vivi_read_free(&softc);

	chan_case("ch_burst_del", chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 7, .p_burst = 1.0, .burst_len = 4,
			.p_burst_del = 1.0 }, 0, 1);
	chan_case("ch_burst_del_ctrl", chdata, sizeof(chdata),
		&(vivi_channel_opts){ .seed = 7, .p_burst = 1.0, .burst_len = 4,
			.p_burst_del = 0.0 }, 0, 1);

	vivi_amp_opts sao = { .p_access = 0.0, .p_cross = 0.01, .seed = 0x1234,
		.p_primer = 0.05, .coverage = 3, .soft = 1 };
	sao.ch.p_sub = 0.01;
	vivi_amp_result sar;
	if (!vivi_pool_amplify(&sar, &pool, 1, &sao, &err)) {
		printf("ch_soft_pool_amp fail\n");
	} else {
		hexline("ch_soft_pool", sar.read.strand.data, sar.read.strand.len);
		printf("ch_soft_pool_n %d %d\n", sar.id, sar.read.dropped);
		vivi_read_free(&sar.read);
	}

	vivi_pool_free(&pool);
	vivi_bytes_free(&pchr);
}

int main(void)
{
	scenario_chromosomes();
	scenario_organisms();
	scenario_research();
	return 0;
}
