/* dna.c -- packed quaternary codec, exact transliteration of dna.lua */
#include "vivi/dna.h"
#include <string.h>

/* absorbs messages when a caller passes err == nullptr */
static const char *dna_err_sink;

static const uint32_t DNA_HDR_CONST = 0x1B87'3593u;

static uint32_t dna_xs32(uint32_t x)
{
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return x;
}

static uint32_t dna_seed_for(size_t len)
{
	uint32_t s = 0x9E37'79B9u ^ (uint32_t)len ^ 0x5BD1'E995u;
	if (s == 0) s = 1; /* xorshift32 is stuck at zero (len 0xC5E6902C) */
	return dna_xs32(s);
}

static uint32_t dna_hdr_seed(void) { return dna_xs32(DNA_HDR_CONST); }

/* keystream bytes: each byte is the next 8 keystream bits, MSB first */
typedef struct { uint32_t state, word; int kleft; } dna_ks;

/* base-4 digits */
enum : uint8_t { BASE_A = 0, BASE_C = 1, BASE_G = 2, BASE_T = 3 };

static void dna_ks_init(dna_ks *k, uint32_t seed)
{
	k->state = seed;
	k->word = 0;
	k->kleft = 0;
}

static uint8_t dna_ks_byte(dna_ks *k)
{
	if (k->kleft == 0) {
		k->state = dna_xs32(k->state);
		k->word = k->state;
		k->kleft = 32;
	}
	k->kleft -= 8;
	return (uint8_t)((k->word >> k->kleft) & 255u);
}

/* byte fast path: entry 0 = slow path, else 1 + nrun*128 + nlast*16 + ngc*2 */
static uint32_t *dna_ft_cache[13];

static const uint32_t *dna_fast_table(int h)
{
	if (h < 1 || h > 12) return nullptr;
	if (!dna_ft_cache[h]) {
		size_t count = (size_t)(h - 1) * 4 * 256;
		uint32_t *t = vivi_zalloc(count ? count : 1, sizeof(uint32_t));
		if (!t) return nullptr;
		for (int run = 1; run <= h - 1; run++) {
			for (int last = 0; last <= 3; last++) {
				for (int b = 0; b < 256; b++) {
					int r = run, p = last, g = 0;
					int good = 1;
					for (int sh = 6; sh >= 0; sh -= 2) {
						int d = (b >> sh) & 3;
						if (r >= h) { good = 0; break; }
						if (d == 1 || d == 2) g++;
						if (d == p) r++; else { r = 1; p = d; }
					}
					if (good)
						t[((run - 1) * 4 + last) * 256 + b] =
							(uint32_t)(1 + r * 128 + p * 16 + g * 2);
				}
			}
		}
		dna_ft_cache[h] = t;
	}
	return dna_ft_cache[h];
}

/* whitened byte at 1-based stream position w (header is fixed, payload
 * byte j is fetched as data[j] ^ keystream byte, strictly sequential) */
typedef struct {
	const uint8_t *data;
	size_t len;
	const uint8_t *WHB;
	dna_ks *kp;
	uint32_t pj, pbc;
} dna_wf;

static uint8_t dna_wfetch(dna_wf *f, size_t w)
{
	if (w <= 4) return f->WHB[w - 1];
	size_t j = w - 4;
	if (j > f->len) return 0;
	if (j != f->pj) {
		f->pj = (uint32_t)j;
		f->pbc = (uint8_t)(f->data[j - 1] ^ dna_ks_byte(f->kp));
	}
	return (uint8_t)f->pbc;
}

static int dna_gc_digit(int d) { return (d == 1 || d == 2) ? 1 : 0; }

bool dna_encode(vivi_bytes *out, const uint8_t *data, size_t len,
	const dna_opts *opts, const char **err)
{
	const char *unused_err;
	if (!err) err = &unused_err;
	*out = (vivi_bytes){ 0 };
	int h = (!opts || opts->h == 0) ? 3 : opts->h;
	double eps = (!opts || opts->gc_eps == 0) ? 0.05 : opts->gc_eps;
	if (!data && len) { *err = "expected string"; return false; }
	if (h < 1 || h > 12) { *err = "invalid h (integer 1..12)"; return false; }
	if (!(eps >= 0.005 && eps < 0.5)) { *err = "invalid gc_eps ([0.005, 0.5))"; return false; }
	if (len > 0xFFFF'FFFFull) { *err = "payload exceeds 4GB limit"; return false; }
	uint64_t nbits = 32 + 8ull * len;
	int hp = h - 1;
	int pd0 = hp % 4;
	int pd1 = ((hp / 4) + pd0 + 1) % 4;
	dna_ks kh, kp;
	dna_ks_init(&kh, dna_hdr_seed());
	dna_ks_init(&kp, dna_seed_for(len));
	uint8_t WHB[4];
	uint32_t lw = (uint32_t)len;
	for (int i = 0; i < 4; i++)
		WHB[i] = (uint8_t)(((lw >> (24 - 8 * i)) & 255u) ^ dna_ks_byte(&kh));
	dna_wf wf = { data, len, WHB, &kp, 0, 0 };
	const uint32_t *ft = dna_fast_table(h);
	if (!ft) { *err = "out of memory"; return false; }

	vivi_buf o = { 0 };
	if (!vivi_buf_reserve(&o, (size_t)(nbits / 2) / 4 + 64)) {
		*err = "out of memory";
		return false;
	}
	int acc = pd0 * 4 + pd1, nib = 2;
	int gc = dna_gc_digit(pd0) + dna_gc_digit(pd1);
	size_t wi = 1;
	uint8_t wb = WHB[0];
	int boff = 0;
	int last = pd1, run = h;   /* saturated: first base must differ from prefix */
	uint64_t bi = 0;
	int oom = 0;
	while (bi < nbits && !oom) {
		/* remaining is per segment (header | payload) so the
		 * dangling-bit rule mirrors the decoder */
		uint64_t remaining = (bi < 32) ? (32 - bi) : (nbits - bi);
		int d;
		int fast = 0;
		if (run < h && boff == 0 && nib == 0 && (bi >= 32 || bi + 8 <= 32)) {
			uint32_t e = ft[((run - 1) * 4 + last) * 256 + wb];
			if (e) {
				/* whole byte codes in 2-bit mode: emit it as is */
				if (!vivi_buf_push(&o, wb)) { oom = 1; break; }
				int q = (int)e - 1;
				run = (int)((q >> 7) & 15);
				last = (int)((q >> 4) & 3);
				gc += (int)((q >> 1) & 7);
				bi += 8;
				wi++;
				wb = dna_wfetch(&wf, wi);
				fast = 1;
			}
		}
		if (!fast) {
			if (run >= h || remaining == 1) {
				/* 1-bit alphabet over the base pair opposite to last */
				int bit = (wb >> (7 - boff)) & 1;
				boff++;
				if (boff == 8) { wi++; wb = dna_wfetch(&wf, wi); boff = 0; }
				if (last == 0 || last == 3) d = (bit == 0) ? 1 : 2;
				else d = (bit == 0) ? 0 : 3;
				run = 1;
				last = d;
				bi += 1;
			} else {
				if (boff <= 6) {
					d = (wb >> (6 - boff)) & 3;
					boff += 2;
					if (boff == 8) { wi++; wb = dna_wfetch(&wf, wi); boff = 0; }
				} else {
					uint8_t nb = dna_wfetch(&wf, wi + 1);
					d = (int)(((wb & 1u) << 1) | (nb >> 7));
					wi++;
					wb = nb;
					boff = 1;
				}
				bi += 2;
				if (d == last) run++; else { run = 1; last = d; }
			}
			acc = acc * 4 + d;
			nib++;
			if (nib == 4) {
				if (!vivi_buf_push(&o, (uint8_t)acc)) { oom = 1; break; }
				acc = 0;
				nib = 0;
			}
			if (d == 1 || d == 2) gc++;
		}
	}
	if (!oom) {
		/* finalize GC content: fill digits + heavy/safe padding so the
		 * whole strand, fills included, stays inside the GC band */
		size_t n = o.len * 4 + (size_t)nib;
		double lo_b = 0.5 - eps, hi_b = 0.5 + eps;
		const int AT[2] = { 0, 3 }, GCP[2] = { 1, 2 };
#define DNA_PAD(cnt, pair) \
	do { \
		for (int pc = 0; pc < (cnt) && !oom; pc++) { \
			int dd = (last == (pair)[0]) ? (pair)[1] : (pair)[0]; \
			acc = acc * 4 + dd; \
			nib++; \
			if (nib == 4) { \
				if (!vivi_buf_push(&o, (uint8_t)acc)) { oom = 1; break; } \
				acc = 0; nib = 0; \
			} \
			n++; last = dd; \
			if (dd == 1 || dd == 2) gc++; \
		} \
	} while (0)
		while (!oom) {
			int nfp = (nib > 0) ? (4 - nib) : 0;
			int done = 0;
			int cands[2][2] = { { AT[0], AT[1] }, { GCP[0], GCP[1] } };
			int ncand = (nfp > 0) ? 2 : 1;
			for (int ci = 0; ci < ncand && !done; ci++) {
				int addgc = (ncand == 2 && ci == 1) ? nfp : 0;
				double r = ((double)gc + addgc) / (double)(n + nfp);
				if (r >= lo_b && r <= hi_b) {
					if (nfp > 0) DNA_PAD(nfp, cands[ci]);
					done = 1;
				}
			}
			if (done) break;
			int below = ((double)gc / (double)n) < 0.5;
			int kk = (nfp > 0) ? nfp : 4;
			double r2 = ((double)gc + (below ? kk : 0)) / (double)(n + kk);
			if (below && r2 > hi_b) {
				DNA_PAD(1, GCP);   /* heavy step would overshoot; creep */
			} else if (!below && r2 < lo_b) {
				DNA_PAD(1, AT);
			} else {
				DNA_PAD(kk, below ? GCP : AT);
			}
		}
#undef DNA_PAD
	}
	if (!oom && nib > 0) {
		/* partial final byte: alternating fill bases that never extend
		 * the last run (not payload) */
		int f = (last == 0) ? 3 : 0;
		while (nib < 4 && !oom) {
			acc = acc * 4 + f;
			nib++;
			if (nib == 4) {
				if (!vivi_buf_push(&o, (uint8_t)acc)) oom = 1;
				else { acc = 0; nib = 0; }
			}
			f = 3 - f;
		}
	}
	if (oom) {
		vivi_buf_free(&o);
		*err = "out of memory";
		return false;
	}
	out->data = o.p;
	out->len = o.len;
	vivi_buf_release(&o);
	return true;
}

bool dna_decode(vivi_bytes *out, const uint8_t *strand, size_t slen,
	const char **err)
{
	*out = (vivi_bytes){ 0 };
	const char *unused_err; if (!err) err = &unused_err;
	*err = nullptr;
	if (!strand || slen < 1) { *err = "invalid strand"; return false; }
	int b0 = strand[0];
	int pd0 = b0 >> 6;
	int pd1 = (b0 >> 4) & 3;
	int fc = (((pd1 - pd0 - 1) % 4) + 4) % 4;
	if (fc > 2) { *err = "invalid strand prefix"; return false; }
	int h = pd0 + 4 * fc + 1;
	const uint32_t *ft = dna_fast_table(h);
	if (!ft) { *err = "out of memory"; return false; }
	dna_ks kh, kp;
	dna_ks_init(&kh, dna_hdr_seed());
	dna_ks_init(&kp, 0); /* real seed set when the header completes */
	vivi_buf bytes = { 0 };
	if (!vivi_buf_reserve(&bytes, 64)) { *err = "out of memory"; return false; }
	int acc8 = 0, nacc = 0;
	uint64_t nbitsout = 0;
	size_t need_bytes = 0;
	int need_known = 0;
	size_t total_bases = slen * 4;
	size_t k = 2;                 /* digits 0,1 are the prefix */
	size_t ci = 0;                /* 0-based byte cursor */
	int coff = 4;                 /* next digit offset within cb */
	uint8_t cb = strand[0];
	int last = pd1, run = h;      /* saturated, mirrors encode */
	int oom = 0;
#define DNA_FLUSH(w) \
	do { \
		uint8_t ksb = (bytes.len < 4) ? dna_ks_byte(&kh) : dna_ks_byte(&kp); \
		if (!vivi_buf_push(&bytes, (uint8_t)((w) ^ ksb))) { oom = 1; break; } \
		acc8 = 0; nacc = 0; \
		if (!need_known && bytes.len == 4) { \
			size_t L = ((size_t)bytes.p[0] << 24) | ((size_t)bytes.p[1] << 16) | \
				((size_t)bytes.p[2] << 8) | (size_t)bytes.p[3]; \
			need_bytes = 4 + L; \
			need_known = 1; \
			dna_ks_init(&kp, dna_seed_for(L)); \
		} \
	} while (0)
	while (!oom) {
		if (need_known && bytes.len >= need_bytes) break;
		if (k >= total_bases) { *err = "truncated strand"; break; }
		uint64_t before = nbitsout;
		uint64_t need_bits = need_known ? (uint64_t)need_bytes * 8 : 32;
		uint64_t remaining = (before < 32) ? (32 - before) : (need_bits - before);
		int fast = 0;
		if (run < h && nacc == 0 && coff == 0 && remaining >= 8
			&& (before >= 32 || before + 8 <= 32)) {
			uint32_t e = ft[((run - 1) * 4 + last) * 256 + cb];
			if (e) {
				/* whole byte decodes in 2-bit mode: one whitened byte */
				DNA_FLUSH(cb);
				int q = (int)e - 1;
				run = (int)((q >> 7) & 15);
				last = (int)((q >> 4) & 3);
				nbitsout = before + 8;
				k += 4;
				/* 4 digits always end at the next byte boundary */
				ci++;
				cb = (ci < slen) ? strand[ci] : 0;
				coff = 0;
				fast = 1;
			}
		}
		if (!fast) {
			int d = (cb >> (6 - coff)) & 3;
			coff += 2;
			if (coff == 8) {
				ci++;
				cb = (ci < slen) ? strand[ci] : 0;
				coff = 0;
			}
			k++;
			before = nbitsout;
			need_bits = need_known ? (uint64_t)need_bytes * 8 : 32;
			remaining = (before < 32) ? (32 - before) : (need_bits - before);
			if (run >= h || remaining == 1) {
				if (d == last) { *err = "homopolymer violation"; break; }
				int bt;
				if (last == 0 || last == 3) {
					if (d != 1 && d != 2) { *err = "invalid transition"; break; }
					bt = d - 1;
				} else {
					if (d != 0 && d != 3) { *err = "invalid transition"; break; }
					bt = (d == 3) ? 1 : 0;
				}
				acc8 = acc8 * 2 + bt;
				nacc++;
				nbitsout = before + 1;
				run = 1;
				last = d;
				if (nacc == 8) DNA_FLUSH(acc8);
			} else {
				if (nacc == 7) {
					/* 1 bit completes the byte, 1 bit starts the next */
					acc8 = acc8 * 2 + ((d >> 1) & 1);
					DNA_FLUSH(acc8);
					acc8 = d & 1;
					nacc = 1;
				} else {
					acc8 = acc8 * 4 + d;
					nacc += 2;
					if (nacc == 8) DNA_FLUSH(acc8);
				}
				nbitsout = before + 2;
				if (d == last) run++; else { run = 1; last = d; }
			}
		}
	}
#undef DNA_FLUSH
	if (oom) {
		vivi_buf_free(&bytes);
		*err = "out of memory";
		return false;
	}
	if (*err) {
		vivi_buf_free(&bytes);
		return false;
	}
	if (!need_known || bytes.len < need_bytes) {
		vivi_buf_free(&bytes);
		*err = "truncated strand";
		return false;
	}
	size_t paylen = need_bytes - 4;
	if (paylen > 0) {
		out->data = vivi_alloc(paylen);
		if (!out->data) {
			vivi_buf_free(&bytes);
			*err = "out of memory";
			return false;
		}
		memcpy(out->data, bytes.p + 4, paylen);
	}
	out->len = paylen;
	vivi_buf_free(&bytes);
	return true;
}

double dna_gc_content(const uint8_t *s, size_t n)
{
	size_t gc = 0;
	for (size_t i = 0; i < n; i++) {
		uint8_t b = s[i];
		for (int j = 0; j < 4; j++) {
			int d = (int)((b >> (2 * (3 - j))) & 3u);
			if (d == 1 || d == 2) gc++;
		}
	}
	if (n == 0) return 0.0;
	return (double)gc / (double)(n * 4);
}

size_t dna_max_homopolymer(const uint8_t *s, size_t n)
{
	size_t best = 0;
	int prevd = -1;
	size_t prevr = 0;
	for (size_t i = 0; i < n; i++) {
		uint8_t b = s[i];
		int d0 = (int)(b >> 6);
		int d1 = (int)((b >> 4) & 3u);
		int d2 = (int)((b >> 2) & 3u);
		int d3 = (int)(b & 3u);
		int run = 1, prev = d0, mr = 1;
		if (d1 == prev) run++; else { run = 1; prev = d1; }
		if (run > mr) mr = run;
		if (d2 == prev) run++; else { run = 1; prev = d2; }
		if (run > mr) mr = run;
		if (d3 == prev) run++; else run = 1;
		if (run > mr) mr = run;
		int sr = 1;
		if (d1 == d0) {
			sr = 2;
			if (d2 == d0) {
				sr = 3;
				if (d3 == d0) sr = 4;
			}
		}
		int er = 1;
		if (d3 == d2) {
			er = 2;
			if (d2 == d1) {
				er = 3;
				if (d1 == d0) er = 4;
			}
		}
		if ((size_t)mr > best) best = (size_t)mr;
		if (d0 == prevd && prevd >= 0) {
			size_t r = prevr + (size_t)sr;
			if (r > best) best = r;
			/* a uniform byte (sr == 4) continues the run */
			prevr = (sr == 4) ? prevr + 4 : (size_t)er;
		} else {
			prevr = (size_t)er;
		}
		prevd = d3;
	}
	return best;
}

bool dna_complement(vivi_bytes *out, const uint8_t *s, size_t n, const char **err)
{
	if (!err) err = &dna_err_sink;
	*out = (vivi_bytes){ 0 };
	if (!s && n) { *err = "expected string"; return false; }
	out->data = vivi_alloc(n ? n : 1);
	if (!out->data) { *err = "out of memory"; return false; }
	for (size_t i = 0; i < n; i++) out->data[i] = (uint8_t)(s[i] ^ 0xFFu);
	out->len = n;
	return true;
}

bool dna_reverse_complement(vivi_bytes *out, const uint8_t *s, size_t n, const char **err)
{
	if (!err) err = &dna_err_sink;
	*out = (vivi_bytes){ 0 };
	if (!s && n) { *err = "expected string"; return false; }
	out->data = vivi_alloc(n ? n : 1);
	if (!out->data) { *err = "out of memory"; return false; }
	for (size_t i = 0; i < n; i++)
		out->data[i] = (uint8_t)(s[n - 1 - i] ^ 0xFFu);
	out->len = n;
	return true;
}

static const char DNA_CHARS[4] = { 'A', 'C', 'G', 'T' };

bool dna_to_ascii(vivi_bytes *out, const uint8_t *s, size_t n, long nbases, const char **err)
{
	if (!err) err = &dna_err_sink;
	*out = (vivi_bytes){ 0 };
	if (!s && n) { *err = "expected string"; return false; }
	size_t total;
	if (nbases < 0) total = n * 4;
	else {
		if ((unsigned long)nbases > (unsigned long)(n * 4)) total = n * 4;
		else total = (size_t)nbases;
	}
	out->data = vivi_alloc(total ? total : 1);
	if (!out->data) { *err = "out of memory"; return false; }
	for (size_t i = 0; i < total; i++) {
		size_t bytei = i / 4;
		int off = (int)(i % 4);
		int d = (int)((s[bytei] >> (2 * (3 - off))) & 3u);
		out->data[i] = (uint8_t)DNA_CHARS[d];
	}
	out->len = total;
	return true;
}

bool dna_from_ascii(vivi_bytes *out, const char *s, size_t n, const char **err)
{
	if (!err) err = &dna_err_sink;
	*out = (vivi_bytes){ 0 };
	if (!s) { *err = "expected string"; return false; }
	size_t cap = (n + 3) / 4;
	uint8_t *t = vivi_alloc(cap ? cap : 1);
	if (!t) { *err = "out of memory"; return false; }
	size_t tn = 0;
	int acc = 0, k = 0;
	for (size_t i = 0; i < n; i++) {
		uint8_t c = (uint8_t)s[i];
		if (c == ' ' || c == '\t' || c == '\r' || c == '\n') continue;
		if (c >= 'a' && c <= 'z') c = (uint8_t)(c - 32);
		int d;
		switch (c) {
		case 'A': d = 0; break;
		case 'C': d = 1; break;
		case 'G': d = 2; break;
		case 'T': d = 3; break;
		default:
			vivi_dealloc(t);
			*err = "invalid base";
			return false;
		}
		acc = acc * 4 + d;
		k++;
		if (k == 4) {
			t[tn++] = (uint8_t)acc;
			acc = 0;
			k = 0;
		}
	}
	if (k > 0) {
		/* trailing partial group is left-aligned, like the Lua reference */
		t[tn++] = (uint8_t)(acc << (2 * (4 - k)));
	}
	out->data = t;
	out->len = tn;
	return true;
}

bool dna_validate(const uint8_t *s, size_t n, int h, double eps, const char **err)
{
	if (!err) err = &dna_err_sink;
	if (!s || n < 1) { *err = "invalid strand"; return false; }
	if (h < 1 || h > 12) { *err = "invalid h (integer 1..12)"; return false; }
	if (!(eps >= 0.005 && eps < 0.5)) { *err = "invalid gc_eps ([0.005, 0.5))"; return false; }
	size_t m = dna_max_homopolymer(s, n);
	if ((int)m > h) { *err = "homopolymer run exceeds h"; return false; }
	double g = dna_gc_content(s, n);
	if (g < 0.5 - eps || g > 0.5 + eps) { *err = "gc content outside band"; return false; }
	return true;
}

void dna_free_caches(void)
{
	for (int h = 1; h <= 12; h++) {
		vivi_dealloc(dna_ft_cache[h]);
		dna_ft_cache[h] = nullptr;
	}
}











