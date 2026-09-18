/* rs.c -- systematic Reed-Solomon over GF(256), byte-oriented research codec */
#include "vivi/rs.h"
#include <string.h>

#define GF_POLY 0x11Du
#define RS_MAX_T 128

static uint8_t gf_exp[512];
static uint8_t gf_log[256];
static int gf_ready;

static void gf_init(void)
{
	if (gf_ready) return;
	uint32_t x = 1;
	for (int i = 0; i < 255; i++) {
		gf_exp[i] = (uint8_t)x;
		gf_log[x] = (uint8_t)i;
		x <<= 1;
		if (x & 0x100u) x ^= GF_POLY;
	}
	for (int i = 255; i < 512; i++) gf_exp[i] = gf_exp[i - 255];
	gf_ready = 1;
}

static uint8_t gf_mul(uint8_t a, uint8_t b)
{
	if (!a || !b) return 0;
	int t = (int)gf_log[a] + (int)gf_log[b];
	if (t >= 255) t -= 255;
	return gf_exp[t];
}

static uint8_t gf_div(uint8_t a, uint8_t b)
{
	int t = (int)gf_log[a] - (int)gf_log[b];
	if (t < 0) t += 255;
	return gf_exp[t];
}

static uint8_t gf_inv(uint8_t a)
{
	return gf_exp[255 - (int)gf_log[a]];
}

static uint8_t gf_pow(uint8_t a, unsigned e)
{
	if (a == 0) return 0;
	return gf_exp[(int)((unsigned)gf_log[a] * e % 255u)];
}

/* Horner evaluation of a big-endian vector: sum v[i] * x^(n-1-i) */
static uint8_t poly_eval(const uint8_t *v, size_t n, uint8_t x)
{
	uint8_t y = 0;
	for (size_t i = 0; i < n; i++) y = (uint8_t)(gf_mul(y, x) ^ v[i]);
	return y;
}

/* out = a * b, low-to-high; na, nb are coefficient counts */
static void poly_mul(const uint8_t *a, size_t na, const uint8_t *b, size_t nb, uint8_t *out)
{
	size_t nc = na + nb - 1;
	for (size_t i = 0; i < nc; i++) out[i] = 0;
	for (size_t i = 0; i < na; i++)
		for (size_t j = 0; j < nb; j++)
			out[i + j] ^= gf_mul(a[i], b[j]);
}

/* unique degree < m polynomial through (xs[j], ys[j]) */
static void lagrange_interp(const uint8_t *xs, const uint8_t *ys, size_t m, uint8_t *out)
{
	uint8_t basis[RS_MAX_N], tmp[RS_MAX_N + 2];
	for (size_t i = 0; i < m; i++) out[i] = 0;
	for (size_t j = 0; j < m; j++) {
		size_t nb = 1;
		basis[0] = 1;
		uint8_t den = 1;
		for (size_t l = 0; l < m; l++) {
			if (l == j) continue;
			uint8_t lin[2] = { xs[l], 1 };
			poly_mul(basis, nb, lin, 2, tmp);
			nb++;
			memcpy(basis, tmp, nb);
			den = gf_mul(den, (uint8_t)(xs[j] ^ xs[l]));
		}
		uint8_t scale = gf_mul(ys[j], gf_inv(den));
		for (size_t t = 0; t < m; t++) out[t] ^= gf_mul(basis[t], scale);
	}
}

bool rs_encode(vivi_bytes *out, const uint8_t *data, size_t k, size_t m,
	const char **err)
{
	static const char *sink;
	if (!err) err = &sink;
	if (out) *out = (vivi_bytes){ 0 };
	if (!out || (!data && k)) { *err = "expected buffers"; return false; }
	if (m < 1 || k + m > RS_MAX_N) { *err = "invalid rs geometry"; return false; }
	gf_init();
	uint8_t xs[RS_MAX_N], ys[RS_MAX_N], parity[RS_MAX_N];
	for (size_t j = 0; j < m; j++) {
		xs[j] = gf_exp[j];
		ys[j] = gf_mul(poly_eval(data, k, xs[j]), gf_pow(xs[j], (unsigned)m));
	}
	lagrange_interp(xs, ys, m, parity);
	out->data = vivi_alloc(k + m);
	if (!out->data) { *err = "out of memory"; return false; }
	if (k) memcpy(out->data, data, k);
	for (size_t r = 0; r < m; r++) out->data[k + r] = parity[m - 1 - r];
	out->len = k + m;
	return true;
}

static void syndromes(const uint8_t *cw, size_t n, size_t m, uint8_t *S)
{
	for (size_t j = 0; j < m; j++) S[j] = poly_eval(cw, n, gf_exp[j]);
}

static void berlekamp_massey(const uint8_t *S, size_t m, uint8_t *C, size_t *out_L)
{
	uint8_t B[RS_MAX_N], T[RS_MAX_N];
	for (size_t i = 0; i < m; i++) { C[i] = 0; B[i] = 0; }
	C[0] = 1;
	B[0] = 1;
	size_t L = 0, shift = 1;
	uint8_t b = 1;
	for (size_t nn = 0; nn < m; nn++) {
		uint8_t d = S[nn];
		for (size_t i = 1; i <= L && i <= nn; i++) d ^= gf_mul(C[i], S[nn - i]);
		if (d == 0) { shift++; continue; }
		uint8_t coef = gf_div(d, b);
		if (2 * L <= nn) {
			memcpy(T, C, m * sizeof(uint8_t));
			for (size_t i = 0; i + shift < m; i++)
				C[i + shift] ^= gf_mul(coef, B[i]);
			L = nn + 1 - L;
			memcpy(B, T, m * sizeof(uint8_t));
			b = d;
			shift = 1;
		} else {
			for (size_t i = 0; i + shift < m; i++)
				C[i + shift] ^= gf_mul(coef, B[i]);
			shift++;
		}
	}
	*out_L = L;
}

/* solve S[j] = sum_l e_l * X_l^j for j < L by Gauss-Jordan */
static bool solve_values(const uint8_t *S, const uint8_t *X, size_t L, uint8_t *e)
{
	uint8_t A[RS_MAX_T][RS_MAX_T + 1];
	for (size_t l = 0; l < L; l++) {
		uint8_t pw = 1;
		for (size_t j = 0; j < L; j++) { A[j][l] = pw; pw = gf_mul(pw, X[l]); }
	}
	for (size_t j = 0; j < L; j++) A[j][L] = S[j];
	for (size_t col = 0; col < L; col++) {
		size_t piv = L;
		for (size_t r = col; r < L; r++) if (A[r][col]) { piv = r; break; }
		if (piv == L) return false;
		if (piv != col)
			for (size_t c = col; c <= L; c++) {
				uint8_t t = A[col][c]; A[col][c] = A[piv][c]; A[piv][c] = t;
			}
		uint8_t inv = gf_inv(A[col][col]);
		for (size_t c = col; c <= L; c++) A[col][c] = gf_mul(A[col][c], inv);
		for (size_t r = 0; r < L; r++) {
			if (r == col || !A[r][col]) continue;
			uint8_t f = A[r][col];
			for (size_t c = col; c <= L; c++) A[r][c] ^= gf_mul(f, A[col][c]);
		}
	}
	for (size_t i = 0; i < L; i++) e[i] = A[i][L];
	return true;
}

bool rs_decode(uint8_t *cw, size_t n, size_t k, size_t *corrected, const char **err)
{
	static const char *sink;
	if (!err) err = &sink;
	if (corrected) *corrected = 0;
	if (!cw || n > RS_MAX_N || k < 1 || k >= n) { *err = "invalid rs geometry"; return false; }
	size_t m = n - k;
	gf_init();
	uint8_t S[RS_MAX_N];
	syndromes(cw, n, m, S);
	int nonzero = 0;
	for (size_t j = 0; j < m; j++) if (S[j]) nonzero = 1;
	if (!nonzero) return true;
	uint8_t C[RS_MAX_N];
	size_t L = 0;
	berlekamp_massey(S, m, C, &L);
	if (L == 0 || L > m / 2) { *err = "too many errors"; return false; }
	uint8_t pos[RS_MAX_T], X[RS_MAX_T], e[RS_MAX_T];
	size_t found = 0;
	for (size_t i = 0; i < n; i++) {
		uint8_t x = gf_exp[(unsigned)(n - 1 - i)];
		uint8_t xinv = gf_inv(x);
		uint8_t val = 0, pw = 1;
		for (size_t t = 0; t <= L; t++) { val ^= gf_mul(C[t], pw); pw = gf_mul(pw, xinv); }
		if (val == 0) {
			if (found >= L) { *err = "too many roots"; return false; }
			pos[found] = (uint8_t)i;
			X[found] = x;
			found++;
		}
	}
	if (found != L) { *err = "error locator has no roots"; return false; }
	if (!solve_values(S, X, L, e)) { *err = "singular system"; return false; }
	for (size_t l = 0; l < L; l++) cw[pos[l]] ^= e[l];
	uint8_t S2[RS_MAX_N];
	syndromes(cw, n, m, S2);
	for (size_t j = 0; j < m; j++)
		if (S2[j]) { *err = "miscorrection"; return false; }
	if (corrected) *corrected = L;
	return true;
}
