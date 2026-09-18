/* vivi.c -- execution context, shared buffer, allocator hooks and rng helpers */
#include <string.h>
#include "vivi/vivi.h"

struct vivi_context {
	vivi_alloc_fn alloc_fn;
	vivi_free_fn free_fn;
	void *dna_ft[13];
};

#ifdef VIVI_NO_HOSTED
static vivi_context vivi_default_ctx = { 0 };
#else
#include <stdlib.h>
static vivi_context vivi_default_ctx = { malloc, free, { 0 } };
#endif

static VIVI_THREAD_LOCAL vivi_context *vivi_tls;

static vivi_context *vivi_current_raw(void)
{
	return vivi_tls ? vivi_tls : &vivi_default_ctx;
}

void vivi_set_allocator(vivi_alloc_fn alloc_fn, vivi_free_fn free_fn)
{
	vivi_default_ctx.alloc_fn = alloc_fn;
	vivi_default_ctx.free_fn = free_fn;
}

vivi_context *vivi_context_new(vivi_alloc_fn alloc_fn, vivi_free_fn free_fn)
{
	if (!alloc_fn) return nullptr;
	vivi_context *ctx = alloc_fn(sizeof(vivi_context));
	if (!ctx) return nullptr;
	memset(ctx, 0, sizeof(*ctx));
	ctx->alloc_fn = alloc_fn;
	ctx->free_fn = free_fn;
	return ctx;
}

void vivi_context_release_caches(vivi_context *ctx)
{
	if (!ctx) return;
	vivi_free_fn f = ctx->free_fn;
	for (int h = 0; h < 13; h++) {
		if (ctx->dna_ft[h] && f) f(ctx->dna_ft[h]);
		ctx->dna_ft[h] = nullptr;
	}
}

void vivi_context_free(vivi_context *ctx)
{
	if (!ctx || ctx == &vivi_default_ctx) return;
	vivi_context_release_caches(ctx);
	vivi_free_fn f = ctx->free_fn;
	if (f) f(ctx);
}

vivi_context *vivi_context_enter(vivi_context *ctx)
{
	vivi_context *prev = vivi_tls;
	vivi_tls = ctx;
	return prev;
}

vivi_context *vivi_context_leave(vivi_context *prev)
{
	vivi_context *cur = vivi_tls;
	vivi_tls = prev;
	return cur;
}

vivi_context *vivi_context_current(void)
{
	return vivi_current_raw();
}

void **vivi_context_dna_cache(vivi_context *ctx)
{
	return ctx->dna_ft;
}

void *vivi_alloc(size_t size)
{
	vivi_context *ctx = vivi_current_raw();
	return ctx->alloc_fn ? ctx->alloc_fn(size) : nullptr;
}

void *vivi_zalloc(size_t n, size_t size)
{
	if (size && n > ((size_t)-1) / size) return nullptr;
	void *p = vivi_alloc(n * size);
	if (p) memset(p, 0, n * size);
	return p;
}

void vivi_dealloc(void *p)
{
	vivi_context *ctx = vivi_current_raw();
	if (ctx->free_fn) ctx->free_fn(p);
}

bool vivi_buf_reserve(vivi_buf *b, size_t need)
{
	if (b->cap >= need) return true;
	size_t cap = b->cap ? b->cap : 64;
	while (cap < need && cap <= ((size_t)-1) / 2) cap *= 2;
	if (cap < need) cap = need; /* near SIZE_MAX: exact fit, the alloc fails cleanly */
	uint8_t *p = vivi_alloc(cap);
	if (!p) return false;
	if (b->len) memcpy(p, b->p, b->len);
	vivi_dealloc(b->p);
	b->p = p;
	b->cap = cap;
	return true;
}

bool vivi_buf_push(vivi_buf *b, uint8_t c)
{
	if (!vivi_buf_reserve(b, b->len + 1)) return false;
	b->p[b->len++] = c;
	return true;
}

bool vivi_buf_append(vivi_buf *b, const void *data, size_t len)
{
	if (!vivi_buf_reserve(b, b->len + len)) return false;
	if (len) {
		uint8_t *dst = b->p + b->len;
		const uint8_t *src = (const uint8_t *)data;
		for (size_t i = 0; i < len; i++) dst[i] = src[i];
		b->len += len;
	}
	return true;
}

void vivi_buf_release(vivi_buf *b)
{
	b->p = nullptr;
	b->len = 0;
	b->cap = 0;
}

void vivi_buf_free(vivi_buf *b)
{
	vivi_dealloc(b->p);
	b->p = nullptr;
	b->len = 0;
	b->cap = 0;
}

void vivi_rng_init(vivi_rng *r, uint32_t seed)
{
	r->state = seed ? seed : 1u;
}

uint32_t vivi_rng_next(vivi_rng *r)
{
	/* xorshift32 (13,17,5) -- identical constants to the Lua framework */
	uint32_t x = r->state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	r->state = x;
	return x;
}

void vivi_prng_init(vivi_prng *p, uint64_t seed)
{
	p->state = seed;
}

static uint64_t prng_mix(vivi_prng *p)
{
	uint64_t z = (p->state += 0x9E37'79B9'7F4A'7C15ull);
	z = (z ^ (z >> 30)) * 0xBF58'476D'1CE4'E5B9ull;
	z = (z ^ (z >> 27)) * 0x94D0'49BB'1331'11EBull;
	return z ^ (z >> 31);
}

uint64_t vivi_prng_next64(vivi_prng *p) { return prng_mix(p); }

uint32_t vivi_prng_next(vivi_prng *p) { return (uint32_t)(prng_mix(p) >> 32); }

bool vivi_prng_chance(vivi_prng *p, double prob)
{
	if (prob <= 0.0) return false;
	if (prob >= 1.0) return true;
	uint64_t thr = (uint64_t)(prob * 18446744073709551616.0);
	return prng_mix(p) < thr;
}




