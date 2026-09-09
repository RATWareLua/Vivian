/* vivi.c -- shared buffer, allocator hooks and rng helpers */
#include <string.h>
#include "vivi/vivi.h"

#ifdef VIVI_NO_HOSTED
static vivi_alloc_fn vivi_the_alloc = nullptr; /* must be set via vivi_set_allocator */
static vivi_free_fn vivi_the_free = nullptr;
#else
#include <stdlib.h>
static vivi_alloc_fn vivi_the_alloc = malloc;
static vivi_free_fn vivi_the_free = free;
#endif

void vivi_set_allocator(vivi_alloc_fn alloc_fn, vivi_free_fn free_fn)
{
	vivi_the_alloc = alloc_fn;
	vivi_the_free = free_fn;
}

void *vivi_alloc(size_t size)
{
	return vivi_the_alloc ? vivi_the_alloc(size) : nullptr;
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
	if (vivi_the_free) vivi_the_free(p);
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




