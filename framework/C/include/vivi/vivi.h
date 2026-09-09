/* vivi.h -- common utilities for the genetic storage framework (C23 port).
 *
 * The framework stores everything as packed strands: plain byte arrays
 * with explicit lengths (never NUL-terminated). vivi_bytes is the
 * universal view/owner type: functions that produce data allocate it
 * through the portability layer below, and the caller frees it with
 * vivi_bytes_free().
 *
 * PORTABILITY. The library is environment-agnostic: it never touches
 * the OS, files, console, clock or threads, and needs no hosted CRT.
 * It only asks for raw memory through this hook pair:
 *
 *   - hosted builds (default): malloc/free;
 *   - freestanding / bare-metal: compile with -DVIVI_NO_HOSTED and call
 *     vivi_set_allocator() once before any other vivi_* call, wiring it
 *     to the target's allocator or a static bump pool.
 *
 * Beyond allocation the library needs only memcpy/memset/memcmp, which
 * every freestanding toolchain is required to provide. There is no I/O
 * of any kind: writing a strand to a file, a terminal or a bus is the
 * embedding user's business.
 */
#ifndef VIVI_H
#define VIVI_H

#include <stddef.h>
#include <stdint.h>

typedef struct vivi_bytes {
	uint8_t *data;
	size_t len;
} vivi_bytes;

/* --- portability surface: memory --------------------------------------- */
typedef void *(*vivi_alloc_fn)(size_t);
typedef void (*vivi_free_fn)(void *);

void vivi_set_allocator(vivi_alloc_fn alloc_fn, vivi_free_fn free_fn);

void *vivi_alloc(size_t size);            /* raw memory, may return nullptr */
void *vivi_zalloc(size_t n, size_t size); /* calloc semantics: n*size, zeroed */
void vivi_dealloc(void *p);               /* nullptr is a no-op */

static inline void vivi_bytes_free(vivi_bytes *b)
{
	if (b) {
		vivi_dealloc(b->data);
		b->data = nullptr;
		b->len = 0;
	}
}

static inline void vivi_bytes_free_n(vivi_bytes *b, size_t n)
{
	for (size_t i = 0; i < n; i++)
		vivi_bytes_free(&b[i]);
	vivi_dealloc(b);
}

/* growable byte buffer used internally while building strands */
typedef struct {
	uint8_t *p;
	size_t len, cap;
} vivi_buf;

[[nodiscard]] bool vivi_buf_reserve(vivi_buf *b, size_t need);
[[nodiscard]] bool vivi_buf_push(vivi_buf *b, uint8_t c);
[[nodiscard]] bool vivi_buf_append(vivi_buf *b, const void *data, size_t len);
void vivi_buf_release(vivi_buf *b);   /* moves ownership out */
void vivi_buf_free(vivi_buf *b);

/* deterministic xorshift32 rng (same constants as the Lua framework) */
typedef struct { uint32_t state; } vivi_rng;

void vivi_rng_init(vivi_rng *r, uint32_t seed);
uint32_t vivi_rng_next(vivi_rng *r);

static_assert(sizeof(uint32_t) * 8 == 32, "the codec requires 32-bit uint32_t");
static_assert(sizeof(size_t) >= sizeof(uint32_t), "size_t must hold 32-bit lengths");

#endif /* VIVI_H */

