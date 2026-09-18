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

/* --- machine-readable status ------------------------------------------- */
/* A stable classification of the library's static error strings, so a
 * caller can branch on the kind of failure instead of matching text. */
typedef enum vivi_errc {
	VIVI_OK = 0,
	VIVI_ERR_ARG,           /* invalid argument from the caller */
	VIVI_ERR_OOM,           /* allocation failed */
	VIVI_ERR_FORMAT,        /* malformed input / bad container */
	VIVI_ERR_CORRUPT,       /* data damaged beyond recovery */
	VIVI_ERR_LIMIT,         /* a documented limit was exceeded */
	VIVI_ERR_RANGE,         /* value outside its allowed range */
	VIVI_ERR_SENESCENT,     /* Hayflick limit reached */
	VIVI_ERR_INCOMPATIBLE,  /* mismatched genomes or options */
	VIVI_ERR_DEAD,          /* organism or cell is dead */
	VIVI_ERR_OTHER
} vivi_errc;

const char *vivi_strerror(vivi_errc code);
/* NULL/empty -> VIVI_OK; otherwise a stable code for the message text. */
vivi_errc vivi_error_code(const char *err);

/* --- stable API -------------------------------------------------------- */
/* VIVI_API_VERSION pins the stable surface documented in
 * docs/contracts.md. The research structures (vivi_organism, vivi_cell,
 * chr_record, genome_gene, ...) stay visible for advanced use, but their
 * layout is NOT part of that contract: stable code reads them through the
 * accessors in cell.h / organism.h and never dereferences fields. */
#define VIVI_API_VERSION_MAJOR 1
#define VIVI_API_VERSION_MINOR 0
#define VIVI_API_VERSION ((VIVI_API_VERSION_MAJOR << 16) | VIVI_API_VERSION_MINOR)
uint32_t vivi_api_version(void);

/* --- portability surface: memory --------------------------------------- */
typedef void *(*vivi_alloc_fn)(size_t);
typedef void (*vivi_free_fn)(void *);

/* Legacy global hooks: they configure the process-wide DEFAULT context,
 * used by threads that never call vivi_context_enter(). */
void vivi_set_allocator(vivi_alloc_fn alloc_fn, vivi_free_fn free_fn);

void *vivi_alloc(size_t size);            /* raw memory, may return nullptr */
void *vivi_zalloc(size_t n, size_t size); /* calloc semantics: n*size, zeroed */
void vivi_dealloc(void *p);               /* nullptr is a no-op */

/* --- execution context ------------------------------------------------- */
/* A context owns one allocator pair and the codec's lazy caches, so two
 * threads or two subsystems stay independent. Enter a per-thread context
 * before its allocations and keep it current until they are freed. */
typedef struct vivi_context vivi_context;

[[nodiscard]] vivi_context *vivi_context_new(vivi_alloc_fn alloc_fn, vivi_free_fn free_fn);
void vivi_context_free(vivi_context *ctx);

/* Set ctx as the current context for this thread; returns the previous
 * one to pass back to vivi_context_leave(). */
vivi_context *vivi_context_enter(vivi_context *ctx);
vivi_context *vivi_context_leave(vivi_context *prev);

/* Never NULL: this thread's context, or the process default. */
vivi_context *vivi_context_current(void);

/* Release the codec caches owned by ctx. */
void vivi_context_release_caches(vivi_context *ctx);

/* internal: the per-context DNA fast-table slots (13 entries) */
void **vivi_context_dna_cache(vivi_context *ctx);

/* thread-local storage qualifier: empty on a freestanding, single-threaded
 * build where _Thread_local need not exist */
#if defined(VIVI_NO_HOSTED)
#define VIVI_THREAD_LOCAL
#else
#define VIVI_THREAD_LOCAL _Thread_local
#endif

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

/* SplitMix64 -- the research-layer simulation PRNG (channel, pools).
 * vivi_rng above is frozen because damage patterns are byte-compatible
 * with the Lua reference; SplitMix64 has a 2^64 period and passes the
 * standard test batteries, which xorshift32 does not guarantee once a
 * sweep runs past 2^32 draws. */
typedef struct { uint64_t state; } vivi_prng;

void vivi_prng_init(vivi_prng *p, uint64_t seed);
uint32_t vivi_prng_next(vivi_prng *p);    /* upper 32 bits of one mix step */
uint64_t vivi_prng_next64(vivi_prng *p);
bool vivi_prng_chance(vivi_prng *p, double prob);  /* prob in [0, 1] */

static_assert(sizeof(uint32_t) * 8 == 32, "the codec requires 32-bit uint32_t");
static_assert(sizeof(size_t) >= sizeof(uint32_t), "size_t must hold 32-bit lengths");

#endif /* VIVI_H */

