# Usage contracts

This is the fine print: ownership, error handling, limits, determinism and
the rules that keep the genome model consistent.

## Memory

### Allocator

- Hosted builds default to `malloc`/`free`. Freestanding builds
  (`-DVIVI_NO_HOSTED`) start with **no** allocator: every call fails until
  `vivi_set_allocator()` is called once, before any other `vivi_*` call.
- The allocator pair is global process state. Do not swap it while any
  library allocation is still alive; the free side must match the alloc side.
- `vivi_zalloc(n, size)` checks the `n * size` overflow for you.
- The library is otherwise environment-free: no files, time, console or
  threads; only `memcpy`, `memset`, `memcmp` and the allocator hooks.

### Ownership matrix

| Producer | You must free with |
|---|---|
| any function filling a `vivi_bytes *` | `vivi_bytes_free()` |
| `organism_read` (`vivi_bytes **out`) | `vivi_bytes_free_n(*out, o->nchr)` |
| `chr_parse` → `chr_record` | `chr_record_free(&rec)` (frees genes + data) |
| `genome_gene_scan` result | `genome_scan_free(&res)` |
| `genome_gene_read` (`genome_gene *out`) | `vivi_dealloc(out->data)` (ownership moves to you) |
| `org_mut_result.data` | `vivi_bytes_free(&result.data)` |
| `vivi_buf_release` | the moved-out `b->p` (save `p`/`len` first) |

Borrowed, never freed by you:

- `c->stem_source` / `o->stem_source`: the stem must outlive every user it
  is attached to; `*_attach_stem(x, NULL)` clears the link.
- `rec` and `g` structs passed *into* functions.
- `*err` message strings (static).

### Output contract

Unless documented otherwise, a function that fails either leaves `*out`
untouched or zeroed (`{ NULL, 0 }`). Do not read an output after a `false`
return.

## Errors

- `err` may be `nullptr`; the `bool` return is authoritative.
- `*err` holds a static string and is only meaningful when the call failed.
  Never free it, and don't treat it as stable across calls.
- `[[nodiscard]]` marks calls whose result you should check.

## Limits and validation

| Item | Allowed values |
|---|---|
| payload length (`dna_encode`) | 0 .. 0xFFFFFFFF bytes |
| `dna_opts.h` | 0 (=3) or 1..12 |
| `dna_opts.gc_eps` | 0 (=0.05) or [0.005, 0.5) |
| gene raw / packed length | 0 .. 65535 bytes |
| gene `id`, `usertype`, chromosome `id`, `flags` | 0 .. 255 |
| gene `h` (dense mode) | 3 .. 12 |
| `chr_opts.gene_raw` | 0 (=1024) or >= 16 |
| `chr_opts.units` | 0 (=4) or >= 1 |
| genes per chromosome | <= 256 (one-byte ids) |
| chromosomes per organism | 1 .. 255 |
| `generation` | 0 .. 65535 (centromere field) |
| `max_gen` | 0 = 60; replication at `generation + 1 > max_gen` fails "senescent" |

## Determinism and byte format

- All randomness is xorshift32 with explicit seeds; the same seed and inputs
  produce the same bytes on every platform (byte order is fixed, no padding
  in the formats). `seed == 0` is normalized to 1.
- `dna_encode` guarantees homopolymer runs <= `h` and GC content inside
  `[0.5-eps, 0.5+eps]`; `dna_decode` reconstructs the payload exactly.
- Gene tags are 32-bit Chaskey-12 over the raw data; a payload change is
  detected with probability ~1 - 2^-32. Wobble-base corruption is invisible
  by design. The key is public: integrity, not authenticity.
- Damage/mutation never insert or delete: strand lengths are stable, so
  homologous loci stay aligned.
- The C port is byte-compatible with the Lua reference (`framework/luau/`)
  for codec output, gene/chromosome encoding, damage patterns, cell
  homologs and the `VIV1` container.

## Concurrency

There is none. The codec fast tables (`dna_free_caches()` clears them) and
the allocator hooks are process-global mutable state. Use one lock around
the library, or one instance per thread with a thread-local allocator, if
you need concurrency.

## Model semantics worth knowing

- **Diploid defaults are homozygous.** `organism_new`, `cell_new`,
  `organism_mutate`, `organism_cross` and `organism_deserialize` all create
  identical homologs; divergence only appears after damage to one copy.
- **Checkpoint repairs in place** and reports events:
  `repaired` (gene spliced from the homolog), `structural` (telomere or
  centromere splice), `dead` (locus broken in both homologs),
  `anomaly` (extra/phantom genes).
- **Stem renewal resets generation to 0** and replaces the entire genome
  (the stem may even have a different chromosome count). `renewed` is set
  on the report by `read` / `maintain`.
- **Mitosis ages the mother too.** Both parent and daughter end at the new
  generation, and the daughter inherits the stem link. `replicate` is the
  self-only variant.
- **`organism_cross` is not meiosis of the parents**: the child is a fresh
  per-gene mosaic, homozygous, generation 0, `max_gen` from parent A, no
  stem link. Parents are untouched. The parents must have equal `nchr`,
  matching ids and matching chromosome options.
- **`organism_read` may resize the organism** if auto-renewal from a
  differently sized stem happens; always free the result with the current
  `o->nchr`.

## Serialization details

- `VIV1` is deterministic and stores only homolog 0. Loading makes the
  organism homozygous and resets `max_gen` to 60. Trailing bytes are
  ignored; duplicate/unsorted chromosome ids are not rejected on load.
- The Lua and C implementations accept each other's containers.

## Freestanding checklist

1. Compile with `-DVIVI_NO_HOSTED`.
2. Call `vivi_set_allocator(alloc, free)` before anything else.
3. Provide `memcpy`, `memset`, `memcmp` (standard requires them).
4. A symbol audit of a linked image shows only those three libc names
   (plus toolchain markers such as `_fltused` on MSVC) and `vivi_*`.
