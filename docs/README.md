# Vivian C API — documentation

This folder documents the C implementation in `framework/C`. Start here,
then follow the two references:

1. [`api.md`](api.md) — every public type and function, module by module.
2. [`contracts.md`](contracts.md) — memory ownership, error handling, limits,
   determinism, threading, and the lifecycle rules you must respect.
3. [`research.md`](research.md) — the error channel, `vivi_sim` experiments,
   format evolution plan and the research roadmap.

The Lua implementation in `framework/luau/` is the normative byte-format
reference: for the same inputs the C port produces exactly the same bytes.

## Build and run

```bat
cd framework\C
build.bat lib        rem static library -> vivi.lib
build.bat test       rem 387-check suite
build.bat asan       rem the same suite under AddressSanitizer
build.bat density    rem packing-density report
```

`make lib | test | demo | density | asan` is equivalent on POSIX.
Public headers live in `include/vivi/`; include them as `#include <vivi/organism.h>`.

## Minimal tracked example

```c
#include <stdio.h>
#include <vivi/organism.h>

int main(void)
{
    const char *err = nullptr;
    const int ids[2] = { 0, 1 };
    const chr_opts co = { 1024, 0, 3, 4, 0 };   /* gene_raw, codon, h, units, flags */
    const chr_opts opts[2] = { co, co };
    static const char D0[] = "BRAIN:CURIOUS;SPEED:12;";
    static const char D1[] = "METABOLISM:FAST;STAMINA:100;";
    const uint8_t *datas[2] = { (const uint8_t *)D0, (const uint8_t *)D1 };
    const size_t lens[2] = { sizeof(D0) - 1, sizeof(D1) - 1 };

    vivi_organism *pet = nullptr;
    if (!organism_new(&pet, ids, opts, datas, lens, 2, 60, &err)) {
        printf("new failed: %s\n", err);
        return 1;
    }

    vivi_bytes *traits = nullptr;      /* organism_read allocates the array */
    cell_report rep;
    if (organism_read(&traits, pet, &rep, &err)) {
        printf("repair: repaired=%d structural=%d dead=%d\n",
            rep.repaired, rep.structural, rep.dead);
        vivi_bytes_free_n(traits, pet->nchr);   /* frees the entries AND the array */
    }

    (void)organism_damage(pet, 2, 0xC0FFEEu, &err);
    traits = nullptr;                  /* read replaces the pointer */
    if (organism_read(&traits, pet, &rep, &err)) {
        printf("after damage: renewed=%d repaired=%d\n", rep.renewed, rep.repaired);
        vivi_bytes_free_n(traits, pet->nchr);
    }

    organism_free(pet);
    dna_free_caches();                 /* optional: releases codec fast tables */
    return 0;
}
```

Compile from `framework/C`:

```bat
clang -std=c23 -O2 -Iinclude -o hello hello.c src\vivi.c src\dna.c src\genome.c ^
    src\chromosome.c src\cell.c src\organism.c
```

## Rules of thumb

- Every `vivi_bytes` you receive is heap memory: release it with
  `vivi_bytes_free()`; arrays of them with `vivi_bytes_free_n(array, count)`.
- `err` may be `nullptr` — failures are still reported by the `bool` return.
- No I/O, no threads, no hidden state beyond the codec caches; on bare metal
  call `vivi_set_allocator()` once before anything else.
- Damage never changes strand length, so loci stay aligned across homologs.
- `stem_source` is a borrowed pointer: the stem must outlive its users.
