# Vivian — data that lives

> Encode any bytes into a diploid synthetic genome: homologous repair,
> viable point mutations, crossing-over, and a stem-cell niche for
> rejuvenation. Not a checksum container — an organism with a medical chart.

**Vivian** is a portable C23 framework that stores data the way life does.
Your bytes become the raw payload of a diploid genome: every piece of data
is written twice (as two homologous chromosomes), damage is repaired by
splicing the healthy homolog, integrity is enforced by a Chaskey-12 tag
living in wobbling codons, and the whole thing ages, dies, mutates and
breeds like a population of tiny organisms.

## What it is

Three layers, one idea:

1. **The codec (biochemistry).** Any byte stream is re-spelled in a
   4-letter alphabet (A/C/G/T, 2 bits per base) with biologically
   motivated constraints: bounded homopolymer runs and a target GC
   content. No 4-bit tail is ever lost to padding.
2. **The architecture (anatomy).** Data is wrapped into genes —
   `[PROM][header][payload][Chaskey-12 tag as wobble codons][TERM]` —
   and genes into chromosomes: `[TELOMERE][CENTROMERE][genes][TELOMERE]`.
   Cells are diploid: two homologs of every chromosome, a three-pass
   checkpoint on every replication, senescence at the Hayflick limit.
3. **The model (a way of life).** Organisms can take damage, get healed,
   permanently die, be resurrected from a stem-cell niche, mutate with a
   single honest bit-flip, and reproduce: `cross()` splices two parents
   into a mosaic child gene by gene.

So when something goes wrong, you don't get *"file corrupt"* — you get a
report: `Repaired: 1, Structural: 2, Dead: 1`, plus a genome that
physically survived.

## What it is for

- **Learning & demonstration.** DNA data storage, constrained encoding,
  MAC-tagged self-describing records, diploid redundancy as an erasure
  code — all of it visible and pokeable. A point mutation is not an
  abstract CRC failure: `CARNIVORE` literally turns into `CARNIVNRE`.
- **An experiment in failure optics.** Checksums say binary yes/no.
  Vivian treats damage as a medical event with history: how much was
  repaired, what became structural, what died, whether the stem niche
  rescued the line.
- **Generative pets & art.** The genome *is* the data. Evolve it, breed
  it with wild mates, serialize it, load it years later — still one
  self-contained strand format.

## What it is not

- Not encryption: the Chaskey-12 tag is an integrity MAC and its key is
  public.
- Not a compressor: expect a small, deterministic size overhead (the
  VIV1 container adds ~3%).
- Not production archival software — it is a rigorous toy: 379-check test
  suite, official Chaskey-12 vectors, ASan-clean, but experimental.

## Quick start

C23, clang, no dependencies beyond the CRT's `memcpy/memset/memcmp`:

```bat
cd C
build.bat lib        rem framework only -> vivi.lib
build.bat            rem + test suite + interactive demo
vivi_test.exe        rem pass=379 fail=0
```

Or with `make`: `make lib`, `make test`, `make demo`.

### Hello, organism

```c
#include <stdio.h>
#include <vivi/organism.h>

int main(void)
{
    const char *err = nullptr;
    const int ids[2] = { 0, 1 };
    const chr_opts co[2] = { { 1024, 0, 3, 4, 0 }, { 1024, 0, 3, 4, 0 } };
    static const char D0[] = "BRAIN:CURIOUS;SPEED:12;";
    static const char D1[] = "METABOLISM:FAST;STAMINA:100;";
    const uint8_t *data[2] = { (const uint8_t *)D0, (const uint8_t *)D1 };
    const size_t lens[2] = { sizeof(D0) - 1, sizeof(D1) - 1 };

    vivi_organism *pet = nullptr;
    if (!organism_new(&pet, ids, co, data, lens, 2, 60, &err)) return 1;

    /* diploid read: every chromosome goes through homologous repair */
    vivi_bytes *traits = nullptr;
    cell_report rep;
    if (organism_read(&traits, pet, &rep, &err)) {
        printf("repaired=%d structural=%d dead=%d\n",
            rep.repaired, rep.structural, rep.dead);
        vivi_bytes_free_n(traits, pet->nchr);
    }

    /* take some radiation, read again — the homologs fix it */
    (void)organism_damage(pet, 2, 0xC0FFEEu, &err);
    traits = nullptr;
    if (organism_read(&traits, pet, &rep, &err)) {
        printf("after damage: Chr %d traits intact\n", rep.repaired);
        vivi_bytes_free_n(traits, pet->nchr);
    }

    organism_free(pet);
    return 0;
}
```

Compile: `clang -std=c23 -O2 -Iinclude -o hello hello.c src\vivi.c src\dna.c
src\genome.c src\chromosome.c src\cell.c src\organism.c`

## Bare metal / freestanding

The library is environment-agnostic: no OS calls, no files, no console,
no clock, no threads, no stdio. It only asks for raw memory through
`vivi_set_allocator()` (hosted builds default to `malloc/free`). For
freestanding targets:

```c
vivi_set_allocator(my_bump_alloc, my_noop_free);   /* once, first thing */
```

```bat
clang -std=c23 -O2 -DVIVI_NO_HOSTED -Iinclude -c src\*.c
```

The symbol audit (`llvm-nm --undefined-only`) of a freestanding build
shows exactly three external symbols: `memcpy`, `memset`, `memcmp` — all
required of every freestanding toolchain by the standard. Everything else
is `vivi_*`. Writing a genome to a file, a terminal or a bus is the
embedding user's business; `demo/` and `test/` are hosted examples with
their own entrypoints.

## The model, term by term

| DNA concept | What it does for your data |
|---|---|
| bases A/C/G/T | 2 bits per base, bytes re-spelled, constraint-enforced |
| gene: PROM / TERM | self-describing record framing, findable by scan |
| wobble codons | Chaskey-12 integrity tag hidden in synonymous codons |
| telomeres | end markers that protect and locate the strand edges |
| centromere | codon-coded chromosome header (id, flags, gene count, generation) |
| diploid cell | every strand stored twice on two homologs |
| homologous repair | damage on one copy is spliced out from the other |
| checkpoint (3 passes) | replication verified three ways before it counts |
| Hayflick limit | generation counter; past it, no more replication |
| stem-cell niche | your snapshot: whole lines regrown from it |
| point mutation | a real bit-flip in payload, re-encoded with a fresh tag |
| crossing-over | two parents → one per-gene mosaic child |
| VIV1 container | versioned serialization of the whole organism |

## The VIV1 container

`"VIV1"` magic, generation (2 bytes BE), chromosome count (1 byte), then
per chromosome: id (1), gene_raw (2 BE), mode\|h (1), units (1), flags (1),
strand length (4 BE), strand bytes. Deterministic: same organism, same
bytes. Byte-compatible with the Lua reference below.

## Layout

```
.
├── framework/          Lua reference implementation (the "biochemistry" spec,
│                       byte-compatible with the C port)
├── C/
│   ├── include/vivi/   public API: vivi.h, dna.h, genome.h, chromosome.h,
│   │                   cell.h, organism.h
│   ├── src/            the framework itself (no I/O, no CRT assumptions)
│   ├── test/           379-check self-test (own entrypoint, hosted)
│   ├── demo/           interactive terminal tamagotchi (own entrypoint, hosted)
│   ├── vivi.lib        built with build.bat lib / make lib
│   ├── build.bat       all | lib | test | clean
│   └── Makefile        all | lib | test | demo | clean
```

## Verification

- **379/379** checks: 64 official Chaskey-12 vectors; codec round-trips
  for every parameter combination; constraint edge cases; gene/chromosome
  structure and corruption detection; diploid repair, checkpoints,
  senescence, stem rejuvenation; mutations, crossing-over mosaics; VIV1
  round-trips.
- AddressSanitizer-clean (library and demo), zero leaks (CRT leak check).
- Byte-format equality between the Lua reference and the C port.

## Status

Experimental, but honest: the format is pinned by tests and two
independent implementations agree on every byte. Expect the API to
settle, not the fun to.

---
*Что это по-русски:* C23-фреймворк, который хранит любые данные как диплоидный
синтетический геном: нуклеотидный кодек, гены с тегом Chaskey-12, хромосомы
с теломерами и центромерой, гомологичная репарация повреждений, мутации
реальными бит-флипами, кроссинговер и стволовое омоложение. Не формат
с чек-суммой, а организм с медкартой.
