# API reference

All public functions return `bool` (or a report struct) and take a trailing
`const char **err`. `err` may be `nullptr`; when it is not and a function
returns `false`, `*err` points to a static message string (never free it).

Returned `vivi_bytes` values own heap memory; release them with
`vivi_bytes_free()` / `vivi_bytes_free_n()` (see [contracts.md](contracts.md)).

---

## `vivi.h` — common utilities

```c
typedef struct { uint8_t *data; size_t len; } vivi_bytes;
typedef void *(*vivi_alloc_fn)(size_t);
typedef void (*vivi_free_fn)(void *);

void  vivi_set_allocator(vivi_alloc_fn alloc_fn, vivi_free_fn free_fn);
void *vivi_alloc(size_t size);             /* may return nullptr */
void *vivi_zalloc(size_t n, size_t size);  /* calloc semantics */
void  vivi_dealloc(void *p);               /* nullptr is a no-op */

void vivi_bytes_free(vivi_bytes *b);              /* one buffer */
void vivi_bytes_free_n(vivi_bytes *b, size_t n);  /* n buffers + the array */
```

`vivi_bytes` is a plain pointer+length view; it is never NUL-terminated.

### Growable buffer (for embedders)

```c
typedef struct { uint8_t *p; size_t len, cap; } vivi_buf;
bool vivi_buf_reserve(vivi_buf *b, size_t need);
bool vivi_buf_push(vivi_buf *b, uint8_t c);
bool vivi_buf_append(vivi_buf *b, const void *data, size_t len);
void vivi_buf_release(vivi_buf *b);  /* hands ownership of b->p to the caller */
void vivi_buf_free(vivi_buf *b);     /* frees b->p */
```

### Deterministic RNG

```c
typedef struct { uint32_t state; } vivi_rng;
void     vivi_rng_init(vivi_rng *r, uint32_t seed);  /* seed 0 is normalized to 1 */
uint32_t vivi_rng_next(vivi_rng *r);                 /* xorshift32 (13,17,5) */
```

The library is freestanding (`-DVIVI_NO_HOSTED`): it needs only `memcpy`,
`memset`, `memcmp` and the allocator hooks.

---

## `dna.h` — quaternary codec

Binary data <-> packed A/C/G/T strand: 2 bits per base, 4 bases per byte.
The strand is self-describing (whitened 32-bit big-endian length header,
constraint-coded payload, GC fill).

```c
typedef struct {
    int h;          /* homopolymer limit: 0 = default 3, else 1..12 */
    double gc_eps;  /* GC half-band: 0 = default 0.05, else [0.005, 0.5) */
} dna_opts;

bool dna_encode(vivi_bytes *out, const uint8_t *data, size_t len,
                const dna_opts *opts, const char **err);   /* opts may be NULL */
bool dna_decode(vivi_bytes *out, const uint8_t *strand, size_t slen,
                const char **err);
```

`dna_decode` validates the strand prefix, homopolymer transitions and the
embedded length; it returns the exact original payload.

```c
double dna_gc_content(const uint8_t *s, size_t n);     /* whole strand */
size_t dna_max_homopolymer(const uint8_t *s, size_t n);

bool dna_complement(vivi_bytes *out, const uint8_t *s, size_t n, const char **err);
bool dna_reverse_complement(vivi_bytes *out, const uint8_t *s, size_t n, const char **err);

bool dna_to_ascii(vivi_bytes *out, const uint8_t *s, size_t n, long nbases,
                  const char **err);
bool dna_from_ascii(vivi_bytes *out, const char *s, size_t n, const char **err);
bool dna_validate(const uint8_t *s, size_t n, int h, double eps, const char **err);

void dna_free_caches(void);
```

- `dna_to_ascii`: `nbases < 0` or larger than available = the whole strand;
  a partial trailing byte is rendered as its first `nbases % 4` bases.
- `dna_from_ascii`: whitespace is skipped, lowercase accepted; a trailing
  partial group is left-aligned (matches the Lua reference).
- `dna_validate`: constraint check only (homopolymer <= `h`, GC inside
  `[0.5-eps, 0.5+eps]`), not a decode.
- `dna_free_caches`: releases the lazily built per-`h` fast tables; pure
  cache, rebuilt on demand. Optional housekeeping for leak checks.

---

## `genome.h` — genes and Chaskey-12

### Integrity tag

```c
void     genome_chaskey_tag(uint32_t out[4], const uint8_t *s, size_t len,
                            const uint32_t key[4]);   /* key NULL = public default */
uint32_t genome_chaskey32(const uint8_t *s, size_t len);
```

Chaskey-12 (ISO/IEC 29192-6), tag truncated to 32 bits. The default key is
public: this is corruption detection, not authentication.

### Codon helpers

```c
bool genome_pack_codons(vivi_bytes *out, const int *values, size_t nvalues,
                        const char **err);          /* nvalues % 4 == 0, values 0..15 */
bool genome_unpack_codons(int *values, size_t nvalues, const uint8_t *packed,
                          size_t plen, const char **err);
bool genome_values_from_bytes(int *values, size_t nvalues, const uint8_t *s,
                              size_t len);          /* needs nvalues >= 2*len */
bool genome_bytes_from_values(vivi_bytes *out, const int *values, size_t nvalues,
                              size_t nbytes, const char **err);  /* needs nvalues >= 2*nbytes */
```

A codon is 3 bases `(b1, b2, wobble)`; the value is `b1*4 + b2` (4 bits) and
the wobble base carries no data, so its corruption is tolerated.

### Gene records

```c
bool genome_gene_encode(vivi_bytes *out, int id, int usertype, int codon_mode,
                        int h, const uint8_t *data, size_t len, const char **err);
```

Layout: `[PROM 3][header 9][payload][tag 6][TERM 3]`. `id` and `usertype`
are bytes (0..255); `usertype`'s low bit is reserved for `codon_mode`.
`len <= 65535`; dense mode needs `h` in 3..12. The payload is either the
DNA codec output (dense) or the raw data re-spelled as codons (codon mode).

```c
typedef struct {
    int id, type;
    int codon;          /* 0 = dense, 1 = codon */
    int rawlen, packedlen;
    size_t offset;      /* promoter index in the scanned strand */
    size_t size;        /* total gene bytes: 21 + packedlen */
    int crc_ok;         /* tag verified AND payload decodable */
    uint32_t tag;       /* tag as decoded from the strand */
    uint8_t *data;      /* decoded raw data (NULL unless crc_ok) */
    size_t data_len;
} genome_gene;

typedef struct { genome_gene *genes; size_t count; } genome_scan_result;

void genome_scan_free(genome_scan_result *r);
bool genome_gene_scan(genome_scan_result *out, const uint8_t *strand,
                      size_t slen, const char **err);
bool genome_gene_read(genome_gene *out, const uint8_t *strand, size_t slen,
                      int id, const char **err);   /* id < 0 = first valid gene */
void genome_set_base(uint8_t *strand, size_t n, size_t idx, int digit);
```

- `genome_gene_scan` walks the strand and collects every well-formed gene
  (promoter to terminator); each entry's tag check is recorded in `crc_ok`.
  Free the result with `genome_scan_free()`.
- `genome_gene_read` returns the first valid gene with the requested id and
  **moves ownership** of `out->data` to the caller (release with
  `vivi_dealloc(out->data)`).
- `genome_set_base` rewrites one base (0-based index, digit 0..3) in place;
  out-of-range indices are ignored.

---

## `chromosome.h` — chromosome layer

Layout: `[telomere][centromere][gene 0]..[gene n-1][telomere]`.
`CENBYTES` (18) is the centromere size: marker 3 + 20 codons 15.

```c
typedef struct {
    int gene_raw;   /* 0 = default 1024, else >= 16: raw bytes per gene */
    int codon;      /* 0 = dense payloads, 1 = codon payloads */
    int h;          /* dense gene homopolymer limit, 0 = default 3 */
    int units;      /* telomere repeat units, 0 = default 4 (>= 1) */
    int flags;      /* user byte, stored in the centromere */
} chr_opts;

typedef struct {
    int id, flags, ngenes, generation;
    int telomere_ok, cen_ok;
    size_t telomere_bytes;
    genome_gene *genes;    /* owned; release with chr_record_free */
    size_t gene_count;
} chr_record;
```

```c
bool chr_encode(vivi_bytes *out, int chr_id, const uint8_t *data, size_t len,
                const chr_opts *opts, const char **err);  /* opts may be NULL */
bool chr_parse(chr_record *out, const uint8_t *strand, size_t slen,
               int units_hint, const char **err);          /* < 0 = autodetect */
void chr_record_free(chr_record *r);
bool chr_read(vivi_bytes *out, const chr_record *rec, const char **err);
bool chr_set_generation(vivi_bytes *out, const uint8_t *strand, size_t slen,
                        int generation, const char **err);  /* 0..65535 */
```

- The number of genes is bounded by the one-byte gene id: **256 per
  chromosome** (`data` larger than `256 * gene_raw` fails to encode).
- `chr_parse` finds genes by scanning; `genes[i].offset` is absolute in the
  strand. A damaged centromere leaves `cen_ok == 0` and no genes.
- `chr_read` is strict: requires a valid centromere, gene count match, no
  duplicate/missing ids and every tag valid; then concatenates the payloads.
- `chr_set_generation` re-stamps the centromere generation, preserving the
  strand length byte for byte.

---

## `cell.h` — diploid cell and repair

```c
typedef struct {
    int repaired;    /* genes spliced from the healthy homolog */
    int dead;        /* genes broken in BOTH homologs (or a lost centromere) */
    int structural;  /* telomere/centromere splices */
    int anomaly;     /* extra/phantom genes seen during repair */
    int renewed;     /* set by read()/maintain() when the stem rescued the cell */
} cell_report;

typedef struct vivi_cell {
    int chr_id, generation, max_gen;
    int stem, dead;
    const struct vivi_cell *stem_source;   /* borrowed, may be NULL */
    uint8_t *hom[2];
    size_t hlen[2];
} vivi_cell;
```

### Low-level machinery

```c
bool cell_damage_strand(vivi_bytes *out, const uint8_t *strand, size_t slen,
                        int count, uint32_t seed, const char **err);
bool cell_repair_homolog(vivi_bytes *out, const uint8_t *dst, size_t dlen,
                         const uint8_t *src, size_t slen, cell_report *rep,
                         const char **err);
bool cell_checkpoint_pair(vivi_bytes *a_out, vivi_bytes *b_out,
                          const uint8_t *h1, size_t l1, const uint8_t *h2, size_t l2,
                          cell_report *rep, const char **err);
```

- `cell_damage_strand` flips `count` random digits (never length); `count <= 0`
  returns an unchanged copy; deterministic per `seed`. Empty input yields a
  0-byte result.
- `cell_repair_homolog` repairs `dst` using `src` as the healthy template
  (equal lengths required). `cell_checkpoint_pair` is the three-pass
  symmetric cascade: each homolog repairs the other, then the first gets a
  third pass.

### Lifecycle

```c
bool cell_new(vivi_cell **out, int chr_id, const uint8_t *data, size_t len,
              const chr_opts *opts, const char **err);
bool cell_stem(vivi_cell **out, int chr_id, const uint8_t *data, size_t len,
               const chr_opts *opts, const char **err);
void cell_attach_stem(vivi_cell *c, const vivi_cell *stem);  /* NULL clears */
bool cell_renew(vivi_cell *c, const vivi_cell *stem, const char **err);
void cell_kill(vivi_cell *c);
void cell_free(vivi_cell *c);

bool cell_read(vivi_bytes *out, vivi_cell *c, cell_report *rep, const char **err);
cell_report cell_checkpoint(vivi_cell *c);
bool cell_replicate(vivi_cell *c, const char **err);
bool cell_mitosis(vivi_cell **out, vivi_cell *c, const char **err);
bool cell_damage(vivi_cell *c, int count, uint32_t seed, const char **err);
cell_report cell_maintain(vivi_cell **cells, size_t count, const vivi_cell *stem);
```

- `cell_new` starts both homologs as copies of the same strand (generation 0,
  `max_gen` = 60); `cell_stem` additionally marks the cell as a stem.
- `cell_read` checkpoints (repairs in place), then returns the payload from
  the first readable homolog. If the cell is dead or unreadable and a stem is
  attached, it renews and reports `renewed = 1`.
- `cell_replicate` / `cell_mitosis` checkpoint, reject when past `max_gen`
  (`"senescent"`), advance the generation and re-stamp both centromeres.
  Mitosis also ages the mother (both cells end at the new generation) and
  the daughter inherits the stem link.
- `cell_maintain` runs a population pass: checkpoint, replace unreadable
  cells from the stem, kill what cannot be rescued.

---

## `organism.h` — organism layer

```c
typedef struct vivi_organism {
    int *chr_ids;
    chr_opts *chr_opts;
    size_t nchr;
    uint8_t **hom[2];   /* per homolog: nchr strands */
    size_t *hlen[2];
    int generation, max_gen;
    int stem, dead;
    const struct vivi_organism *stem_source;  /* borrowed */
} vivi_organism;

typedef struct {
    int chr;            /* chromosome id that was mutated */
    vivi_bytes data;     /* mutated raw data (owned) */
} org_mut_result;
```

### Construction and lifetime

```c
bool organism_new(vivi_organism **out, const int *ids, const chr_opts *opts,
                  const uint8_t *const *datas, const size_t *lens, size_t nchr,
                  int max_gen, const char **err);      /* max_gen 0 = 60 */
bool organism_stem(vivi_organism **out, ...same...);
void organism_attach_stem(vivi_organism *o, const vivi_organism *stem);
bool organism_renew(vivi_organism *o, const vivi_organism *stem, const char **err);
void organism_kill(vivi_organism *o);
void organism_free(vivi_organism *o);
```

- `nchr` must be 1..255 and ids unique bytes; every chromosome is stored
  twice (homologs).
- `organism_renew` replaces the whole genome with the stem's and resets the
  generation to 0. The stem need not have the same chromosome count.
- `organism_kill` marks dead and releases the strands; `organism_free` also
  releases the organism struct.

### Evolution and reading

```c
bool organism_read(vivi_bytes **out, vivi_organism *o, cell_report *rep,
                   const char **err);
cell_report organism_checkpoint(vivi_organism *o);
bool organism_replicate(vivi_organism *o, const char **err);
bool organism_mitosis(vivi_organism **out, vivi_organism *o, const char **err);
bool organism_damage(vivi_organism *o, int count, uint32_t seed, const char **err);
bool organism_mutate(org_mut_result *out, vivi_organism *o, int count,
                     uint32_t seed, const char **err);
bool organism_cross(vivi_organism **out, vivi_organism *pa, vivi_organism *pb,
                    uint32_t seed, const char **err);
```

- `organism_read` checkpoints every chromosome, then returns an array of
  `o->nchr` payloads (`out[i]` belongs to `chr_ids[i]`). On success the array
  is freshly allocated — free with `vivi_bytes_free_n(*out, o->nchr)`.
  If a stem is attached, an unreadable organism is renewed automatically and
  `rep->renewed` is set.
- `organism_damage` writes deterministic damage to both homologs of every
  chromosome (per-chromosome different streams).
- `organism_mutate` picks one chromosome, flips `count` random bits in its
  raw data, re-encodes it and **replaces both homologs with the new strand**
  (the mutation becomes homozygous). `out->data` is the mutated payload and
  is owned by the caller. `count <= 0` or an empty payload returns an
  unchanged copy.
- `organism_cross` requires equal `nchr`, matching ids and matching per-
  chromosome options; it builds a per-gene mosaic from the first valid copy
  of each gene in either parent. The child is homozygous, starts at
  generation 0, inherits `max_gen` from `pa`, and has no stem link.
- `organism_replicate` / `organism_mitosis` advance and re-stamp the
  generation; past `max_gen` they fail with `"senescent"`. The daughter in
  mitosis inherits the stem link.

### Serialization

```c
bool organism_serialize(vivi_bytes *out, vivi_organism *o, const char **err);
bool organism_deserialize(vivi_organism **out, const uint8_t *s, size_t len,
                          const char **err);
```

The `VIV1` container:

```
"VIV1" | generation (2 BE) | nchr (1)
per chromosome:
  id (1) | gene_raw (2 BE) | mode|h (1) | units (1) | flags (1)
  | strand length (4 BE) | strand bytes
```

Deterministic: the same organism serializes to the same bytes. Only homolog
0 is stored; deserialization duplicates it into both homologs (a
homozygous organism). `max_gen` is not stored and resets to 60. Trailing
bytes after the last chromosome are ignored.

### Population maintenance

```c
cell_report organism_maintain(vivi_organism **organisms, size_t count,
                              const vivi_organism *stem);
```

One pass over the population: checkpoint, revive from the stem (also when
`generation >= max_gen`), kill what cannot be rescued.

---

## `channel.h` — synthesis/sequencing channel (research layer)

```c
typedef struct {
    double p_sub;   /* substitution probability per original base, [0, 1] */
    double p_ins;   /* insertion probability after each original base, [0, 1] */
    double p_del;   /* deletion probability per original base, [0, 1] */
    double p_drop;  /* probability that the whole read is lost, [0, 1] */
    uint32_t seed;  /* rng seed; 0 behaves as 1 */
} vivi_channel_opts;

typedef struct {
    vivi_bytes strand;  /* damaged packed strand ({ NULL, 0 } when dropped) */
    size_t bases;       /* exact base count of the read */
    int dropped;        /* 1 = the read was lost */
} vivi_read;

[[nodiscard]] bool vivi_channel_read(vivi_read *out, const uint8_t *strand, size_t slen,
    const vivi_channel_opts *opts, const char **err);
```

- `opts` may be `NULL` (identity channel); probabilities outside [0, 1] are
  rejected.
- Draw order per original base is pinned (delete?, substitute?, insert?), so
  a read is a deterministic function of `(strand, rates, seed)`.
- A dropped read comes back as `dropped = 1` with `strand = { NULL, 0 }`.
- The damaged sequence is repacked into whole bytes; `bases` is the exact
  base count (0..3 filler bases are appended as `A`).

See [research.md](research.md) for the `vivi_sim` experiment tool, the CSV
schema and measured success curves.
