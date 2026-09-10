# VIV14 — format revision design

`VIV14` is the first revision after `VIV1`. It keeps the codec, genes,
chromosomes and repair machinery untouched and modernizes two things that
VIV1 gets wrong at the container level:

1. **Versioned magic.** VIV1 identifies itself with exactly four bytes;
   future revisions need a version prefix and a growth path.
2. **Diploid persistence.** `organism_serialize` (VIV1) stores only homolog
   0 and duplicates it on load, so a damaged-but-not-yet-repaired homolog 1
   is silently lost. VIV14 stores both homologs and `max_gen`.

Everything below the container is unchanged, so VIV14 files are still
self-describing genomes with the same repair semantics.

## Revision ladder

The ladder is a promise, not a batch: each step changes the wire format
once, with its own migration.

| revision | scope |
|---|---|
| `VIV1` | current format; stays readable forever |
| `VIV14` | **this document**: versioned container, both homologs, `max_gen` |
| `VIV14N` | parity genes as first-class chromosome members (outer code in the format) — **implemented** |
| `VIV14NB4NSH33` | "Vivian Banshee": on-strand primer sites for physical random access; **implemented, layout frozen** |

Implementation policy from VIV14 on: **the C port is the reference**. The
Lua implementation in `framework/luau/` ports every revision (VIV14 and
VIV14N included) and is verified against the C output byte for byte in
CI (`make xcheck` / `framework/luau/xcheck.lua`).

## Wire format

```
"VIV14" (5)
flags (1)            bit0 = 1 (diploid pair stored); all other bits must be 0
generation (2 BE)
nchr (1)             1..255
max_gen (2 BE)       0 means 60 on load
per chromosome:
  id (1)             byte 0..255
  gene_raw (2 BE)
  mode|h (1)         bit0 codon, bits 1..4 = h-1
  units (1)
  flags (1)
  len0 (4 BE)  strand0 bytes
  len1 (4 BE)  strand1 bytes
```

Deterministic: the same organism and options produce the same bytes. The
per-chromosome header is byte-identical to VIV1's, so a VIV14 reader reuses
the same option validation.

## API

```c
int  organism_container_version(const uint8_t *s, size_t len);  /* 0, 1, 14, 141, 143 */

/* VIV1/VIV14/VIV14N/Banshee writers; organism_serialize() keeps VIV1 */
bool organism_serialize14(vivi_bytes *out, vivi_organism *o, const char **err);
bool organism_serialize14n(vivi_bytes *out, vivi_organism *o, const char **err);
bool organism_serialize14nb(vivi_bytes *out, vivi_organism *o, const char **err);

/* reads every revision */
bool organism_deserialize(vivi_organism **out, const uint8_t *s, size_t len,
                          const char **err);
```

## Validation and the magic ambiguity

`"VIV1"` is a prefix of `"VIV14"`, and VIV1's byte 5 is the generation high
byte, which can be `0x34` (`'4'`). `organism_deserialize` therefore tries
the VIV14 parser first and **falls back to VIV1** when the VIV14 structure
is invalid; a genuine VIV1 file with `'4'` there would have to pass VIV14
structural validation by coincidence.

VIV14 validation rules:

- `flags` must be exactly `1`;
- `max_gen` and `generation` are 16-bit; `nchr` in 1..255;
- per chromosome: `gene_raw >= 16`, `h` in 3..12, `units >= 1`, both strand
  buffers fully contained in the input;
- each homolog must pass `chr_parse` (structural), and **at least one**
  homolog must have a valid centromere with the matching id. The other may
  be damaged — that is the point of storing both.

## Migration

- VIV1 -> VIV14: `organism_deserialize` the old bytes, `organism_serialize14`
  the result. Nothing is lost; a VIV1 load is already homozygous.
- VIV14 -> VIV1: not supported in general (it would drop homolog 1); write
  `organism_serialize` deliberately if that loss is acceptable.
- No codec or genome bytes change, so strands copied between revisions stay
  valid.

## Test plan (implemented in `test_cells`/`test_organisms`)

- version sniffing: VIV1 -> 1, VIV14 -> 14, garbage -> 0;
- VIV14 roundtrip preserves a divergence: damage one homolog, save with
  `organism_serialize14`, load, and check both homologs came back as stored;
- the loaded pair still repairs: `organism_read` returns the original payload;
- truncated/bad-flags VIV14 containers are rejected, not misparsed;
- VIV1 roundtrip is unchanged (regression).

## VIV14N: parity genes in the chromosome

A chromosome with `chr_opts.parity = m > 0` uses the second centromere
revision and appends `m` parity genes after the data genes:

```
[TELOMERE][CEN2][data genes 0..n-1][parity genes n..n+m-1][TELOMERE]

CEN2 = CMARK (3) + 32 codons (24), raw block (12 bytes, tagged):
  id (1) | flags (1) | ngenes (2 BE) | generation (2 BE)
  | parity (1) | reserved (1) | rawlen (4 BE)
```

- Data genes are the payload split by `gene_raw` (last one shorter); parity
  genes carry `gene_raw` bytes each and are marked with bit 1 of the gene
  `type` byte.
- Shards are the data genes zero-padded to `gene_raw`; parity shards are a
  systematic Reed-Solomon code over GF(256) (`vivi/parity.h`): **any `n` of
  the `n + m` genes reconstruct the payload**.
- `chr_read` reads directly when every data gene verifies its tag, and
  falls back to RS reconstruction when some data genes are unreadable; a
  parity gene that survives provides the shard length.
- `chr_parse` detects CEN1 vs CEN2 by validating each centromere tag, so
  VIV1/VIV14 chromosomes are unchanged and remain byte-identical.
- Limits: `parity` in 1..16, `n + m <= 255` (gene ids are one byte).
- `chr_set_generation` rewrites the matching centromere revision, keeping
  parity and `rawlen`.

Container: `organism_serialize14n` writes the `"VIV14N"` magic (6 bytes)
with a parity byte per chromosome; `organism_container_version` returns
`141`. `organism_deserialize` reads VIV1, VIV14 and VIV14N; when a VIV14
container carries CEN2 chromosomes, the parity count is inferred from the
strand.

## VIV14NB4NSH33 "Vivian Banshee": on-strand primer sites

The final revision adds a **primer site on each side of the chromosome**,
so the molecule carries its own random-access address:

```
[TELOMERE][PRIMER_L][CEN1|CEN2][data genes][parity genes][PRIMER_R][TELOMERE]

PRIMER = PMARK (3, bytes 0x55) + 16 codons (12)
  raw block (8 bytes, tagged): barcode (2 BE, 1..65535) | reserved (2 = 0)
  tag: Chaskey32 over the 4 raw bytes
```

- `chr_opts.primer` selects the barcode (0 = no sites, everything below
  off). Both sites carry the same barcode; `chr_parse` validates each by
  its tag and checks they match: `rec.primer`, `rec.primer_ok`,
  `rec.primer_bytes`.
- A damaged forward site loses the layout anchor, so the chromosome does
  not parse (in the lab, a lost forward primer means the molecule cannot
  be amplified). A damaged reverse site keeps the data readable but
  `chr_amplifiable(rec)` turns false: physical access is lost while the
  information is intact.
- `chr_set_generation` rewrites the centromere in place and keeps both
  primer sites byte for byte.
- VIV1/VIV14/VIV14N chromosomes are unchanged: primer detection is
  marker+tag based, and the telomere autodetect subtracts the
  `PRIMERBYTES` block before dividing by the repeat size.

Container: `organism_serialize14nb` writes the 13-byte
`"VIV14NB4NSH33"` magic; `organism_container_version` returns `143`. The
header stores parity **and** the 16-bit barcode per chromosome; earlier
containers infer both from the strand when it carries CEN2/primer blocks.

Access model: the pool layer (`vivi_pool_amplify`) gains `p_primer`, the
probability that the on-strand site itself is hit: such a molecule is
reported dropped even when the synthetic primer pair is available
(`vivi_sim --access --p-primer`). Layout freeze: no further revisions are
planned; any change after this point is a new format family.

## Out of scope for every revision

- compression/checksums at the container level (genes already carry tags),
- error-correction beyond systematic Reed-Solomon and homology.
