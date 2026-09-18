# Research notes

Vivian's research layer simulates a DNA storage pipeline: encode a payload
as a diploid genome, pass the strands through a synthesis/sequencing error
channel, then let the existing repair machinery try to recover the data.
Everything is deterministic per seed, so every number below is reproducible
from the command line.

## Format evolution

The container magic evolves in named revisions; all of them stay readable
forever:

| version | meaning |
|---|---|
| `VIV1` | original container (4-byte magic) |
| `VIV14` | versioned container, both homologs + `max_gen` |
| `VIV14N` | VIV14 + parity genes in the chromosome (CEN2) |
| `VIV14NB4NSH33` | final revision, "Vivian Banshee": on-strand primer sites; layout frozen |

Deserializers dispatch on a version prefix (`organism_container_version`:
1, 14, 141, 143, 0), so longer magics never alias the 4-byte `VIV1` field.

## The channel (`vivi/channel.h`)

The channel is a research layer: it never changes the byte format. The Lua
reference ports it too (`framework/luau/channel.lua`, `pool.lua`) with a
byte-identical SplitMix64, and `xcheck` diffs the channel, consensus and pool
scenarios on both sides.

One call turns a pristine packed strand into one **read**:

```c
vivi_channel_opts ch = { .p_sub = 1e-3, .p_ins = 0, .p_del = 0, .p_drop = 0, .seed = 42 };
vivi_read rd;
vivi_channel_read(&rd, strand, slen, &ch, &err);
/* rd.strand: damaged packed strand, rd.bases: exact base count, rd.dropped: read lost */
```

- Errors are applied per base, in this pinned order: delete?, substitute?
  (then the replacement), insert? (then the digit). This makes a read a
  pure function of `(strand, rates, seed)`.
- `p_drop` models a lost read (PCR dropout): `dropped = 1`, empty strand.
- Realistic extras, all defaulting to 0 and only drawing when enabled, so
  the pinned single-read stream is byte-identical when they are off:
  `p_sub_gc` (extra substitution on G/C), `p_sub_hp` (extra substitution
  inside a homopolymer run), `p_trunc` (the read is cut at a random base),
  and `p_burst`/`burst_len`/`p_burst_del` (a correlated run of `burst_len`
  bases; each is deleted with probability `p_burst_del`, else substituted).
  Combined substitution probability is
  `1 - (1-p_sub)(1-p_sub_gc)(1-p_sub_hp)`; a burst forces `1.0`.
- `vivi_channel_read_soft` also returns a per-base quality (`vivi_read.qual`):
  `VIVI_QUAL_HI` (60) for a base read unchanged, `VIVI_QUAL_LO` (4) for a
  substituted or inserted one. Free any read with `vivi_read_free()`.
- The damaged sequence is repacked to whole bytes; `bases` is the
  authoritative length (0..3 filler bases are appended as `A`).
- An identity channel (`rates = 0`, or `opts = NULL`) is a lossless read.
- The research layer uses SplitMix64 (`vivi_prng`); the format-level
  xorshift32 (`vivi_rng`) stays frozen for Lua byte parity.

## Inner code (`vivi/rs.h`)

The 32-bit gene tag turns any corruption into a whole-gene erasure. The
inner code shortens the erasure radius: `genome_gene_encode_inner` writes
the payload as a systematic Reed-Solomon codeword (data + `inner_m` parity
bytes, GF(256), `n <= 255`) with gene usertype bit 1 set, so the decoder
repairs up to `floor(inner_m/2)` corrupted bytes *before* the tag check.
`chr_opts.inner` selects it per chromosome (`--inner M` in `vivi_sim` and
the CLI); the parity count is self-describing, so a deserialized organism
infers it from the genes and `mutate`/`cross` preserve it. Dense mode
only, and `gene_raw + inner_m <= 255`.

A base substitution only sometimes maps to one bad byte: measured over
20 000 single-base corruption trials on a 64-byte payload, 25% made the
DNA decoder fail outright, 55% produced exactly one byte error, 2% two,
and the rest desynchronized the constraint coder into many. RS repairs
the single- and double-byte cases (the dominant ones) but cannot undo a
desync, so the inner code is a partial, honest win.

## Consensus (coverage)

`vivi_consensus_read` reads the same molecule `coverage` times and
majority-votes every base, the way a real pipeline sequences many copies
of one oligo before calling a base:

```c
vivi_consensus_opts cn = { { .p_sub = 1e-3, .seed = 42 }, .coverage = 9 };
vivi_read rd;
vivi_consensus_read(&rd, strand, slen, &cn, &err);
```

- Reads are seeded `ch.seed + r` for `r = 0..coverage-1`; the vote is per
  2-bit base and a tie picks the lowest digit, so the consensus is a pure
  function of `(strand, ch, coverage)`.
- Substitution-only: with `coverage > 1` every surviving read must keep
  the original base count, so `p_ins`/`p_del` must stay 0.
- `coverage <= 1` is exactly `vivi_channel_read`, so every pinned
  single-read stream is unchanged.
- `--coverage K` applies it in all three modes: `whole` votes each
  homolog, `library` votes each replica, `access` votes each amplicon.

## Random access (`vivi/pool.h`)

A pool stores each gene as a separate molecule, and one amplification picks
one address:

```c
vivi_pool pool = { 0 };
vivi_pool_add_chromosome(&pool, strand, slen, &err);   /* one copy per gene */

vivi_amp_opts ao = { .p_access = 0.05, .p_cross = 0.01, .seed = 7,
                     .ch = { .p_sub = 1e-3 },
                     .p_primer = 0.05 };   /* Banshee: on-strand site loss */
vivi_amp_result ar;
vivi_pool_amplify(&ar, &pool, target_id, &ao, &err);
/* ar.read = damaged amplicon, ar.id = molecule actually amplified
   (-1 on primer failure; a differing id is an off-target product) */
```

Primer sequences live in the primer database, exactly as in the lab; the
pool is that database plus the molecules. From `VIV14NB4NSH33` the primer
sites live on the chromosome itself (`chr_opts.primer`, validated by
`chr_amplifiable`), and `p_primer` models their dropout: the draw happens
only when the rate is non-zero, so the pinned stream of earlier revisions
is unchanged. An off-target amplicon is reported (not silently accepted):
the caller decodes with the requested id and must check.

## Outer code (`vivi/parity.h`)

Data genes are split into `n` equal-length shards; `m` parity shards are
computed with a systematic Reed-Solomon code over GF(256), so **any** `n` of
the `n + m` molecules reconstruct the data. Each shard is stored as a normal
gene molecule (its own tag and id), and `--replicas K` additionally
synthesizes `K` copies of every molecule. `vivi_sim --library` reads all
molecules through the channel, marks a gene present when any replica
survives intact, and calls `vivi_parity_decode`.

```bat
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 ^
             --replicas 2 --parity 4 --header
```

## Running experiments (`vivi_sim`)

```bat
build.bat sim            rem build + a small sweep
vivi_sim.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub 0.001 --header
vivi_sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 ^
             --p-access 0.05 --p-cross 0.01
vivi_sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 ^
             --p-access 0.05 --p-cross 0.01 --p-primer 0.05
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 ^
             --coverage 5
vivi_sim.exe --size 1024 --gene-raw 64 --library --trials 1000 --p-sub 0.001 ^
             --coverage-lambda 5
vivi_sim.exe --in payload.bin --out run.csv --p-sub 0.0005 --p-ins 1e-5 --p-del 1e-5
```

One invocation = one CSV row; the first column selects the schema:

```
whole:   mode,p_sub,p_ins,p_del,p_drop,p_sub_gc,p_sub_hp,p_trunc,p_burst,p_burst_del,coverage,cov_lambda,soft,inner,trials,
         success,wrong,failed,dropped,rate,avg_repaired,avg_structural,avg_dead,avg_anomaly
access:  mode,p_sub,p_ins,p_del,p_drop,p_sub_gc,p_sub_hp,p_trunc,p_burst,p_burst_del,p_access,p_cross,p_primer,
         coverage,cov_lambda,soft,inner,trials,success,wrong,failed,dropped,cross,rate
library: mode,p_sub,p_ins,p_del,p_drop,p_sub_gc,p_sub_hp,p_trunc,p_burst,p_burst_del,coverage,cov_lambda,soft,inner,replicas,
         parity,genes,trials,success,wrong,failed,dropped,rate,avg_present
```

- `success` — the payload (whole) or the target gene (access) came back exact;
- `wrong` — a read succeeded but bytes differ (this must stay 0; the 32-bit
  Chaskey tag makes a false accept ~2^-32 per gene);
- `failed` — no readable product;
- `dropped` — reads lost to `p_drop` (whole), primer failure `p_access`, or
  on-strand site loss `p_primer` (access); `cross` — off-target products.

Sweep by looping the tool, e.g. `for /l %p in (1,1,10) do vivi_sim.exe ...`
or a shell loop on POSIX; concatenate the CSV rows (use `--header` once).

## What we measured (single chromosome, 256 B payload, `h=3`, `gene_raw=1024`)

| p_sub | p_drop | success rate |
|---|---|---|
| 0.0005 | 0 | 0.823 |
| 0.001 | 0 | 0.523 |
| 0.002 | 0 | 0.183 |
| 0.004 | 0.01 | 0.023 |

Reading the curve: with per-gene 32-bit tags and 2× replication only, a
gene is lost unless at least one homolog's copy is untouched. That makes
tolerance scale with gene size (a 1 KiB gene is far easier to hit than a
64 B one) and drops sharply with the error rate. This is exactly the gap
an inner error-correcting code would close, and it is the first result the
research layer exists to show.

### Random access vs whole-chromosome read

1 KiB payload, `gene_raw=64` (16 genes), 1000 trials each:

| read strategy | p_sub | p_access | p_cross | p_primer | success |
|---|---|---|---|---|---|
| whole chromosome + repair | 0.001 | 0 | 0 | 0 | 0.230 |
| one gene by primer | 0.001 | 0.05 | 0.01 | 0 | **0.668** |
| one gene, on-strand site lost 5% | 0.001 | 0.05 | 0.01 | 0.05 | 0.634 |

A targeted read only needs one molecule, so it survives the same per-base
error rate far better than reading every gene — but it pays with primer
dropout and loses the diploid backup. Banshee's on-strand sites add their
own dropout (`p_primer`); off-target products (`cross`) are rejected by the
id check with zero silent miscalls (`wrong = 0`).

### Outer code: replication (K) and parity (m)

1 KiB payload, `gene_raw=64` (n = 16), p_sub = 0.001, 1000 trials:

| redundancy | success | avg present |
|---|---|---|
| none (K=1, m=0) | 0.006 | 11.4 / 16 |
| K = 3 replicas | 0.704 | 15.7 / 16 |
| m = 4 parity | 0.274 | 14.3 / 20 |
| K = 2 + m = 4 | **0.985** | 18.4 / 20 |

Replication is cheap and effective but scales linearly with storage; parity
already at `m = 4` recovers many patterns replication cannot (any 4 erasures
among 20), and the combination is nearly lossless at 2.5× storage. With no
redundancy the baseline is 0.006 — every one of 16 genes must survive, and
each has only ~67% chance at this error rate.

### Consensus: reads per molecule turn erasures into corrections

1000 trials each; coverage = reads majority-voted before the tag check:

| mode | parameters | coverage 1 | coverage 5 | coverage 15 |
|---|---|---|---|---|
| whole | 256 B, p_sub = 0.002 | 0.183 | **1.000** | 1.000 |
| library | 1 KiB, n=16, p_sub = 0.001 | 0.006 | **1.000** | — |
| access | 1 KiB, p_sub = 0.001, p_access = 0.05 | 0.668 | **0.936** | — |

Consensus is the cheapest lever for substitutions: at 0.2% per base a
9-fold read turns a 5-fold loss into a total recovery without spending any
storage redundancy, because a base is only wrong if more than half its
copies are wrong. It composes with the outer code (the `library` row uses
K=1, m=0), but it multiplies read cost, so the trade is reads vs. molecules.
Random access stays bounded by primer dropout (`p_access`), which no amount
of coverage fixes: the 0.936 ceiling is its 5% loss, not a sequencing error.

### Soft-decision (`--soft`)

A hard vote counts every read equally. With `soft != 0` each read uses
`vivi_channel_read_soft`: an unchanged base scores `VIVI_QUAL_HI` (60), a
substituted or inserted one `VIVI_QUAL_LO` (4), and the consensus vote is
weighted by `quality + 1`, so one high-quality base outweighs several
low-quality ones. Hard vs soft at coverage 3, p_sub = 0.4, one strand:

| vote | mismatched bases (of 1024) |
|---|---|
| hard majority | 20 |
| soft-weighted | **5** |

Soft-decision is what real basecallers give you; modelling it closes the
"a failed tag is a hard erasure" caveat for the consensus path.

### Coverage distribution (`--coverage-lambda`)

Fixed `--coverage K` assumes every molecule gets exactly K reads. Real pools
do not: `--coverage-lambda L` draws each molecule's read count from
Poisson(L) (a molecule drawn 0 times is simply absent). 1000 trials:

| mode | fixed K | success | Poisson L | success |
|---|---|---|---|---|
| whole | 5 | 1.000 | 5 | 0.977 |
| library | 5 | 1.000 | 5 | 0.518 |
| access | 5 | 0.936 | 5 | 0.919 |

The idealized fixed count flatters the result: with a realistic spread some
whole trials lose a homolog entirely (15 dropped reads) and, far more
importantly, the library needs all 16 molecules covered, so uneven coverage
drops it from 1.000 to 0.518 even at the same mean. This is the first place
where the "coverage distribution" gap in the assumptions bites.

### Inner code: bounded erasures, and where it stops

| mode | parameters | inner 0 | inner 16 |
|---|---|---|---|
| whole | 256 B, gene_raw 64, p_sub 0.001 | 0.689 | **0.903** |
| library | 1 KiB, n=16, gene_raw 64, p_sub 0.001 | 0.006 | **0.066** |
| access | 1 KiB, gene_raw 64, p_sub 0.001, p_access 0.05 | 0.668 | **0.801** |

The inner code trades capacity (shorter genes) for per-molecule repair, so
it lifts whole and access strongly but library only modestly: with 16 genes
all needing to survive, the unavoidable ~25% decode-desync share dominates.
It composes with parity and coverage.

## Benchmark (`tools/vivi_bench`)

`vivi_bench` compares recovery strategies at their real storage cost for
the same payload and channel (a separate module, like the CLI):

```bat
vivi_bench.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub 0.001 --header
```

1 KiB, 16 genes, p_sub = 0.001, 1000 trials (`storage_x` = stored packed
bytes / payload bytes):

| strategy | storage_x | success |
|---|---|---|
| plain (no code) | 1.44 | 0.006 |
| replica K=2 | 2.87 | 0.247 |
| replica K=3 | 4.31 | 0.704 |
| parity m=8 | 2.15 | 0.776 |
| inner m=4 | 1.50 | 0.076 |
| coverage K=5 | 1.44 | **1.000** |

The headline: coverage buys a full recovery at zero storage cost (only
reads), while replication buys it at 4.3×; parity and the inner code sit in
between. The tool prints one row per strategy so the trade is explicit
rather than asserted.

## Fuzzing

`test/fuzz/` holds one libFuzzer harness per parser: `dna.c` (strand, ASCII
view, channel and consensus), `gene.c`, `rs.c` (the inner codec), `inner.c`
(inner-coded genes/chromosomes), `chr.c`, `viv1.c` (container), `pool.c` and
`parity.c`. Each harness pushes a cached valid sample through the parser
first, then the fuzz input, so the valid path runs on every iteration.

```bat
make fuzz          rem build every harness (clang + libFuzzer, POSIX)
make fuzz-smoke    rem 10 s per harness; this is what CI runs
```

libFuzzer requires `-fsanitize=fuzzer` and is not available under Windows
clang; the harnesses are POSIX-only (the library itself still builds
everywhere).

## Channel assumptions

Modelled: independent per-base substitution, insertion and deletion with
fixed rates, whole-read dropout, whole-read truncation (`p_trunc`), context
bias (extra substitutions on G/C and inside homopolymers), correlated bursts
of substitutions **or** deletions (`p_burst`/`burst_len`/`p_burst_del`), and
per-base quality with soft-decision consensus (`vivi_channel_read_soft`,
`soft`), all deterministic per seed. Per-molecule coverage is
`coverage` independent reads with majority voting (`vivi_consensus_read`) or
a Poisson count per molecule (`--coverage-lambda`), and the inner code
(`vivi/rs.h`) corrects byte errors inside a molecule.

Not modelled (yet): PCR amplification bias, per-oligo synthesis dropout
separate from `p_drop`, adapter contamination, read quality calibrated to a
real basecaller (our quality is the ideal unchanged/changed split), and
per-platform error rates. Treat absolute
numbers as comparative, not as predictions for a specific platform.

## Related work

- G. M. Church, Y. Gao, S. Kosuri, *Science* 337 (2012) — DNA as a storage medium.
- N. Goldman et al., *Nature* 494 (2013) — the first practical-scale DNA archive.
- R. Grass et al., *Angew. Chem. Int. Ed.* 54 (2015) — chemical preservation.
- S. M. H. T. Yazdi et al., *Sci. Rep.* 5 (2015) — rewritable random-access DNA.
- Y. Erlich, D. Zielinski, *Science* 355 (2017) — DNA Fountain; `vivi/parity.h`
  is the systematic Reed-Solomon baseline that such schemes build on.
- L. Organick et al., *Nat. Biotechnol.* 36 (2018) — random access with PCR
  primers; the access model in `vivi/pool.h`.

Vivian is a simulation and teaching testbed with two byte-exact
implementations; it does not claim to compete with these systems.

## Threats to validity

- The channel is synthetic; rates are chosen for curve shape, not measured
  on hardware.
- A failed tag is treated as an erasure; real pipelines carry soft
  information and non-zero false-accept rates.
- Parity decoding assumes correct shard classification and equal-length
  shards (zero padded).
- Absolute success depends on payload size, `gene_raw` and `h`; only
  like-for-like rows are comparable.

## Reproducible protocol

- All randomness is integer SplitMix64 seeded per read; the same flags
  produce the same CSV bytes on any platform.
- `make research` emits `research_whole.csv`, `research_access.csv` and
  `research_library.csv`; every row carries the full parameter set.
- Report numbers together with the commit hash and the CSV rows.

## Format evolution

`VIV14` shipped: a versioned container that stores **both homologs** and
`max_gen`, with VIV1 still readable. `VIV14N` shipped: parity genes as
first-class chromosome members (CEN2 centromere, RS reconstruction in
`chr_read`). `VIV14NB4NSH33` ("Banshee") shipped: on-strand primer sites
and the layout freeze. Full design and wire format: [viv14.md](viv14.md).

The Lua reference is at byte parity with the C port for every container
revision, enforced by the `xcheck` scenario pair in CI. The inner code
(gene usertype bit 1) is ported too (`framework/luau/rs.lua`), and the
`xcheck` scenarios cover it: the RS codeword, a two-byte correction, a
three-byte rejection, an inner-coded chromosome read, and inner inference
through a Banshee container.

## Roadmap

1. **Channel + sim + CSV** — done.
2. **Random access** — done (`vivi/pool.h`, Banshee primer sites).
3. **Outer code** — done (`vivi/parity.h`).
4. **Fuzzing + benchmarks** — done (`test/fuzz/`, `make fuzz-smoke`).
5. **Formalization** — done.
6. **Consensus / coverage** — done (`vivi_consensus_read`, `--coverage K`).
7. **Inner code** — done in C and Lua (`vivi/rs.h`, `framework/luau/rs.lua`,
   gene usertype bit 1, `chr_opts.inner`, `--inner M`); `xcheck` covers it.
8. **Realistic channel** — done (context rates, truncation, bursts:
   `--p-sub-gc`, `--p-sub-hp`, `--p-trunc`, `--p-burst`).
9. **Equal-redundancy benchmark** — done (`tools/vivi_bench`).
10. **Coverage distribution** — done (`--coverage-lambda`, Poisson reads).
11. **Model layer** — done (`vivi/model.h`: fitness, population, selection).
12. **Lua research parity** — done (`framework/luau/channel.lua`, `pool.lua`;
    byte-identical SplitMix64, `xcheck` covers channel/consensus/pool).
13. **Soft-decision** — done (`vivi_channel_read_soft`, `--soft`; quality-
    weighted consensus).
14. **Correlated indels** — done (`p_burst_del`: deletions inside a burst).

Still open (honest limits, not missing code): a basecaller-calibrated quality
model, PCR amplification bias, per-oligo synthesis dropout, adapter
contamination, and the inner-code desync (~25% of single substitutions make
`dna_decode` fail before the RS decoder can help).

All results are simulation only; no wet-lab claims are made.
