# Research notes

Vivian's research layer simulates a DNA storage pipeline: encode a payload
as a diploid genome, pass the strands through a synthesis/sequencing error
channel, then let the existing repair machinery try to recover the data.
Everything is deterministic per seed, so every number below is reproducible
from the command line.

## Format evolution

The container magic evolves in named revisions; `VIV1` is the current
format and stays readable forever:

| version | meaning |
|---|---|
| `VIV1` | current container (4-byte magic) |
| `VIV14` | next revision (planned) |
| `VIV14N` | ... |
| `VIV14NB4NSH33` | final planned revision, "Vivian Banshee" |

Future magics are longer than 4 bytes, so deserializers will dispatch on a
version prefix instead of comparing a fixed 4-byte field. Until a format
change actually ships, the magic remains `VIV1` and the byte format does
not move.

## The channel (`vivi/channel.h`)

The channel is a C-only research layer: it never changes the byte format,
so the Lua port needs no counterpart.

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
- The damaged sequence is repacked to whole bytes; `bases` is the
  authoritative length (0..3 filler bases are appended as `A`).
- An identity channel (`rates = 0`, or `opts = NULL`) is a lossless read.
- The research layer uses SplitMix64 (`vivi_prng`); the format-level
  xorshift32 (`vivi_rng`) stays frozen for Lua byte parity.

## Random access (`vivi/pool.h`)

A pool stores each gene as a separate molecule, and one amplification picks
one address:

```c
vivi_pool pool = { 0 };
vivi_pool_add_chromosome(&pool, strand, slen, &err);   /* one copy per gene */

vivi_amp_opts ao = { .p_access = 0.05, .p_cross = 0.01, .seed = 7,
                     .ch = { .p_sub = 1e-3 } };
vivi_amp_result ar;
vivi_pool_amplify(&ar, &pool, target_id, &ao, &err);
/* ar.read = damaged amplicon, ar.id = molecule actually amplified
   (-1 on primer failure; a differing id is an off-target product) */
```

Primer sequences live in the primer database, exactly as in the lab; the
pool is that database plus the molecules. On-strand primer sites are planned
for the `VIV14` revision. An off-target amplicon is reported (not silently
accepted): the caller decodes with the requested id and must check.

## Running experiments (`vivi_sim`)

```bat
build.bat sim            rem build + a small sweep
vivi_sim.exe --size 1024 --gene-raw 64 --trials 1000 --p-sub 0.001 --header
vivi_sim.exe --size 1024 --gene-raw 64 --access --trials 1000 --p-sub 0.001 ^
             --p-access 0.05 --p-cross 0.01
vivi_sim.exe --in payload.bin --out run.csv --p-sub 0.0005 --p-ins 1e-5 --p-del 1e-5
```

One invocation = one CSV row; the first column selects the schema:

```
whole:  mode,p_sub,p_ins,p_del,p_drop,trials,success,wrong,failed,dropped,
        rate,avg_repaired,avg_structural,avg_dead,avg_anomaly
access: mode,p_sub,p_ins,p_del,p_drop,p_access,p_cross,trials,success,
        wrong,failed,dropped,cross,rate
```

- `success` — the payload (whole) or the target gene (access) came back exact;
- `wrong` — a read succeeded but bytes differ (this must stay 0; the 32-bit
  Chaskey tag makes a false accept ~2^-32 per gene);
- `failed` — no readable product;
- `dropped` — reads lost to `p_drop` (whole) or primer failure `p_access`
  (access); `cross` — off-target products in access mode.

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

| read strategy | p_sub | p_access | p_cross | success |
|---|---|---|---|---|
| whole chromosome + repair | 0.001 | 0 | 0 | 0.230 |
| one gene by primer | 0.001 | 0.05 | 0.01 | **0.668** |

A targeted read only needs one molecule, so it survives the same per-base
error rate far better than reading every gene — but it pays with primer
dropout (5%) and loses the diploid backup. Off-target products (`cross = 8`)
were rejected by the id check with zero silent miscalls (`wrong = 0`).

## Roadmap

1. **Channel + sim + CSV** — done (this document).
2. **Random access** — done (`vivi/pool.h`, `vivi_sim --access`); on-strand
   primer sites move to the `VIV14` revision.
3. **Outer code** — k homologs / parity genes so erasures are repairable,
   with success-vs-k curves.
4. **Fuzzing + benchmarks** — libFuzzer harnesses for every parser and a
   `make research` sweep runner.
5. **Formalization** — channel assumptions, related work (DNA Fountain,
   Goldman et al.), reproducible protocol.

All results are simulation only; no wet-lab claims are made.
