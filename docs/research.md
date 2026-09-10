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

## Running experiments (`vivi_sim`)

```bat
build.bat sim            rem build + a 4-point sweep
vivi_sim.exe --size 256 --trials 1000 --p-sub 0.001 --header
vivi_sim.exe --in payload.bin --out run.csv --p-sub 0.0005 --p-ins 1e-5 --p-del 1e-5
```

One invocation = one CSV row:

```
p_sub,p_ins,p_del,p_drop,trials,success,wrong,failed,dropped,
rate,avg_repaired,avg_structural,avg_dead,avg_anomaly
```

- `success` — `organism_read()` returned the exact payload;
- `wrong` — a read succeeded but bytes differ (this must stay 0; the 32-bit
  Chaskey tag makes a false accept ~2^-32 per gene);
- `failed` — no readable homolog;
- `dropped` — reads lost to `p_drop` across all trials.

Sweep by looping the tool, e.g. `for /l %p in (1,1,10) do vivi_sim.exe ...`
or a shell loop on POSIX; concatenate the CSV rows (use `--header` once).

## What we measured (single chromosome, 256 B payload, `h=3`, `gene_raw=1024`)

| p_sub | p_drop | success rate |
|---|---|---|
| 0.0005 | 0 | 0.817 |
| 0.001 | 0 | 0.518 |
| 0.002 | 0 | 0.168 |
| 0.004 | 0.01 | 0.027 |

Reading the curve: with per-gene 32-bit tags and 2× replication only, a
gene is lost unless at least one homolog's copy is untouched. That makes
tolerance scale with gene size (a 1 KiB gene is far easier to hit than a
64 B one) and drops sharply with the error rate. This is exactly the gap
an inner error-correcting code would close, and it is the first result the
research layer exists to show.

## Roadmap

1. **Channel + sim + CSV** — done (this document).
2. **Random access** — read one gene by id ("PCR primer"), model primer
   dropout and localized amplification.
3. **Outer code** — k homologs / parity genes so erasures are repairable,
   with success-vs-k curves.
4. **Fuzzing + benchmarks** — libFuzzer harnesses for every parser and a
   `make research` sweep runner.
5. **Formalization** — channel assumptions, related work (DNA Fountain,
   Goldman et al.), reproducible protocol.

All results are simulation only; no wet-lab claims are made.
