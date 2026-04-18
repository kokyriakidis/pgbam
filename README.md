# pgbam

> **Development status:** `pgbam` is still under active development and should not be used yet for production or relied on for stable results.

pgbam bridges pangenome graph alignments and the BAM ecosystem. Given a set of reads aligned to a pangenome graph (GAF format) and a haplotype index (GBWT or GBZ), it annotates each read in the BAM file with the haplotype paths it traverses and writes a compact sidecar file that maps those annotations back to named samples and haplotypes.

The typical workflow is:
1. Align reads to a pangenome graph (e.g. with `vg giraffe`) → produces a GAF file.
2. Run `pgbam annotate` with the BAM, GAF, and graph index → produces an annotated BAM and a `.pgbam` sidecar.
3. If you want a human-readable view of the thread identities stored in the sidecar, run `pgbam decode` to export them as a TSV.

---

## End-to-end example

**Scenario:** HiFi reads from sample HG002 aligned to the HPRC pangenome graph with `vg giraffe`. The goal is to annotate each read with the donor haplotypes it is most consistent with, and export a TSV mapping those annotations to named samples.

**Inputs:**

```
HG002.bam                        # coordinate-sorted BAM from your aligner
HG002.gam                        # graph alignments from vg giraffe (GAM format)
hprc-v2.1-mc-chm13-eval.gbz     # GBZ pangenome graph index
hprc-v2.1-mc-chm13-eval.gbwt    # GBWT haplotype index
hprc-v2.1-mc-chm13-eval.ri      # r-index (optional, speeds up annotation)
```

**Step 1 — Convert GAM to GAF** (skip if your aligner already produced a GAF):

```bash
vg convert -G HG002.gam hprc-v2.1-mc-chm13-eval.gbz > HG002.gaf
```

**Step 2 — Sort inputs by query name:**

```bash
samtools sort -n -@ 8 -o HG002.qname.bam HG002.bam
sort -k1,1V -k12,12nr HG002.gaf > HG002.qname.mapq.gaf
```

**Step 3 — Annotate:**

```bash
pgbam annotate \
    --bam      HG002.qname.bam \
    --gaf      HG002.qname.mapq.gaf \
    --gbwt     hprc-v2.1-mc-chm13-eval.gbwt \
    --r-index  hprc-v2.1-mc-chm13-eval.ri \
    --out-bam  HG002.annotated.bam \
    --out-sets HG002.pgbam \
    --threads  8
```

Each read in `HG002.annotated.bam` now carries three BAM tags:

```
hs:B:I,4,7    <- set IDs: two subpath matches
hb:B:I,0,391  <- start offset within node for each match
he:B:I,391,520 <- end offset (exclusive) for each match
```

`hs` indexes into `HG002.pgbam`. A read with a single `hs` entry maps entirely within one group of haplotype paths; multiple entries mean the read spans graph regions where different haplotypes diverge.

**Step 4 — Decode to TSV:**

```bash
pgbam decode \
    --sets HG002.pgbam \
    --gbwt hprc-v2.1-mc-chm13-eval.gbwt \
    --out  HG002.threads.tsv
```

`HG002.threads.tsv` maps every set ID to named haplotypes:

```
set_id  thread_id  path_id  sample   haplotype  locus         path_name
4       54312      54312    GRCh38   0          chr9          GRCh38#0#chr9[68220832]
4       54372      54372    HG00140  1          CM087114.1    HG00140#1#CM087114.1#59589124
7       61204      61204    HG00513  2          CM088003.1    HG00513#2#CM088003.1#71304981
```

Each row is one haplotype thread that passes through the graph subpath a read was aligned to. Joining on `hs` from the BAM to `set_id` in this TSV gives you the full haplotype identity for every read.

---

## Dependencies

### Required

| Dependency | Version | Install |
|---|---|---|
| CMake | ≥ 3.20 | `brew install cmake` / `apt install cmake` |
| C++20 compiler | GCC ≥ 11 or Clang ≥ 13 | system toolchain |
| HTSlib | any recent | `brew install htslib` / `apt install libhts-dev` |
| OpenSSL | any recent | `brew install openssl` / `apt install libssl-dev` |
| Zstandard | any recent | `brew install zstd` / `apt install libzstd-dev` |

### Optional

| Dependency | Purpose | Install |
|---|---|---|
| patchelf | **Linux only** — self-contained install (see [Installing](#installing)) | auto-built if not found on `PATH` |

### Bundled (via git submodules)

- [sdsl-lite](https://github.com/vgteam/sdsl-lite) — compressed data structures
- [gbwt](https://github.com/jltsiren/gbwt) — graph BWT and r-index
- [gbwtgraph](https://github.com/jltsiren/gbwtgraph) — graph query layer
- [libhandlegraph](https://github.com/vgteam/libhandlegraph) — pangenome graph interface

---

## Building

Clone with submodules:

```bash
git clone --recurse-submodules https://github.com/kokyriakidis/pgbam.git
cd pgbam
```

Configure and build:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The `pgbam` binary is at `build/pgbam`.

---

## Installing

```bash
cmake --install build --prefix /path/to/install
```

This installs:

```
/path/to/install/
  bin/pgbam
  lib/libpgbam.dylib   (macOS) / libpgbam.so  (Linux)
  lib/libhts.*         ← bundled at install time
  lib/libzstd.*        ← bundled at install time
  include/pgbam/
  lib/cmake/pgbam/
```

The install step automatically copies all non-system runtime libraries into `lib/` and rewrites their load paths so the installation is self-contained — no Homebrew, conda, or package manager prefix needs to be present on the target machine.

**Linux note:** `patchelf` is required for the bundling step. If it is not on `PATH`, it is automatically built from source during `cmake --build` and used transparently — no manual install needed.

---

## Usage

### `pgbam annotate`

Annotate a BAM file with graph thread information.

> **GAM input:** pgbam accepts GAF only. If your aligner produced a GAM file (vg's binary protobuf format), convert it first:
> ```bash
> vg convert -G alignments.gam graph.gbz > alignments.gaf
> ```
> The GBZ is required — `vg convert` needs the graph to compute target lengths and resolve node sequences for the GAF records.

```
pgbam annotate
    --bam      <in.bam>
    --gaf      <in.gaf>
    --out-bam  <out.bam>
    --out-sets <out.pgbam>
  ( --gbz    <graph.gbz>
  | --gbwt   <graph.gbwt> )
  [ --r-index <graph.ri> ]
  [ --threads <N> ]
  [ --primary-only ]
```

| Flag | Required | Description |
|---|---|---|
| `--bam` | yes | Input BAM file |
| `--gaf` | yes | Input GAF file (read-to-graph alignments) |
| `--out-bam` | yes | Output annotated BAM file |
| `--out-sets` | yes | Output sidecar metadata file (`.pgbam`) |
| `--gbz` | one of | GBZ graph index |
| `--gbwt` | one of | GBWT graph index |
| `--r-index` | no | R-index (`.ri`) for the GBWT — enables faster locate queries |
| `--threads` | no | Worker threads (default: 1) |
| `--primary-only` | no | Use only the first GAF alignment in each `qname` block; annotate every BAM record with that `qname` from that primary alignment's set |

The r-index accelerates the locate step — recovering which haplotype paths pass through each aligned subpath. For large cohorts or deeply sequenced samples, supplying it meaningfully reduces runtime.

**Memory model:**

- The BAM and GAF inputs are streamed during annotation.
- The GBZ or GBWT index is loaded into memory before annotation starts.
- If an r-index is provided, it is also loaded into memory.
- In practice, RAM usage is therefore driven primarily by the graph index files, not by the BAM/GAF join itself.

**Input requirements:**

- `pgbam annotate` requires both BAM and GAF to be non-decreasing by `qname`.
- Sort BAM by query name with `samtools sort -n -o reads.qname.bam reads.bam`.
- Sort GAF by `qname` and then decreasing MAPQ with `sort -k1,1V -k12,12nr reads.gaf > reads.qname.mapq.gaf`.
- This ordering is recommended in general and is required for `--primary-only`, because `pgbam` uses the first GAF alignment in each `qname` block as the primary one.
- Sorting by `qname` groups identical names together, but it does not create a reliable per-record order within a duplicate-name block.
- By default, `pgbam` aggregates all GAF alignments in a `qname` block, combines their graph matches, and applies the same aggregated annotation to every BAM record in the matching `qname` block.
- With `--primary-only`, `pgbam` uses only the first GAF alignment in each `qname` block and applies that annotation to every BAM record in the matching `qname` block.
- If your pipeline needs one-to-one pairing between multiple BAM records and multiple GAF records with the same read name, `qname` alone is not a sufficient join key.

**Example:**

```bash
pgbam annotate \
    --bam      reads.qname.bam \
    --gaf      reads.qname.mapq.gaf \
    --gbz      hprc.gbz \
    --out-bam  reads.annotated.bam \
    --out-sets reads.pgbam \
    --r-index  hprc.ri \
    --threads  8
```

Primary-only example:

```bash
pgbam annotate \
    --bam          reads.qname.bam \
    --gaf          reads.qname.mapq.gaf \
    --gbz          hprc.gbz \
    --out-bam      reads.annotated.bam \
    --out-sets     reads.pgbam \
    --primary-only
```

**BAM tags written to each read:**

`hs`, `hb`, and `he` are parallel arrays — entry `i` across all three describes one subpath match:

| Tag | Type | Content |
|---|---|---|
| `hs` | `B:I` | Set ID for match `i` — index into the sidecar file to retrieve the matching thread IDs |
| `hb` | `B:I` | Begin node offset for match `i` — where the alignment starts within the node sequence |
| `he` | `B:I` | End-exclusive node offset for match `i` — where the alignment ends |

A read that maps to multiple graph threads (or to the same thread at different positions) gets one entry per match, so all three arrays have the same length. `hs` alone tells you *which* haplotype paths a read belongs to; `hb`/`he` tell you *where* it sits within them. All three are needed to reconstruct the exact alignment coordinate inside the graph.

Tags are omitted on reads with no valid graph mapping. The BAM header gains a `@PG` record of the exact command and a `@CO` comment explaining the tags.

---

### `pgbam decode`

Decode a sidecar file into a TSV of thread identities.

```
pgbam decode
    --sets <in.pgbam>
    --out  <out.tsv>
  ( --gbz   <graph.gbz>
  | --gbwt  <graph.gbwt> )
```

| Flag | Required | Description |
|---|---|---|
| `--sets` | yes | Sidecar metadata file produced by `annotate` |
| `--out` | yes | Output TSV file |
| `--gbz` | one of | GBZ graph index (must match the one used at annotation) |
| `--gbwt` | one of | GBWT graph index |

**Example:**

```bash
pgbam decode \
    --sets reads.pgbam \
    --gbz  hprc.gbz \
    --out  threads.tsv
```

**Output columns:**

```
set_id  thread_id  path_id  sample  haplotype  locus  path_name
```

---

## Using pgbam as a library

If you installed with `cmake --install`, downstream CMake projects can find pgbam via:

```cmake
find_package(pgbam REQUIRED)
target_link_libraries(my_target PRIVATE pgbam::pgbam)
```

Pass `-Dpgbam_DIR=/path/to/install/lib/cmake/pgbam` if the install prefix is not in CMake's default search path.
