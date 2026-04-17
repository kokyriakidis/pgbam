# pgbam

Annotate BAM files with pangenome graph path information.

`pgbam annotate` takes a BAM file, a GAF file mapping reads to a pangenome graph, and a GBWT/GBZ graph index. It writes an annotated BAM where each read carries auxiliary tags describing which graph threads (haplotype paths) it traverses, plus a sidecar metadata file for decoding.

`pgbam decode` converts the sidecar file back into a human-readable TSV of thread identities.

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

### Optional but recommended

| Dependency | Purpose | Install |
|---|---|---|
| OpenMP | multi-threaded annotation | `brew install libomp` / `apt install libomp-dev` |
| patchelf | **Linux only** — self-contained install (see [Installing](#installing)) | `apt install patchelf` / `conda install patchelf` |

### Bundled (via git submodules)

- [sdsl-lite](https://github.com/vgteam/sdsl-lite) — compressed data structures
- [gbwt](https://github.com/jltsiren/gbwt) — graph BWT index
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
  lib/libomp.*         ← bundled at install time (if used)
  ...
  include/pgbam/
  lib/cmake/pgbam/
```

The install step automatically copies all non-system runtime libraries into `lib/` and rewrites their load paths so the installation is self-contained — no Homebrew, conda, or package manager prefix needs to be present on the target machine.

**Linux note:** `patchelf` must be on `PATH` at install time for the bundling step to run. If it is not found, a warning is printed and the install proceeds without bundling (the binary will depend on system-installed libraries at runtime).

---

## Usage

### `pgbam annotate`

Annotate a BAM file with graph thread information.

```
pgbam annotate --bam <in.bam> --gaf <in.gaf> \
               (--gbz <graph.gbz> | --gbwt <graph.gbwt>) \
               --out-bam <out.bam> --out-sets <out.pgbam> \
               [--r-index <graph.ri>] [--threads <N>]
```

| Flag | Required | Description |
|---|---|---|
| `--bam` | yes | Input BAM file |
| `--gaf` | yes | Input GAF file (read-to-graph alignments) |
| `--gbz` | one of | GBZ graph index |
| `--gbwt` | one of | GBWT graph index |
| `--out-bam` | yes | Output annotated BAM file |
| `--out-sets` | yes | Output sidecar metadata file (`.pgbam`) |
| `--r-index` | no | R-index for the graph (speeds up path lookup) |
| `--threads` | no | Worker threads (default: 1) |

**Example:**

```bash
pgbam annotate \
  --bam reads.bam \
  --gaf reads.gaf \
  --gbz hprc.gbz \
  --out-bam reads.annotated.bam \
  --out-sets reads.pgbam \
  --threads 8
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
pgbam decode --sets <in.pgbam> \
             (--gbz <graph.gbz> | --gbwt <graph.gbwt>) \
             --out <out.tsv>
```

| Flag | Required | Description |
|---|---|---|
| `--sets` | yes | Sidecar metadata file produced by `annotate` |
| `--gbz` | one of | GBZ graph index (must match the one used at annotation) |
| `--gbwt` | one of | GBWT graph index |
| `--out` | yes | Output TSV file |

**Example:**

```bash
pgbam decode \
  --sets reads.pgbam \
  --gbz hprc.gbz \
  --out threads.tsv
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
