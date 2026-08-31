---
icon: lucide/rocket
---

# Welcome to BroCOLI's documentation!

## Introduction

**BroCOLI** (**Bro**n-Kerbosch **c**alibrator **o**f **L**ong-read **I**soform) leverages Bron-Kerbosch maximal-clique enumeration for transcript identification and quantification from long-read RNA-Seq data, supporting bulk, single-cell and spatial applications, while maintaining low memory usage and fast performance for large-scale datasets.

The toolkit ships as **two executables** that cover the two halves of a long-read experiment:

| Executable | Role | When you need it |
| :--- | :--- | :--- |
| **`prebrocoli`** | Cell barcode and UMI extraction from raw long-read FASTQ | Single-cell and spatial libraries only |
| **`brocoli`** | Isoform detection and quantification from aligned reads | Every workflow — bulk, single-cell and spatial |

``` mermaid
%%{init: {'themeVariables': { 'fontSize': '18px' }}}%%
graph LR
  classDef input fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,color:black;
  classDef process fill:#fff9c4,stroke:#fbc02d,stroke-width:2px,color:black;
  classDef result fill:#e8f5e9,stroke:#2e7d32,stroke-width:2px,color:black;

  A[raw fastq.gz]:::input
  C[demultiplexed fastq]:::process
  B{sorted bam}:::process
  E[gene &amp; transcript counts]:::result

  A --> |"Bulk: minimap2"| B
  A --> |"SC / Spatial: prebrocoli"| C
  C --> |"minimap2"| B
  B --> |"brocoli"| E

  linkStyle 0 stroke:#ff9800,stroke-width:4px;
  linkStyle 1,2 stroke:#9c27b0,stroke-width:4px;
  linkStyle 3 stroke:#2e7d32,stroke-width:4px;
```

---

## Requirements

BroCOLI requires a **C++17** compatible compiler (g++ 7 or later) and the following libraries:

| Dependency | Used by | Purpose |
| :--- | :--- | :--- |
| **htslib** | `brocoli` | SAM/BAM parsing (located through `pkg-config`) |
| **WFA2-lib** | `prebrocoli` | gap-affine alignment, C++ bindings (`bindings/cpp/WFAligner.hpp`) |
| **edlib** | `prebrocoli` | fast edit-distance search for the adapter/flank |
| **zlib** | both | gzipped FASTQ input |

All of them are declared in `environment.yml`, so the recommended route is to let conda provide the whole toolchain.

!!! note "edlib is vendored, not downloaded"
    `edlib.h` and `edlib.cpp` must sit in `third_party/edlib/`. If they are missing, the build still produces `brocoli`, prints what to copy, and exits cleanly — you simply will not get `prebrocoli`.

---

## :rocket: Installation

BroCOLI is available on [`BroCOLI GitHub`][BroCOLInew].

  [BroCOLInew]:https://github.com/gyjames/BroCOLI

``` shell title="download"
git clone https://github.com/gyjames/BroCOLI.git
cd BroCOLI
```

### Option A — conda (recommended)

One command creates (or updates) the environment, runs the unit tests, and builds both executables:

``` shell title="one-shot setup"
bash setup_conda.sh
```

``` shell title="verify"
conda activate brocoli
./brocoli -h
./prebrocoli -h
```

The script is idempotent: re-running it after a source change simply rebuilds.

### Option B — make

Inside an environment that already provides the dependencies:

``` shell title="manual build"
make check-env      # htslib / WFAligner.hpp / edlib must all report "found"
make                # builds brocoli and prebrocoli
make test           # unit tests only
```

`make check-env` reads `$CONDA_PREFIX`, so activate the environment first. Both routes leave the executables in the repository root.

??? info "Repository layout"
    ```
    BroCOLI/
    ├── src/                  BroCOLI — detection and quantification
    ├── pre/                  preBroCOLI — barcode and UMI extraction
    ├── tests/                unit tests
    ├── third_party/edlib/    place edlib.h and edlib.cpp here
    ├── environment.yml       conda dependencies
    ├── setup_conda.sh        one-shot environment setup and build
    ├── Makefile
    └── CMakeLists.txt
    ```
    `src/` and `pre/` are two independent source trees compiled separately. They each have their own `main.cpp` **and their own `common.hpp`**, so their files must never be mixed into one directory.

??? failure "Build troubleshooting"
    | Symptom | Cause | Fix |
    | :--- | :--- | :--- |
    | `pkg-config cannot find htslib` | environment not activated, or htslib missing | `conda activate brocoli`, then `conda list htslib` |
    | `could not compile against the WFA2-lib C++ bindings` | `wfa2-lib` not installed in this environment | `conda env update -n brocoli -f environment.yml --prune` |
    | Only `brocoli` appears, no `prebrocoli` | `third_party/edlib/` is empty | copy `edlib.h` **and** `edlib.cpp` there, re-run |
    | `expected 11 headers under src/, found N` | incomplete checkout, or files moved between `src/` and `pre/` | re-clone, do not rearrange the tree |

---

## :sparkles: Supported sequencing data

| Platform / chemistry | `prebrocoli -q` | Status |
| :--- | :--- | :---: |
| ONT cDNA / ONT dRNA (bulk) | — | :material-check-circle:{ .green } |
| PacBio Iso-Seq (bulk) | — | :material-check-circle:{ .green } |
| MAGIC-seq | `magicseq` | :material-check-circle:{ .green } |
| 10x Visium | `visium` | :material-check-circle:{ .green } |
| 10x Genomics 3′ v3 | `10x3v3` | :material-progress-wrench:{ .orange } in development |

!!! warning "`10x3v3` is recognised but not yet implemented"
    The flag is accepted on the command line so that scripts do not break, but the run stops with an explanatory error instead of silently returning zero barcodes. If your 10x library has the same structure as Visium (primer + BC16 + UMI12 + polyT), you can process it with `-q visium` **provided you pass the whitelist through `-w`** — see the [Tutorials](https://weiwei4396.github.io/BroCOLI/Tutorials/).

## Supported reference data

* Reference genome should be provided in **FASTA** format.
* Reference annotation is not mandatory, but providing it will yield better results. Please provide it in **GTF** format.
* Chromosome names in the FASTA and the GTF must match exactly.

---

## Usage

> :memo: Go to [Tutorials](https://weiwei4396.github.io/BroCOLI/Tutorials/) for step-by-step workflows, and to [Output files](https://weiwei4396.github.io/BroCOLI/output/) for a description of every column BroCOLI writes.

## :tada: Citation

*Manuscript in preparation.* Please cite this repository in the meantime.

## Reference

1. [C++11 ThreadPool](https://github.com/progschj/ThreadPool)
2. Li H. [Minimap2](https://github.com/lh3/minimap2): pairwise alignment for nucleotide sequences[J]. Bioinformatics, 2018, 34(18): 3094-3100.
3. Marco-Sola S, Moure J C, Moreto M, Espinosa A. Fast gap-affine pairwise alignment using the wavefront algorithm[J]. Bioinformatics, 2021, 37(4): 456-463. [WFA2-lib](https://github.com/smarco/WFA2-lib)
4. Šošić M, Šikić M. [Edlib](https://github.com/Martinsos/edlib): a C/C++ library for fast, exact sequence alignment using edit distance[J]. Bioinformatics, 2017, 33(9): 1394-1395.
5. Cheng O, et al. [Flexiplex](https://github.com/DavidsonGroup/flexiplex): a versatile demultiplexer and search tool for omics data[J]. Bioinformatics, 2024, 40(3): btae102.

## Feedback and bug reports

If you come across any issues or have suggestions, please feel free to contact Wei Pan (weipan4396@gmail.com), or open an issue if you find bugs.
