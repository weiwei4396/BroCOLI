---
icon: material/lightning-bolt
---

# BroCOLI Tutorials

Welcome to the BroCOLI Tutorials. This page walks through a complete run, from raw FASTQ to transcript- and gene-level counts, for bulk, single-cell and spatial data.

``` mermaid
%%{init: {'themeVariables': { 'fontSize': '20px' }}}%%
graph LR
  classDef input fill:#e3f2fd,stroke:#1565c0,stroke-width:2px,color:black;
  classDef process fill:#fff9c4,stroke:#fbc02d,stroke-width:2px,color:black;
  classDef result fill:#e8f5e9,stroke:#2e7d32,stroke-width:2px,color:black;

  A[raw fastq.gz]:::input
  C[matched fastq]:::process
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

## BroCOLI input files

To run BroCOLI, you should provide:

* FASTQ (FASTQ.gz) processed into a **sorted SAM** with [minimap2].
* Reference sequence in **FASTA** format.
* Optionally, a reference gene annotation in **GTF** format (recommended).

[minimap2]:https://github.com/lh3/minimap2

!!! note "Sorting is mandatory"
    The mapping **SAM files** need to be **sorted** by `samtools` before running BroCOLI. Chromosome names in the FASTA and the GTF must match.

---

## :book: BroCOLI General Usage

### Bulk Data

* [x] **Step 1** — Mapping of the FASTQ files with minimap2

```shell
minimap2 -ax splice -ub --secondary=no -t 20 ref.fasta raw.fastq.gz > raw.sam
samtools sort -@ 20 -o raw_sorted.sam raw.sam
```

??? info "Noisy cDNA data recommended parameter"
    For **noisy 1D cDNA Nanopore data** the developer of Minimap2 suggests adding **-k 14** and **-w 4**:
    ```shell
    minimap2 -ax splice -ub -k14 -w 4 --secondary=no -t 20 ref.fasta raw.fastq.gz > raw.sam
    ```

??? info "A BED file can be provided to assist the mapping"
    ```shell
    paftools.js gff2bed anno.gtf > junctions.bed
    minimap2 -ax splice -ub -k14 -w 4 --junc-bed junctions.bed --secondary=no -t 20 ref.fasta raw.fastq.gz > raw.sam
    ```

* [x] **Step 2** — Transcript identification and quantification

The `-s` parameter accepts three different shapes of input:

* **(i)** a single SAM file — give its absolute path;
* **(ii)** a directory containing all the sorted SAM files to be processed together;
* **(iii)** a TXT/TSV file listing one absolute path per line. The output column order will follow the order in this file.

```
(i)                     (ii)                    (iii)
input_reads.sam     ─── input_directory      ─── input.txt(.tsv)
                        ├── sample0.sam          ├── sample0.sam
                        └── sample1.sam          ├── sample1.sam
                                                 └── sample2.sam
```

```shell
./brocoli -t 8 -s sam_files_path -f ref.fasta -g GTF.gtf -o output_path
```

!!! tip "Which sample is which column?"
    Whichever shape you use, BroCOLI writes `file_explain.txt` in the output directory, mapping each column index back to its SAM file. See [Output files](https://weiwei4396.github.io/BroCOLI/output/).

---

### Single Cell and Spatial Data

* [x] **Step 1** — Barcode and UMI extraction with preBroCOLI

```shell
./prebrocoli -q visium -p 32 -w barcode_whitelist.txt -o demux/ raw.fastq.gz
```

Unlike version 1.x, **nothing is written to stdout** — every output goes into `--outdir`, so no shell redirection is needed and the log can never end up mixed into the data. The read file to process is the positional argument (or `-I`), and `-` reads from stdin.

The essential parameters:

| Flag | Meaning |
| :--- | :--- |
| `-q, --chemistry` | library type: `magicseq`, `visium` (`10x3v3` is in development) |
| `-o, --outdir` | output directory (required) |
| `-p, --thread` | worker threads |
| `-w, --whitelist` | cell barcode whitelist — see below |

The demultiplexed reads you will feed to minimap2 are `demux/preBroCOLI_matched.fastq`. The full option list is in [All arguments](#all-arguments).

#### Preparing the whitelist

`-w` takes a whitelist file whose **first column is the barcode**; anything after it is optional.

* **(i)** A TSV containing only barcodes — for example the filtered whitelist produced by CellRanger.
* **(ii)** A file containing **both barcodes and UMIs**, produced by the bundled `ext_bc_and_umi.py`. In this mode UMI correction is performed automatically, and the `umi_corrected` counter in the summary becomes meaningful.

!!! note "ext_bc_and_umi.py"
    ```python
    python ext_bc_and_umi.py --bam cellranger_processed.bam -f cellranger_filter_barcodes.tsv -o sc_bc_umi.txt
    ```

!!! warning "Use `-w`, not `-x`, for Visium"
    In the `visium` code path the final barcode assignment is guarded by a non-empty whitelist. If you pass only `-x/--barcodeX`, the whitelist stays empty, no barcode is ever assigned, and the run finishes successfully with zero matched reads. `-x/-y/-z` exist for the three segments of a MAGIC-seq barcode and are not a substitute for `-w`.

??? tip "Sanity-check on a subset first"
    Run the first ~10 000 reads and look at `Reads with barcode` in the summary before committing to the full file:
    ```shell
    zcat raw.fastq.gz | head -40000 | ./prebrocoli -q visium -w whitelist.txt -o testrun/ -
    ```
    A very low assignment rate usually means the wrong `-q`, a whitelist that does not match the library, or a `-f` threshold that is too strict.

??? info "You also can use [Flexiplex] for preprocessing"
    You can visit its GitHub page to learn more about its detailed usage.

    **First**, assign reads — short reads or single-cell long reads — to cellular barcodes
    ```shell
    flexiplex -d 10x3v3 -p 20 -k cellRangerbarcodes.tsv raw.fastq > new_reads.fastq
    ```
    **Second**, mapping.
    ```shell
    minimap2 -ax splice -ub -k14 -w 4 --secondary=no -t 20 ref.fasta new_reads.fastq > new_reads.sam
    samtools sort -@ 20 -o new_reads_sorted.sam new_reads.sam
    ```

[Flexiplex]:https://github.com/DavidsonGroup/flexiplex

??? info "You also can use [Sicelore-2.1] for preprocessing"
    You can visit its GitHub page to learn more about its detailed usage. Before you run sicelore, you need to set up the required JAVA environment for it.

    **First**, scan Nanopore reads — assign cell barcodes.
    ```shell
    java -jar -Xmx80g <path>/NanoporeBC_UMI_finder-2.1.jar scanfastq -d <directory to start recursive search for fastq files> -o outPutDirectory --bcEditDistance 1 --cellRangerBCs cellRangerbarcodes.tsv
    ```
    The `--cellRangerBCs` parameter is optional. If Illumina data is available, a TSV file containing cell barcodes (e.g. from Cell Ranger) can be provided, which will improve the accuracy of barcode identification.

    **Second**, mapping.
    ```shell
    minimap2 -ax splice -ub -k14 -w 4 --junc-bed junctions.bed --sam-hit-only --secondary=no -t 20 ref.fasta <fastq.gz path> > raw.sam
    samtools view -bS -@ 20 raw.sam > raw.bam
    samtools sort -@ 20 -o raw_sorted.bam raw.bam
    samtools index -@ 20 raw_sorted.bam raw_sorted_index
    ```

    **Third**, UMI assignment.
    ```shell
    java -jar -Xmx80g <path>/NanoporeBC_UMI_finder-2.1.jar assignumis --inFileNanopore raw_sorted.bam --outfile raw_sorted_umi.bam --ONTgene GE --annotationFile GTF.gtf
    ```
    The output BAM generated by the cell barcode and UMI assignment is converted to a SAM file for BroCOLI's input.

[Sicelore-2.1]:https://github.com/ucagenomix/sicelore-2.1

* [x] **Step 2** — Mapping of the demultiplexed FASTQ with minimap2

The processing is identical to Step 1 in the bulk workflow.

```shell
minimap2 -ax splice -ub --secondary=no -t 20 ref.fasta demux/preBroCOLI_matched.fastq > preBroCOLI.sam
samtools sort -@ 20 -o preBroCOLI_sorted.sam preBroCOLI.sam
```

!!! note "Keep the barcode tags"
    preBroCOLI rewrites each read ID as `barcode_umi#readid` and appends `CB:Z:` / `UB:Z:` tags. minimap2 carries the comment field through, which is how BroCOLI recovers the cell of origin — do not strip it with tools that rewrite read names.

* [x] **Step 3** — Transcript identification and quantification

The input data is specified exactly as in the bulk workflow.

```
(i)                     (ii)                    (iii)
input_reads.sam     ─── input_directory      ─── input.txt(.tsv)
                        ├── sample0.sam          ├── sample0.sam
                        └── sample1.sam          ├── sample1.sam
                                                 └── sample2.sam
```

```shell
./brocoli -t 8 -s sam_files_path -f ref.fasta -g GTF.gtf -o output_path
```

---

## :eyes: Examples

### Simple test

* **Bulk: SIRV4 dataset**

```shell
./brocoli -t 1 -s example/example_SIRV.sam -g example/example_SIRV.gtf -f example/example_SIRV.fasta -o TestResult
```

* **Single cell**

```shell
./brocoli -t 8 -s MouseBrain_sorted.sam -g mouse.gtf -f mouse.fasta -o TestResult
```

### End-to-end spatial run

```shell
conda activate brocoli

# 1. barcode + UMI extraction
./prebrocoli -q visium -p 32 -w barcode_whitelist.txt -o demux/ raw.fastq.gz

# 2. mapping
minimap2 -ax splice -ub --secondary=no -t 20 ref.fasta \
         demux/preBroCOLI_matched.fastq > visium.sam
samtools sort -@ 20 -o visium_sorted.sam visium.sam

# 3. quantification
./brocoli -t 20 -s visium_sorted.sam -f ref.fasta -g anno.gtf -o result/
```

---

## All Arguments

```shell
./brocoli -h
./prebrocoli -h
```

### `brocoli`

```
Arguments:
-s, --sam
      SAM file path. We recommend using absolute paths. If you have a single file, you can
      directly provide its absolute path. If you have multiple files, you can specify the path
      to a folder that contains all the sorted SAM files you want to process. (required)

-f, --fasta
      FASTA file path. FASTA file requires the chromosome names to match the GTF file. (required)

-o, --output
      output folder path. (required)

-g, --gtf
      input annotation file in GTF format. (optional, Recommendation provided)

-n, --support
      min perfect read count for all splice junctions of novel isoform. (optional, default: 2)

-j, --SJDistance
      the minimum distance determined as intron. (optional, default: 18)

-e, --single_exon_boundary
      belongs to the isoform scope of a single exon. (optional, default: 60)

-d, --graph_distance
      the distance threshold for constructing the isoform candidate distance graph. (optional, default: 60)

-t, --thread
      thread number. (optional, default: 8)

-h, --help
      show this help information.
```

### `prebrocoli`

```
Usage: prebrocoli -q <chemistry> -w <whitelist> -o <outdir> <reads.fastq[.gz]>

  reads may be .fastq, .fastq.gz, .fasta, .fasta.gz, or - for stdin.
```

**Required**

| Flag | Description |
| :--- | :--- |
| `-q, --chemistry` | `magicseq` \| `visium` (`10x3v3` not implemented yet) |
| `-o, --outdir` | output directory |

**Barcode source (one of)**

| Flag | Description |
| :--- | :--- |
| `-w, --whitelist` | whitelist file: barcode in column 1, optional UMI after it |
| `-x, --barcodeX` | barcode list (MAGIC-seq: first segment) |
| `-y, --barcodeY` | MAGIC-seq only, second segment |
| `-z, --barcodeZ` | MAGIC-seq only, third segment |

**Common**

| Flag | Default | Description |
| :--- | :--- | :--- |
| `-I, --reads` | positional arg | reads file (alternative to the positional argument) |
| `-n, --prefix` | `preBroCOLI` | output filename prefix |
| `-p, --thread` | `1` | worker threads |
| `-f, --flank_editd` | `20` (magicseq) / `8` (visium) | max adapter edit distance |
| `-i, --trim` | `true` | rewrite the read ID and trim adapters |
| `-s, --split` | `false` | write one file per barcode |
| `-c, --chimeric` | `false` | mark chimeric reads with `_C` |
| `--keep_unmatched` | off | also write reads with no barcode |
| `--batch` | `2000` | reads per thread per batch |
| `-h, --help` | | show this help information |

!!! tip "Tuning notes"
    * **`-f`** is the edit-distance budget for the *adapter/flank*, not for the barcode. Raising it recovers reads with a degraded primer at the cost of runtime and a few more false anchors; `adapter_not_found` in the summary tells you whether it is the limiting step.
    * **`-p` × `--batch`** is how many reads are resident at once. The default 2000 is a good balance; lower it if memory is tight, raise it if the threads are starving on very short reads.
    * **`-s/--split`** is ignored automatically when there are more than 50 barcodes — a Visium or 10x run would otherwise open tens of thousands of files.

---

## :question: FAQ

??? question "The run finished successfully but no barcodes were found"
    Three usual causes, in order of frequency: the whitelist was passed with `-x` instead of `-w` (see the warning above); the chemistry does not match the library; or the reads are already demultiplexed by another tool and no longer carry the adapter. Check `adapter_not_found` and `whitelist_no_match` in `<prefix>_summary.txt` to tell them apart.

??? question "`Records written` is larger than `Reads with barcode`"
    That is expected. `Records written` counts **records in the matched FASTQ**, and one read contributes one record per (barcode × strand). A read carrying two barcodes produces two records, suffixed `1of2` and `2of2`. The number of *reads* is `Reads with barcode`.

??? question "Can I process a 10x 3′ v3 library today?"
    If its structure matches Visium (primer + BC16 + UMI12 + polyT) you can run `-q visium` with the 10x whitelist passed through `-w`. Validate on a small subset first: the candidate count and edit-distance thresholds in that code path were tuned for Visium's smaller whitelist, so the assignment rate may be lower than a native implementation would give.

??? question "Do I have to sort the SAM file?"
    Yes. BroCOLI walks the alignments in coordinate order and will not produce correct clusters otherwise.
