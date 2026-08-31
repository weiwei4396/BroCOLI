---
icon: simple/markdown
---

# Introduction to BroCOLI Output Files

Every file BroCOLI and preBroCOLI produce is plain text and tab-separated, so they can be read directly with `awk`, `pandas` or R.

| Stage | Directory | Files |
| :--- | :--- | :--- |
| `prebrocoli` (single-cell / spatial only) | `--outdir` | `_matched.fastq`, `_unmatched.fastq`, `_reads_barcodes.tsv`, `_barcode_counts.tsv`, `_summary.txt` |
| `brocoli` (all workflows) | `-o` | `file_explain.txt`, `counts_transcript.txt`, `counts_gene.txt`, `updated_annotitions.gtf`, `compatible_isoform.tsv` |

---

## preBroCOLI output files

All five files are prefixed with `-n/--prefix` (default `preBroCOLI`) and written under `-o/--outdir`. Nothing is written to stdout, so the log can never be mixed into the data.

* [x] 1. **`<prefix>_matched.fastq`**

The reads that received a barcode, ready for minimap2. By default (`-i true`) everything up to and including the adapter is trimmed away and the read ID is rewritten:

```
@ACGTACGTACGTACGT_TTGCAAGGCTAC#read0001_+1of1	CB:Z:ACGTACGTACGTACGT	UB:Z:TTGCAAGGCTAC
TCAGGTTAACGGCTTACAGGATCCTTAGGCAT...
+
IIIIIIIIIIHHHHIIIIIIIIIIIIHHGGFF...
```

The header is built as `barcode _ UMI # original_read_id strand n of N [_C]`:

```
- barcode / UMI : the assigned cell barcode and its (corrected) UMI.
- strand        : _+ for the forward strand, _- for the reverse complement.
- n of N        : this record is the n-th of N barcodes found on the read.
- _C            : appended only with -c true, marking a chimeric read
                  (barcodes found on both strands).
- CB:Z: / UB:Z: : SAM-style tags in the comment field, carried through by
                  minimap2 so that BroCOLI can recover the cell of origin.
```

!!! note "One read can produce more than one record"
    A read with two barcodes yields two records (`1of2`, `2of2`), and a chimeric read contributes on both strands. This is why `Records written` in the summary is normally larger than `Reads with barcode`. With `-i false` the read is emitted untrimmed under its original ID instead.

* [x] 2. **`<prefix>_unmatched.fastq`**

Written **only** when `--keep_unmatched` is given. Reads for which no barcode could be called, kept verbatim so the loss can be audited.

* [x] 3. **`<prefix>_reads_barcodes.tsv`**

One line per (read × assigned barcode), with a header row.

| Column | Name | Description |
| :---: | :--- | :--- |
| 1 | `Read` | read ID as it appears in the input FASTQ |
| 2 | `CellBarcode` | the assigned barcode, after whitelist correction |
| 3 | `FlankEditDist` | edit distance between the expected adapter/flank and the read; governed by `-f/--flank_editd` |
| 4 | `BarcodeEditDist` | edit distance between the observed barcode and the whitelist entry that was called |
| 5 | `UMI` | the UMI, corrected against the whitelist when the whitelist file carries UMIs |

```
Read	CellBarcode	FlankEditDist	BarcodeEditDist	UMI
read0001	ACGTACGTACGTACGT	3	0	TTGCAAGGCTAC
read0002	GGTTCCAAGGTTCCAA	5	1	ACCTGGTTACGA
```

!!! tip "Reading this file"
    Only calls with `BarcodeEditDist` ≤ 2 are kept, so this column is a quality signal rather than a filter you still need to apply. A read ID appearing twice means two barcodes were found on it — either a genuine chimera or a concatemer. Count unique reads with:
    ```shell
    awk 'NR>1{print $1}' preBroCOLI_reads_barcodes.tsv | sort -u | wc -l
    ```

* [x] 4. **`<prefix>_barcode_counts.tsv`**

Barcode → read count, sorted by count descending and then by barcode ascending. This is the file to use for a knee plot or for choosing which cells to keep.

```
Barcode	Reads
ACGTACGTACGTACGT	18342
GGTTCCAAGGTTCCAA	17905
...
```

* [x] 5. **`<prefix>_summary.txt`**

The same block that is printed to the terminal at the end of the run, saved for the record. It is the first place to look when the assignment rate is not what you expected.

```
-----------------------------------------------------
Summary:
Chemistry:             visium
Reads in:              100000
Reads with barcode:    86525 (86.53%)
Reads without barcode: 13475 (13.47%)  [discarded, use --keep_unmatched to save]
Reads chimeric:        842
Reads >1 barcode:      1217
Records written:       87788
-- where barcode searches failed (per search attempt, both strands) --
read_too_short:        12
adapter_not_found:     9856
bc1_unresolved:        2104
...
```

**Read accounting**

| Field | Meaning |
| :--- | :--- |
| `Reads in` | reads read from the input file |
| `Reads with barcode` | reads for which at least one barcode was called |
| `Reads without barcode` | the remainder; discarded unless `--keep_unmatched` |
| `Reads chimeric` | barcodes were found on **both** strands of the same read |
| `Reads >1 barcode` | more than one barcode was found on the read |
| `Records written` | **records** in the matched FASTQ, not reads — see the note above |

**Where barcode searches failed** *(counted per search attempt, i.e. each read is attempted on both strands)*

| Field | Meaning |
| :--- | :--- |
| `read_too_short` | the read is shorter than the search pattern itself |
| `adapter_not_found` | the adapter/flank could not be aligned within `-f`; raise `-f` if this dominates |
| `bc1_unresolved` | the first barcode segment had no candidate (the only segment in `visium`) |
| `bc2_unresolved` / `bc3_unresolved` | second and third MAGIC-seq segments; always 0 for `visium` |
| `whitelist_no_match` | candidates existed but none scored well enough against the whitelist |
| `no_confident_call` | a barcode was found but the call was ambiguous |
| `barcode_editd_high` | the best call exceeded the edit-distance ceiling of 2 |

**UMI**

| Field | Meaning |
| :--- | :--- |
| `umi_exact` | the observed UMI was already in the whitelist for that barcode |
| `umi_corrected` | the UMI was rescued to a whitelist UMI within 3 edits |
| `umi_uncorrected` | no whitelist UMI was close enough; the raw sequence is kept |

!!! note
    The three UMI counters are only informative when the whitelist file carries UMIs — for example one produced by `ext_bc_and_umi.py`. With a barcode-only whitelist every UMI is passed through unchanged.

**CPU time** — the last block sums time across threads, so it exceeds the wall time on a multi-threaded run. It is there to show which stage dominates (`adapter search` almost always does).

---

## Bulk output files

* [x] 1. **`file_explain.txt`**

This file establishes a mapping between samples and a numerical index, ranging from 0 to (number of samples − 1).

```
- Column 1 contains the index assigned to each sample by BroCOLI.
- Column 2 contains the absolute path to the corresponding SAM file for that sample.
```

* [x] 2. **`counts_transcript.txt`**

This file contains the quantitative read counts for all transcripts across all samples.

```
- Column 1: Ensembl transcript ID.
- Column 2: Ensembl gene ID of the corresponding gene. For novel transcripts with unclear
            or unmapped gene associations, BroCOLI outputs NA in this column.
- Columns 3 to the end: Read counts of the transcript in each sample (one column per sample;
            the total number of columns equals the number of samples).
```

!!! tip "Column order"
    Columns follow the index in `file_explain.txt`, which for input shape **(iii)** is the order of the lines in your `input.txt`.

* [x] 3. **`counts_gene.txt`**

This file contains the quantitative read counts for all genes across all samples.

```
- Column 1: Ensembl gene ID.
- Columns 2 to the end: Read counts of the gene in each sample (one column per sample;
            the total number of columns equals the number of samples).
```

* [x] 4. **`updated_annotitions.gtf`**

An updated GTF annotation that incorporates both known (annotated) and novel isoforms for the detected transcripts. It can be used directly as the annotation for a re-run or for visualisation in IGV.

```
- The source column indicates the origin of each isoform (novel for newly discovered
  isoforms, or annotated for known isoforms).
- Each isoform is described on a single line containing its feature information, followed
  by one or more subsequent lines detailing its exon coordinates.
```

* [x] 5. **`compatible_isoform.tsv`**

This file reports the assignment of each read to a specific isoform across all sample files — the read-level evidence behind every number in the count tables.

| Column | Name | Description |
| :---: | :--- | :--- |
| 1 | `read_id` | the read identifier as it appears in the original SAM file |
| 2 | `category` | classification of the read–isoform match, see below |
| 3 | `isoform_id` | Ensembl transcript ID of the assigned isoform |
| 4 | `gene_id` | Ensembl gene ID associated with the assigned isoform |
| 5 | `file` | numerical index of the sample the read came from, as listed in `file_explain.txt` |

Categories:

| Code | Meaning |
| :--- | :--- |
| `FSM` | Full splice match — complete match to a known isoform |
| `ISM` | Incomplete splice match — partial match to a known isoform |
| `SE` | The isoform consists of a single exon |

---

## Single cell and spatial output files

The single-cell and spatial run produces the same set of files as the bulk run. The difference is the unit of a column: **counts are reported per cell barcode instead of per sample**, recovered from the `CB:Z:` tag that preBroCOLI wrote into each read.

* [x] 1. **`file_explain.txt`**

This file establishes a mapping between samples and a numerical index, ranging from 0 to (number of samples − 1).

```
- Column 1 contains the index assigned to each sample by BroCOLI.
- Column 2 contains the absolute path to the corresponding SAM file for that sample.
```

* [x] 2. **`counts_transcript.txt` / `counts_gene.txt`**

The count matrices, with one row per feature (transcript or gene) and one column per cell barcode. UMI-aware collapsing is applied so that multiple reads from the same (barcode, UMI) molecule are counted once.

* [x] 3. **`compatible_isoform.tsv`**

As in the bulk workflow, with the read identifier still carrying the `barcode_umi#readid` prefix written by preBroCOLI, so any assignment can be traced back to its cell.
