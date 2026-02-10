---
icon: simple/markdown
---

# Introduction to BroCOLI Output Files

## Bulk output files


* [x] 1. file_explain.txt

This file establishes a mapping between samples and a numerical index, ranging from 0 to (number of samples − 1).

```
- Column 1 contains the index assigned to each sample by BroCOLI.
- Column 2 contains the absolute path to the corresponding SAM file for that sample.
```


* [x] 2. **counts_transcript.txt**

This file contains the quantitative read counts for all transcripts across all samples.

```
- Column 1: Ensembl transcript ID.
- Column 2: Ensembl gene ID of the corresponding gene. For novel transcripts with unclear or unmapped gene associations, BroCOLI outputs NA in this column.
- Columns 3 to the end: Read counts of the transcript in each sample (one column per sample; the total number of columns equals the number of samples).
```


* [x] 3. **counts_gene.txt**

This file contains the quantitative read counts for all genes across all samples.

```
- Column 1: Ensembl gene ID.
- Columns 2 to the end: Read counts of the gene in each sample (one column per sample; the total number of columns equals the number of samples).
```


* [x] 4. updated_annotitions.gtf

This is an updated GTF annotation file that incorporates both known (annotated) and novel isoforms for the detected transcripts.

```
- The source column indicates the origin of each isoform (novel for newly discovered isoforms or annotated for known isoforms).
- Each isoform is described on a single line containing its feature information, followed by one or more subsequent lines detailing its exon coordinates.
```

* [x] 5. compatible_isoform.tsv

This file reports the assignment of each read to a specific isoform across all sample files.

```
- Column 1 (read_id): The read identifier as it appears in the original SAM file.
- Column 2 (category): Classification of the read–isoform match. BroCOLI categorizes reads into four types:
  * FSM: Full splice match (complete match to a known isoform).
  * ISM: Incomplete splice match (partial match to a known isoform).
  * SE: The isoform consists of a single exon.

Column 3 (isoform_id): Ensembl transcript ID of the assigned isoform.
Column 4 (gene_id): Ensembl gene ID associated with the assigned isoform.
Column 5 (file): Numerical index of the sample from which the read originates. The mapping between indices and actual sample files is provided in file_explain.txt.
```





## single cell and spatial output files
* [x] 1. file_explain.txt

This file establishes a mapping between samples and a numerical index, ranging from 0 to (number of samples − 1).

```
- Column 1 contains the index assigned to each sample by BroCOLI.
- Column 2 contains the absolute path to the corresponding SAM file for that sample.
```





