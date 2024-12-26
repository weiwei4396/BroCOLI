# BroCOLI: Bron-Kerbosch calibrator of Long-read Isoform
BroCOLI is 
## Content
- [Requirements](#Requirements)
- [Installation](#Installation)
- [Supported sequencing data](#Supported-sequencing-data)
- [Supported reference data](#Supported-reference-data)
- [General usage](#General-usage)
    + [Bulk data](#Bulk-data)
        * [Step1 Mapping of the fastq files with minimap2](#Step1-Mapping-of-the-fastq-files-with-minimap2)
        * [Step2 Transcript identification and quantification](#Step2-Transcript-identification-and-quantification)
    + [Single cell data](#Single-cell-data)
- [Output files](#Output-files)
- [All Arguments](#All-Arguments)
- [Reference](#Reference)
- [Contact](#Contact)




## Requirements
**C++11** compatible compiler (e.g. **g++ 4.8** or later, clang++ 3.3 or later).

If you are a Linux user, you can check your g++ version with the following code.
```shell
g++ --version
```


## Installation
In order to compile the BroCOLI source in this GitHub repository the following steps can be taken:
```
git clone https://github.com/weiwei4396/BroCOLI.git
cd BroCOLI
sh build.sh
```
Once compiled, two executables ( **BroCOLI_bulk** and **BroCOLI_sc** ) will appear in the **BroCOLI folder**. You can use either the -h (--help) argument or the demo data to verify if the program runs successfully.
```
./BroCOLI_bulk -h
./BroCOLI_sc -h
```

## Supported sequencing data
BroCOLI support all kinds of long RNA data:
- ONT cDNA/ONT dRNA
- PacBio

Please align your reads to the reference genome and provide the resulting **sorted SAM file** to BroCOLI.

## Supported reference data
- Reference genome should be provided in **FASTA** format.
- Reference annotation is not required, but providing it will yield better results. If a reference gene annotation is available, please provide it in **GTF** format.

## General usage
### Bulk data
#### Step1 Mapping of the fastq files with minimap2
One alignment tool that could be considered is [Minimap2](https://github.com/lh3/minimap2).

The input **SAM files** need to be **sorted** by samtools before running BroCOLI.
```shell
samtools sort -o sorted.sam unsorted.sam
```
#### Step2 Transcript identification and quantification

### Single cell data


## Output files
1. After BroCOLI finishes processing the bulk data, a total of five files will be generated.
- `counts_transcript.txt`: Quantitative results of all transcripts contained in all samples.
    + Column 1 indicates Ensembl ID of each transcript.
    + Column 2 indicates Ensembl ID of each gene. Note that some novel transcripts are located in genes that are unclear, so BroCOLI represents its gene id as NA when output.
    + Column 3 to the end indicates the read count of transcript in each sample. The number of columns is equal to the number of samples.
- `counts_gene.txt`: Quantitative results of all genes contained in all samples. The result does not contain rows whose gene_id is NA.
    + Column 1 indicates Ensembl ID of each gene.
    + Column 2 to the end indicates the read count of gene in each sample. The number of columns is equal to the number of samples.
- `updated_annotitions.gtf`: An updated GTF annotation that includes both annotated and novel isoforms for the detected transcripts. 
    + The source column indicates for each isoform whether it is a `novel isoform` and `annotated isoform`.
    + The information for each isoform is recorded as one line, and the exon information as the next few lines.
- `compatible_isoform.tsv`: The result of each read assigned to the transcript in each sample file.
    + Column 1 (read_id) represents the read id presented in sam file. In the case of multiple sample files, column 5 specifies the file associated with each read.
    + Column 2 (category) specifies the classification of the isoform to which each read is associated.
    + Column 3 (isoform_id) provides the Ensembl ID of the isoform associated with each read.
    + Column 4 (gene_id) provides the Ensembl ID of the gene associated with each read.
    + Column 5 (file) contains a mapping between samples and a numerical index, which starts at 0 and ends at the number of samples minus one.
- `file_explain.txt`:
    + Column 1
    + Column 2 
2.After BroCOLI finishes processing the single cell data, a total of four files will be generated.
- `counts_transcript_index.txt`: A matrix of ***transcripts × cells***, with quantitative data obtained from each sample. The index denotes the serial number of the sample. The mapping between the sample and its number is provided in the `file_explain.txt` file or in BroCOLI's output history.
    + Column 1 indicates Ensembl ID of each transcript.
    + Column 2 indicates Ensembl ID of each gene. Note that some novel transcripts are located in genes that are unclear, so BroCOLI represents its gene id as NA when output.
    + Column 3 to the end indicates 

## All Arguments
```c++
./BroCOLI -h

```




## Reference
1. [C++11 ThreadPool](https://github.com/progschj/ThreadPool)
2. Li H. Minimap2: pairwise alignment for nucleotide sequences[J]. Bioinformatics, 2018, 34(18): 3094-3100.



## Contact
If you come across any issues or have suggestions, please feel free to contact Wei Pan (weipan4396@gmail.com), or open an issue if you find bugs.








