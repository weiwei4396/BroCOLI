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




## All Arguments
```c++
./BroCOLI -h

```




## Reference
1. [C++11 ThreadPool](https://github.com/progschj/ThreadPool)
2. Li H. Minimap2: pairwise alignment for nucleotide sequences[J]. Bioinformatics, 2018, 34(18): 3094-3100.



## Contact
If you come across any issues or have suggestions, please feel free to contact Wei Pan. weipan4396@gmail.com

Please open an issue if you find bugs.








