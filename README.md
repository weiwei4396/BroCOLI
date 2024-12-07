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


#### Step2 Transcript identification and quantification

### Single cell data





## Reference
1. c++



## Contact
Wei Pan. weipan4396@gmail.com

Please open an issue if you find bugs.








