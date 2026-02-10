---
icon: lucide/rocket
---

# Welcome to BroCOLI’s documentation!

## Introduction
BroCOLI (Bron-Kerbosch calibrator of Long-read Isoform) leverages BK algorithms for transcript identification and quantification from long-read RNA-Seq data, supporting bulk, single-cell and spatial applications, while maintaining low memory usage and fast performance for large-scale datasets.

## Requirements
BroCOLI requires a **C++11** compatible compiler (e.g., g++ 4.8 or later).

## :rocket: Installation
BroCOLI is available on [`BroCOLI GitHub`][BroCOLInew].
  [BroCOLInew]:https://github.com/gyjames/BroCOLI

BroCOLI can be installed via the source code:
``` shell title="download and install"
git clone https://github.com/gyjames/BroCOLI.git
cd BroCOLI
sh build.sh
```
Once compiled, two executables ( BroCOLI_bulk and BroCOLI_sc ) will appear in the BroCOLI folder. You can use either the -h (--help) argument or the demo data to verify if the program runs successfully.
``` shell title="help"
./BroCOLI_bulk -h
./BroCOLI_sc -h
```


## :sparkles: Supported sequencing data
* ONT cDNA/ONT dRNA
* PacBio
* 10x3v3
* 10X Visium
* MAGIC-seq


## Supported reference data
* Reference genome should be provided in **FASTA** format.
* Reference annotation is not mandatory, but providing it will yield better results. Please provide it in **GTF** format.

## Usage
> :memo: Go to [Tutorials](http://localhost:8000/Tutorials/)



## :tada: Citation



## Reference

1. [C++11 ThreadPool](https://github.com/progschj/ThreadPool)
2. Li H. [Minimap2](https://github.com/lh3/minimap2): pairwise alignment for nucleotide sequences[J]. Bioinformatics, 2018, 34(18): 3094-3100.
3. [Flexiplex](https://github.com/DavidsonGroup/flexiplex) 


## Feedback and bug reports

If you come across any issues or have suggestions, please feel free to contact Wei Pan (weipan4396@gmail.com), or open an issue if you find bugs.
