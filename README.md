# HomologSearch
R code for identic model and nonmodel organism genes searching.

The HomologSearch serves for searching similarities between two sequences (protein or nucleotide) in a fasta format. The evaluation is computed using the BLAST algorithm and output is a ta-ble of homologous locus tags and genes (if available). The HomologSearch is available in the GitHub repository under the name ‘JanaMus/HomologSearch’ and can be installed as a library into the R language. 

### DEPENDENCIES 
1. R (free from https://www.r-project.org/)
2. BLAST+ (free from https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download, NCBI/BLAST+ or from the NCBI website) along with the R path set to its folder
3. R packages: 
a) devtools: 
```shell
install.packages("devtools")
library(devtools)
```
b) Biostrings
```shell
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
library(Biostrings)
```
c) rBLAST
```shell
install_github("mhahsler/rBLAST")
library(rBLAST)
```
### INSTALATION
```shell
install_github("JanaMus/HomologSearch")
library(HomologSearch)
```
### USAGE
To launch the package three inputs in a string format are required. The first one is the sequence of the model organism, which will be used to create the local database, labelled as ‘model_organism’. Second input is the sequence of the non-model organism and will be used for the comparison with the local database labelled as ‘nonmodel_organism’ Both sequences must be uploaded in the fasta format. The last one is one of two values according to the type of uploaded sequences: nucleotide or protein. This input is labelled as ‘fasta_type’. 

### EXAMPLE USAGE
Example sequences were downloaded from the NCBI GenBank database. The model organism is the _Escherichia coli_ ATCC 8739, available from https://www.ncbi.nlm.nih.gov/nuccore/CP000946.1. The selected non-model organism is the _Clostridium beijerinckii_ NRRL B-598 (https://www.ncbi.nlm.nih.gov/nuccore/CP011966.3). Sequences were downloaded as follows. Send to: Coding Sequences, Format: FASTA protein and stored into the same folder as the BLAST+ files. 

1.	Set the path to the folder that contains BLAST + and both sequences: 
```shell
setwd(“docs/Rpackages/HomologSearch/BLAST+sequences”)
```
2.	Load the library:
```shell
library(HomologSearch)
```
3.	Run the HomologSearch(model_organism, nonmodel_organism, fasta_type):
```shell
HomologSearch(‘model_organism.fasta’, ‘nonmodel_organism.fasta’, ‘protein’)
```
4.	The output file containing names of homologous locus tags or gene IDs "homologs.csv" is stored in the folder into path was set in the first step.

