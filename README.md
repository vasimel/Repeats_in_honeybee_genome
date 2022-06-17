# Analysis of repertoire and expression of repeated elements in Apis Mellifera (honeybee) genome
## Introduction
The honeybee (A. mellifera) is an eusocial hymenopteran insect; its families are divided into castes, each with a specific role in the hive. Only the mothers are capable of reproduction, the rest of the females are worker bees, i.e. they take care of the offspring and forage. 

Approximately 11% of the honeybee genome consists of repeated elements - tandem repeats and transposons; they play a central role in chromosome stability, the cell cycle, regulation of gene expression, and are an important substrate for genome evolution. Although the processes and genes of caste determination and the genes regulating this process in honeybees are well studied, the role of non-coding elements is not fully understood. 

The aim of the project is to determine the role of repeated elements in the caste development.

## Methods
### Repeated elements de-novo identification
Apis Mellifera [representative genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_003254395.2/) was analysed to identify repeat models using [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler).
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
gunzip GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
conda activate repeatmodeler
nohup BuildDatabase -name honeybee GCF_003254395.2_Amel_HAv3.1_genomic.fna >& building_database.out &
nohup RepeatModeler -database honeybee >& run.out &
```
To classify the repeat models a merged Dfam 3.5 and RepBase library was constructed.
### Repeat masking
[RepeatMasker](https://www.repeatmasker.org/) is a tool that screens DNA for interspersed repeats and low complexity sequences and masks repeated regions. We used it with an option -a to obtain an .align file to count divergence between particular copies and consensus and to build a Kimura plot.
```
nohup RepeatMasker -a -gff -dir ../repeatmasker/ -lib honeybee-families.fa honeybee.fna >& ./masker.out &
```
### Kimura plot
To obtain table of divergence % we ran a [Perl script](https://github.com/4ureliek/Parsing-RepeatMasker-Outputs) for parsing RepeatMasker .align output file and built Kimura plots of honeybee repeat families.
![image](https://user-images.githubusercontent.com/92578463/174393478-2c0c474f-3e63-4e6b-abf9-ef7588155359.png)

