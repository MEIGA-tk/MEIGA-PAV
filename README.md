## MEIGA-PAV
MEIGA-PAV is a computational method for the annotation of mobile element insertions (MEI) in a VCF containing sequence-resolved structural variation (SV) calls. Although it has been designed to process SVs derived from the Phased Assembly Variant Caller (PAV; https://github.com/EichlerLab/pav), it should be compatible with the output of any long-read variant callers, as long as the sequence for insertion events is included in the VCF. 

MEIGA-PAV was developed for the Human Genome Structural Variation Consortium (HGSVC)

Ebert et al., “Haplotype-Resolved Diverse Human Genomes and Integrated Analysis of Structural Variation,” Science, February 25, 2021, eabf7117, https://doi.org/10.1126/science.abf7117.

## Download 
Two different ways:

* Go to the releases tab and download the latest release. 

* Clone the git repository in case you want the latest version of the code:

```
# Move into the folder in which you want to clone the repositoy.
$ cd ~/apps
# Clone it.
$ git clone https://github.com/Chimera-tools/ChimPipe.git
```

MEIGA-PAV does not require any further installation step. It is written in Python and can be run as a standalone application on diverse Linux systems. 

## Requirements
1. Hardware:

    * 64-bits CPU

2. Software:

    * 64-bit Linux System
    * Python v3.6 or higher
    * paftools v.r755 (https://github.com/lh3/minimap2/tree/master/misc)
    * bwa-mem v0.7.17 (https://github.com/lh3/bwa)
    * minimap2 v2.1 (https://github.com/lh3/minimap2)

3. Python libraries 
    * pysam 
    * cigar
    * numpy
    * itertools
    * Biopython
   
## Input
script=/Users/brodriguez/Research/Projects/MEIGA-tk/MEIGA-PAV/MEIGA-PAV.py
vcf=/Users/brodriguez/Research/Projects/HGSVC2/Data/SV/Freeze4/freeze4.sv-alt.vcf
consensus=/Users/brodriguez/Research/Projects/MEIGA/MEIGA/databases/H.Sapiens/hg19/consensusDb.fa
reference=/Users/brodriguez/Research/References/Genomes/H.sapiens/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set_no_chr_MT.fna
mtGenome=/Users/brodriguez/Research/References/Genomes/H.sapiens/GRCh38/chrM.fa
sampleId='PAV_MEI_MEIGA.freeze4'
outDir=/Users/brodriguez/Research/Projects/HGSVC2/Analysis/Main/MEIGA_calls_pav/freeze4/
mkdir $outDir

source activate py36
time python $script $vcf $consensus $reference $mtGenome $sampleId -o $outDi

## Output
MEIGA-PAV produces as output a standard VCF file containing the subset of insertions from the input VCF corresponding to MEI. MEI annotation information is included in the following additional INFO fields:

    * ITYPE: Type of insertion. Solo retrotransposon insertion, partnered or orphan transduction. 
    * FAM: Repeat family (Alu, L1 or SVA)
    * 3PRIME: Flag indicating a 3-prime partnered transduction
    * 5PRIME: Flag indicating a 5-prime partnered transduction
    * CYTOID: Source element cytoband identifier (for transductions)   
    * RETRO_LEN: Length for the inserted retrotransposon
    * TRUNCATION_5_LEN: Length for 5-prime truncation 
    * TRUNCATION_3_LEN: Length for 3-prime truncation
    * INVERSION_LEN: Length for 5-prime inversion
    * RETRO_COORD: Coordinates with respect to the consensus L1 for the retrotransposon insertion
    * IS_FULL: Flag indicating a full length mobile element insertion  
    * ORF1: Flag indicating the existance of a funtional ORF1 for L1s
    * ORF2: Flag indicating the existance of a funtional ORF2 for L1s  
    * COMPETENT: Flag indicating that the insertion is a potential competent full-length L1 with intact ORFs    
    * TDCOORD_5PRIME: Genomic coordinates for 5-prime transduction
    * TDCOORD_3PRIME: Genomic coordinates for 3-prime transduction
    * TDLEN_5PRIME: Length of 5-prime transduction 
    * TDLEN_3PRIME: Length of 3-prime transduction 
    * STRAND: Insertion DNA strand (+ or -)

## License
MEIGA-PAV is distributed under GPL-3.0 License

## Cite MEIGA-PAV
Ebert et al., “Haplotype-Resolved Diverse Human Genomes and Integrated Analysis of Structural Variation,” Science, February 25, 2021, eabf7117, https://doi.org/10.1126/science.abf7117 (PMID: 33632895).

## Contact
Please open a case on the Github page for problems.

You may also contact Bernardo Rodriguez-Martin directly (e-mail omitted to reduce SPAM). MEIGA-PAV was developed in the lab of Dr. José Tubio at CiMUS and Dr. Jan Korbel at EMBL.

