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
MEIGA-PAV takes as input 6 mandatory arguments:

   1. VCF: Input VCF file containing sequence-resolved structural variation calls. 
   2. consensus: Fasta file containing consensus sequences for retrotransposon subfamilies (Alu, L1 and SVA).
   3. reference: Fasta file for the human genome reference. Use same chromosome labeling nomenclature as for the input VCF. BWA-mem index need to be located in the same folder
   4. mt_reference: Fasta file for the mitochrondrial genome reference. 
   5. sampleID: Output VCF file will be named accordingly.
   6. outDir: Output directory. 

## Execution
python MEIGA-PAV.py freeze3.sv.alt.vcf consensusDb.fa GCA_000001405.15_GRCh38_no_alt_analysis_set_no_chr_MT.fna chrM.fa 'PAV_MEI_MEIGA.freeze3' -o freeze3

File availability:
  * Example VCF file for testing purposes can be downloaded from the 1000 genomes data portal (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/release/v1.0/integrated_callset/)
  * The database of consensus retrotransposon sequences is provided in the /execution folder. 
  * Reference genome sequences can be downloaded from the UCSC browser (ftp://hgdownload.soe.ucsc.edu/goldenPath).

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

An example of VCF output (PAV_MEI_MEIGA.freeze3.vcf) is provided in the /execution folder. This results from the annotation of the input VCF (freeze3.sv.alt.vcf) included in the folder as well as an example 

## License
MEIGA-PAV is distributed under GPL-3.0 License

## Cite MEIGA-PAV
Ebert et al., “Haplotype-Resolved Diverse Human Genomes and Integrated Analysis of Structural Variation,” Science, February 25, 2021, eabf7117, https://doi.org/10.1126/science.abf7117 (PMID: 33632895).

## Contact
Please open a case on the Github page for problems.

You may also contact Bernardo Rodriguez-Martin directly (e-mail omitted to reduce SPAM). MEIGA-PAV was developed in the lab of Dr. José Tubio at CiMUS and Dr. Jan Korbel at EMBL.

