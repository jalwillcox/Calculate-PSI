Calculate Percent Spliced-In (PSI)
==================================

This script compares the PSI between case and control samples for RNA-Seq data.

PSI is calculated on an individual basis and then averaged, NOT by comparing all cases merged vs. all controls merged.

If you just want to calculate PSI (no case/control comparison), label all samples as "case" or exclude the -c flag.

Required Programs
-----------------

This script uses the following programs; alternate versions may also work:

  * [bedtools/2.27.1](https://github.com/arq5x/bedtools2)
  * [gcc/6.2.0](https://linuxfromscratch.org/blfs/view/7.10/general/gcc.html)
  * [bowtie2/2.2.9](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  * [tophat/2.1.1](https://ccb.jhu.edu/software/tophat/index.shtml)
  * [R/3.6.1](https://www.r-project.org/)
  * [samtools/1.3.1](http://www.htslib.org/doc/1.3.1/samtools.html)

The following Rscripts are included with this pipeline:

  * **configure-psi-byIndiv.R** *combine individual PSI output and calculate average PSI across samples* 
  * **comp-PSI.R** *compare cases to control (includes uncorrected p-values)*
  * **plot-psi.R** *generates a plot PSI plot for the transcript with the most exons, comparing cases and controls; the R package ggplot2 is required*
  * **plot-psi-noComp.R** *generates a plot for the transcript with the most exons, WITHOUT a comparison; the R package ggplot2 is required*


Flags
-----

> -b (arg, required)     a file with col1=sample_id, col2=bam, and col2=case/control (example BAM list below)<br />
> -e (arg, required)     G/T - \"G\" to calculate for a gene or \"T\" to calculate for a transcript<br />
> -g (arg, required)     gene/transcript name<br />
><br />
> -c                     compare cases and controls (must be specified in bamlist)<br />
> -f (arg)               the full path to the GTF file (default: ./hg38.gencode.v27.primary_assembly.annotation.gtf)<br />
> -m (arg)               a temporary directory to store intermediate files<br />
> -n (arg)               the full path to the reference genome (default: ./Homo_sapiens_assembly38.fasta)<br />
> -o (arg)               the basename for the output file (default: GENE_PSI)<br />
> -p                     generate a plot of the PSI for the transcript with the highest exon count (requires ggplot2 R package)<br />
> -r (arg)               the length of the read. Use the shortest sample read length, longer reads are trimmed (default: 50)<br />
><br />
> -h                     print usage


Example Usage
-------------

|For the TTN gene with NO comparison:<br />
|  *get-psi-byIndi-github.sh -e G -g TTN -b bamlist.txt -o TTN_PSI*<br />

|For the MYH6-201 transcript comparing cases and controls and generating a plot:<br />
|  *get-psi-byIndi-github.sh -e T -g MYH6-201 -b bamlist.txt -cp -o MYH6-201_PSI*


Example BAM List
----------------

> sample1	/path/to/case1/bam/sample1.bam	case<br />
> sample2	/path/to/case2/bam/sample2.bam	case<br />
> sample3	/path/to/control1/bam/sample3.bam	control<br />
> sample4	/path/to/control2/bam/sample4.bam	control<br />


Output
------

|  * outname.fa<br />		
|          a pseudo-genome extending 10kb on either side of the gene of interest<br />
|    (can be loaded as genome into IGV)<br />
  * outname.gff<br />			
    a GFF file for all exons used<br />
  * outname.gff<br />			
    a GTF file for all exons used <br />
    (can be loaded into IGV with the outname.fa genome)<br />
  * *.bam<br />				
    reads aligned to outname.fa <br />
    (1 BAM/sample; can be loaded into IGV with the outname.fa genome)<br />
  * *junctions.bed<br />			
    a junctions file for each sample<br />
    (1 BED/sample; can be loaded into IGV with the outname.fa genome)<br />
  * *-case.psi<br />			
    a PSI file for each case sample <br />
    (columns: exon ID, exon length, N included reads, N excluded reads, and PSI)<br /> 
  * *-control.psi<br />			
    a PSI file for each control sample, requires -c flag<br />
    (columns: exon ID, exon length, N included reads, N excluded reads, and PSI)<br /> 
  * gene.bed<br />			
    a bed file with the coordinates used to generate outname.fa<br />
  * case-tot.psi<br />			
    a file containing the PSI for all cases along with the average<br />
  * control-tot.psi<br />		
    a file containing the PSI for all controls along with the average, requires -c flag<br /> 
  * psi-comparison.txt<br />		
    a file comparing cases and controls, requires -c flag<br />
    (columns: exon ID, average case PSI, average control PSI, exon, transcript, uncorrected t-test p-value)<br />
  * *pdf<br />				
    a plot of the PSI, requires -p flag<br />
  * outname-ref_coords.txt<br /> 	
    the exon coordinates relative to the reference genome (as opposed to outname.fa)<br />



