#!/bin/bash

echo $0 $@

### USAGE

print_usage() {
  printf "

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 CALCULATE PERCENT SPLICED IN (PSI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 This script compares the PSI between case and control samples
  for RNA-Seq data.

 PSI is calculated on an individual basis and then averaged,
  NOT by comparing all cases merged vs. all controls merged.

 If you just want to calculate PSI (no case/control comparison),
  label all samples as \"case\"

 -------------------
 Required Programs
 -------------------

 This script uses the following programs:

  bedtools/2.23.0 *This version is needed to run coverageBed correctly*
  gcc/6.2.0
  bowtie2/2.2.9
  tophat/2.1.1
  R/3.6.1
  samtools/1.3.1



 This script also uses the following Rscripts:

  configure-psi-byIndiv.R
  comp-PSI.R
  plot-psi.R
  plot-psi-noComp.R

 -------------------
 Flags
 -------------------

 -b (arg, required)	a file with col1=sample_id, col2=bam, and col2=case/control (example BAM list below)
 -e (arg, required)	G/T - \"G\" to calculate for a gene or \"T\" to calculate for a transcript
 -g (arg, required)	gene/transcript name

 -c			compare cases and controls (must be specified in bamlist)
 -f (arg)		the full path to the GTF file (default: ./hg38.gencode.v27.primary_assembly.annotation.gtf)
 -m (arg)		a temporary directory to store intermediate files
 -n (arg)		the full path to the reference genome (default: ./Homo_sapiens_assembly38.fasta)
 -o (arg)		the basename for the output file (default: GENE_PSI)
 -p 			generate a plot of the PSI for the transcript with the highest exon count (requires ggplot2 R package)
 -r (arg)		the length of the read. Use the shortest sample read length, longer reads are trimmed (default: 50)
 -u			unpaired reads
 -x			remove bam files

 -h			print usage

 -------------------
 Example Usage
 -------------------

 get-psi-byIndi.sh -e G -g TTN -b bamlist.txt -o TTN_PSI

 -------------------
 Example BAM List
 -------------------

 sample1	/path/to/case1/bam/sample1.bam	case
 sample2	/path/to/case2/bam/sample2.bam	case
 sample3	/path/to/control1/bam/sample3.bam	control
 sample4	/path/to/control2/bam/sample4.bam	control

 -------------------
 Output
 -------------------

 outname.fa			a pseudo-genome extending 10kb on either side of the gene of interest
				 (can be loaded as genome into IGV)
 outname.gff			a GFF file for all exons used
 outname.gff			a GTF file for all exons used
				 (can be loaded into IGV with the outname.fa genome)
 *.bam				reads aligned to outname.fa
				 (1 BAM/sample; can be loaded into IGV with the outname.fa genome)
 *junctions.bed			a junctions file for each sample
				 (1 BED/sample; can be loaded into IGV with the outname.fa genome)
 *-case.psi			a PSI file for each case sample
				 (columns: exon ID, exon length, N included reads, N excluded reads, and PSI)
 *-control.psi			a PSI file for each control sample, requires -c flag
				 (columns: exon ID, exon length, N included reads, N excluded reads, and PSI)
 gene.bed			a bed file with the coordinates used to generate outname.fa
 case-tot.psi			a file containing the PSI for all cases along with the average
 control-tot.psi		a file containing the PSI for all controls along with the average, requires -c flag
 psi-comparison.txt		a file comparing cases and controls, requires -c flag
				 (columns: exon ID, average case PSI, average control PSI, exon, transcript, uncorrected t-test p-value)
 *pdf				a plot of the PSI, requires -p flag
 outname-ref_coords.txt 	the exon coordinates relative to the reference genome (as opposed to outname.fa)


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



"
}

### Get Arguments

bamlist=''
got=''
gene=''

comp=''
plot=''
gtf="./gencode.v27.primary_assembly.annotation.gtf"
tmp="./"
genome=./Homo_sapiens_assembly38.fasta
outname=''
rl=50
unp=''
rm_bam=''

while getopts 'b:e:g:cf:m:n:o:pr:uxh' flag; do
  case "${flag}" in
    b) bamlist="${OPTARG}" ;;
    e) got="${OPTARG}" ;;
    g) gene="${OPTARG}" ;;
    c) comp=true ;;
    f) gtf="${OPTARG}" ;;
    m) tmp="${OPTARG}" ;;
    n) genome="${OPTARG}" ;;
    o) outname="${OPTARG}" ;;
    p) plot=true ;;
    r) rl="${OPTARG}" ;;
    u) unp=true ;;
    x) rm_bam=true ;;
    h) print_usage ; exit ;;
    *) print_usage
       exit 1 ;;
  esac
done

if [ "${got,,}" != "g" ] && [ "${got,,}" != "gene" ] && [ "${got,,}" != "t" ] && [ "${got,,}" != "transcript" ]; then
  print_usage
  echo "ERROR: -e argument missing or invalid, please enter a \"G\" or a \"T\""
  exit 1
elif [ -z $gene ]; then
  print_usage
  echo "ERROR: -g argument missing, please enter a gene or transcript ID"
  exit 1
elif [ -z $bamlist ]; then
  print_usage
  echo "ERROR: -b argument missing, please enter a list of bam files to analyze"
  exit 1
elif [ ! -s $bamlist ]; then
  print_usage
  echo "ERROR: list of bam files not found"
  exit 1
elif [ -z $outname ]; then
  echo "output can be found in ${gene}_PSI-out"
  outname=${gene}_PSI
fi

### Determine the comp variable

if [ -z $comp ] && [ $(grep "control"$ $bamlist | wc -l) -gt 0 ] ; then

  printf "

Some samples are listed as controls.
No comparison will be made because the -c flag was not included.
If you would like to make a comparison, please run again using the -c flag.

"

elif [ -n "$comp" ] && [ $(grep "control"$ $bamlist | wc -l) -eq 0 ] ; then

  printf "

You requested a comparison (-c), but none of your files are listed as controls.
Please either run again excluding the -c flag or specify controls in your bam list.

"

fi




### Generate run-specific files

dirout="${outname}-out"
mkdir $dirout

bed=${dirout}/${gene}.bed
fasta=${dirout}/${outname}
ogtf=${dirout}/${outname}.gtf


## Make GTF file

grep ^"#" $gtf > ${ogtf}
if [ ${got,,} = "g" ] || [ ${got,,} = "gene" ]
then
  chr=$(grep "gene_name \"${gene}\"" $gtf | cut -f1 | uniq)
  start=$(( $(grep "gene_name \"${gene}\"" $gtf | cut -f4 | sort -g | head -n1) - 10000 ))
  end=$(( $(grep "gene_name \"${gene}\"" $gtf | cut -f5 | sort -g | tail -n1) + 10000 ))
  echo -e "$chr\t$start\t$end" > $bed

  if [ -n "$(echo $chr | grep chr)" ] ; then
    echo -e "$(echo $chr | cut -c4- )\t$start\t$end" >> $bed
    bedtools intersect -a $gtf -b <(grep ^"chr" ${bed}) | awk -v start="$start" -F"\t" '{OFS=FS} {$4=$4-start;$5=$5-start; print}' | sed "s/$chr/$gene/g" | grep "gene_name \"${gene}\"" >> ${ogtf}
  else
    echo -e "chr$chr\t$start\t$end" >> $bed
      bedtools intersect -a $gtf -b <(grep -v ^"chr" ${bed}) | awk -v start="$start" -F"\t" '{OFS=FS} {$4=$4-start;$5=$5-start; print}' | sed "s/$chr/$gene/g" | grep "gene_name \"${gene}\"" >> ${ogtf}
  fi

elif [ ${got,,} = "t" ] || [ ${got,,} = "transcript" ]
then
  chr=$(grep "transcript_name \"${gene}\"" $gtf | cut -f1 | uniq)
  start=$(( $(grep "transcript_name \"${gene}\"" $gtf | cut -f4 | sort -g | head -n1) - 10000 ))
  end=$(( $(grep "transcript_name \"${gene}\"" $gtf | cut -f5 | sort -g | tail -n1) + 10000 ))
  echo -e "$chr\t$start\t$end" > $bed

  if [ -n "$(echo $chr | grep chr)" ] ; then
    echo -e "$(echo $chr | cut -c4- )\t$start\t$end" >> $bed
    bedtools intersect -a $gtf -b <(grep ^"chr" ${bed}) | awk -v start="$start" -F"\t" '{OFS=FS} {$4=$4-start;$5=$5-start; print}' | sed "s/$chr/$gene/g" | grep "transcript_name \"${gene}\"" >> ${ogtf}
  else
    echo -e "chr$chr\t$start\t$end" >> $bed
    bedtools intersect -a $gtf -b <(grep -v ^"chr" ${bed}) | awk -v start="$start" -F"\t" '{OFS=FS} {$4=$4-start;$5=$5-start; print}' | sed "s/$chr/$gene/g" | grep "transcript_name \"${gene}\"" >> ${ogtf}
  fi

fi


## Make GFF file

gff="${dirout}/${outname}.gff"

#grep -v ^"#" $ogtf | cut -f 1-7 -d ";" | cut -f 3-9 | grep ^"exon" | sed 's/ /\t/g' | sed 's/[";]//g' | awk -v gene=$gene 'BEGIN{OFS="\t"}{print gene,".","exon_part",$2,$3,".",".",".",$20";"$18}' > $gff

# Exons

grep -v ^"#" $ogtf | cut -f 1-7 -d ";" | cut -f 3-9 | grep ^"exon" | sed 's/ /\t/g' | sed 's/[";]//g' | awk -v gene=$gene 'BEGIN{OFS="\t"}{print gene,".","exon_part",$2,$3,".",".","."}' | paste - <(cut -f 3- $ogtf | grep ^"exon" | grep -o "exon_number .*" | cut -f1 -d\; | grep -o [1234567890]\* | paste -d\; - <(cut -f 3- $ogtf | grep ^"exon" | grep -o "transcript_id.*" | cut -f1 -d\; | awk '{print $2}' | sed 's/"//g')) > $gff

# Introns

grep -v ^"#" $ogtf | cut -f 1-7 -d ";" | cut -f 3-9 | grep ^"intron" | sed 's/ /\t/g' | sed 's/[";]//g' | awk -v gene=$gene 'BEGIN{OFS="\t"}{print gene,".","intron_part",$2,$3,".",".","."}' | paste - <(cut -f 3- $ogtf | grep ^"intron" | grep -o "intron_number .*" | cut -f1 -d\; | grep -o [1234567890]\* | sed "s/$/;intron/g" | paste -d\; - <(cut -f 3- $ogtf | grep ^"intron" | grep -o "transcript_id.*" | cut -f1 -d\; | awk '{print $2}' | sed 's/"//g')) >> $gff


## Make translation file for reference genome

trans="${dirout}/${outname}-ref_coords.txt"
awk -v chr=$chr -v start=$start 'BEGIN{OFS="\t"}{print $1,$9,chr,$4+start,$5+start}' $gff > $trans


## Make fasta file for the provided coordinates

bedtools getfasta -fi $genome -bed $bed -fo ${fasta}.fa
sed "1c\>${gene}" ${fasta}.fa | sed "2 s/.\{50\}/&\n/g" > ${fasta}-temp.fa
mv -f ${fasta}-temp.fa ${fasta}.fa
bowtie2-build --large-index -f ${fasta}.fa ${fasta}





### Get fastq files for all cases and realign with TopHat

dir=${tmp}/${outname}
mkdir $dir

grep "case"$ $bamlist | cut -f1-2 | while read name i; do

  sam=${dir}/${name}.sam
  samtools view -h -L $bed $i > $sam
  bam=${dir}/${name}.bam
  samtools sort -n -o $bam $sam

  f1=${dir}/${name}.r1.fq
  f2=${dir}/${name}.r2.fq

  if [ -z "$unp" ]; then
    bedtools bamtofastq -i $bam -fq $f1 -fq2 $f2
  else
    bedtools bamtofastq -i $bam -fq $f1
  fi


  rm $sam $bam

  ## Trim fastqs to be correct read length and exclude read pairs shorter than provided read length

  grep -x "+" -B1 -A1 --no-group-separator $f1 | grep -v -x "+" | cut -c 1-${rl} | sed "s/^\(.\{$rl\}\)/&~~/g" | paste -d "~" <(grep -x "+" -B2 --no-group-separator $f1 | grep ^"[+@]") - > ${dir}/temp.fq; mv -f ${dir}/temp.fq $f1
  if [ -z "$unp" ]; then
    grep -x "+" -B1 -A1 --no-group-separator $f2 | grep -v -x "+" | cut -c 1-${rl} | sed "s/^\(.\{$rl\}\)/&||/g" | paste -d "~" <(grep -x "+" -B2 --no-group-separator $f2 | grep ^"[+@]") - > ${dir}/temp.fq; mv -f ${dir}/temp.fq $f2
    paste $f1 $f2 | grep -P "~~\t" | grep "||"$ > ${dir}/temp.fq
    cut -f1 ${dir}/temp.fq | sed 's/~~//g' | sed 's/~/\n/g' > $f1
    cut -f2 ${dir}/temp.fq | sed 's/||//g' | sed 's/~/\n/g' > $f2
    rm ${dir}/temp.fq
  else
    grep "~~"$ $f1 | sed 's/~~//g' | sed 's/~/\n/g' > ${dir}/temp.fq
    mv -f ${dir}/temp.fq $f1
  fi


  ## Use a different command depending on if there are paired or unpaired reads
  if [ -s $f1 ] && [ -s $f2 ]; then
    tophat -g 1 -N 1 -G $ogtf --output-dir ${dir}/case-tophat ${fasta} $f1 $f2
  elif [ -s $f1 ]; then
    tophat -g 1 -N 1 -G $ogtf --output-dir ${dir}/case-tophat ${fasta} $f1
  fi

  ## Move important files out of temporary directory
  mv ${dir}/case-tophat/accepted_hits.bam ${dirout}/${name}-case.bam
  mv ${dir}/case-tophat/junctions.bed ${dirout}/${name}-case-junctions.bed

  samtools index ${dirout}/${name}-case.bam

done

if [ -n "$comp" ]; then

  echo "Starting alignment of controls"

  # Get fastq files for all controls and realign with TopHat

  grep "control"$ $bamlist | cut -f1-2 | while read name i; do

    sam=${dir}/${name}.sam
    samtools view -h -L $bed $i > $sam
    bam=${dir}/${name}.bam
    samtools sort -n -o $bam $sam

    f1=${dir}/${name}.r1.fq
    f2=${dir}/${name}.r2.fq

    if [ -z "$unp" ]; then
      bedtools bamtofastq -i $bam -fq $f1 -fq2 $f2
    else
      bedtools bamtofastq -i $bam -fq $f1
    fi

    rm $sam $bam

    ## Trim fastqs to be correct read length and exclude read pairs shorter than provided read length

    grep -x "+" -B1 -A1 --no-group-separator $f1 | grep -v -x "+" | cut -c 1-${rl} | sed "s/^\(.\{$rl\}\)/&~~/g" | paste -d "~" <(grep -x "+" -B2 --no-group-separator $f1 | grep ^"[+@]") - > ${dir}/temp.fq; mv -f ${dir}/temp.fq $f1
    if [ -z "$unp" ]; then
      grep -x "+" -B1 -A1 --no-group-separator $f2 | grep -v -x "+" | cut -c 1-${rl} | sed "s/^\(.\{$rl\}\)/&||/g" | paste -d "~" <(grep -x "+" -B2 --no-group-separator $f2 | grep ^"[+@]") - > ${dir}/temp.fq; mv -f ${dir}/temp.fq $f2
      paste $f1 $f2 | grep -P "~~\t" | grep "||"$ > ${dir}/temp.fq
      cut -f1 ${dir}/temp.fq | sed 's/~~//g' | sed 's/~/\n/g' > $f1
      cut -f2 ${dir}/temp.fq | sed 's/||//g' | sed 's/~/\n/g' > $f2
      rm ${dir}/temp.fq
    else
      grep "~~"$ $f1 | sed 's/~~//g' | sed 's/~/\n/g' > ${dir}/temp.fq
      mv -f ${dir}/temp.fq $f1
    fi

  ## Use a different command depending on if there are paired or unpaired reads
  if [ -s $f1 ] && [ -s $f2 ]; then
    tophat -g 1 -N 1 -G $ogtf --output-dir ${dir}/control-tophat ${fasta} $f1 $f2
  elif [ -s $f1 ]; then
    tophat -g 1 -N 1 -G $ogtf --output-dir ${dir}/control-tophat ${fasta} $f1
  fi

  # Move important files out of scratch

  mv ${dir}/control-tophat/accepted_hits.bam ${dirout}/${name}-control.bam
  mv ${dir}/control-tophat/junctions.bed ${dirout}/${name}-control-junctions.bed

  samtools index ${dirout}/${name}-control.bam

  done

fi

rm -r $dir


### Calculate PSI for cases

for i in ${dirout}/*-case.bam; do

  casebam=${i}
  casebed=$(echo $i | sed 's/\.bam/-junctions.bed/g')

  coverageBed -split -abam ${casebam} -b $gff | awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$5-$4+1,$9,$10}' | sort -g -k 5 > ${dirout}/case_exonic_parts.inclusion

  sed 's/,/\t/g' ${casebed} | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+$13,$4,$5,$6}' > ${dirout}/left.bed
  sed 's/,/\t/g' ${casebed} | awk 'BEGIN{OFS="\t"}{print $1,$3-$14,$3,$4,$5,$6}' > ${dirout}/right.bed

  intersectBed -u -a ${dirout}/left.bed -b $gff > ${dirout}/left.overlap
  intersectBed -u -a ${dirout}/right.bed -b $gff > ${dirout}/right.overlap

  cat ${dirout}/left.overlap ${dirout}/right.overlap | cut -f4 | sort | uniq -c | awk '{ if($1 == 2) print $2 }' > ${dirout}/case_filtered_junctions.txt

  grep -F -f ${dirout}/case_filtered_junctions.txt ${casebed} > ${dirout}/filtered_junctions.bed

  sed 's/,/\t/g' ${dirout}/filtered_junctions.bed | grep -v description | awk '{OFS="\t"}{print $1,$2+$13,$3-$14,$4,$5,$6}' > ${dirout}/intron.bed

  intersectBed -wao -f 1.0 -a $gff -b ${dirout}/intron.bed | awk 'BEGIN{OFS="\t"}{$16 == 0? s[$9] += 0:s[$9] += $14}END{for (i in s) {print i,s[i]}}' | sort -g -k 1 > ${dirout}/case_exonic_parts.exclusion

  awk 'BEGIN{OFS="\t"}{print $5,$0}' ${dirout}/case_exonic_parts.inclusion | sort -g | cut -f 2- | paste - ${dirout}/case_exonic_parts.exclusion | awk -v "len=$rl" 'BEGIN{OFS="\t"; print "exon_ID" , "length" , "inclusion" , "exclusion" , "PSI"}{NIR=$6/($4+len-1) ; NER=$8/(len-1)}{print $5,$4,$6,$8,(NIR+NER<=0)? "NA":NIR / (NIR + NER)}' > ${dirout}/case_exonic_parts.psi

  sed 1d ${dirout}/case_exonic_parts.psi | sed 's/^e//g' | sort -g | grep -v intron| awk '{print "e"$0}' > ${dirout}/$(basename $i .bam).psi
  sed 1d ${dirout}/case_exonic_parts.psi | sed 's/^e//g' | sort -g | grep intron | awk '{print "i"$0}' | sed "s/;intron//g" >> ${dirout}/$(basename $i .bam).psi

  header="1iexon_ID\tlength\tinclusion\texclusion\t$(basename $i .bam)-psi"
  sed -i $header ${dirout}/$(basename $i .bam).psi

  rm ${dirout}/*exclusion ${dirout}/*inclusion ${dirout}/*exonic_parts.psi ${dirout}/*_filtered_junctions.txt ${dirout}/left* ${dirout}/right* ${dirout}/filtered_junctions.bed ${dirout}/intron.bed

  if [ -n "$rm_bam" ] ; then
    rm $casebam*
  fi

done

cut -f1 $(ls ${dirout}/*-case.psi | head -n1) > ${dirout}/case-tot.psi
for i in ${dirout}/*-case.psi ; do
  cut -f5 $i | paste ${dirout}/case-tot.psi - > ${dirout}/temp ; mv -f ${dirout}/temp ${dirout}/case-tot.psi
done

Rscript configure-psi-byIndiv.R ${dirout}/case-tot.psi



if [ -n "$comp" ] ; then

  # Calculate PSI for controls

  for i in ${dirout}/*-control.bam; do

    conbam=${i}
    conbed=$(echo $i | sed 's/\.bam/-junctions.bed/g')


    coverageBed -split -abam ${conbam} -b $gff | awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$5-$4+1,$9,$10}' | sort -g -k 5 > ${dirout}/control_exonic_parts.inclusion

    sed 's/,/\t/g' ${conbed} | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+$13,$4,$5,$6}' > ${dirout}/left.bed
    sed 's/,/\t/g' ${conbed} | awk 'BEGIN{OFS="\t"}{print $1,$3-$14,$3,$4,$5,$6}' > ${dirout}/right.bed

    intersectBed -u -a ${dirout}/left.bed -b $gff > ${dirout}/left.overlap
    intersectBed -u -a ${dirout}/right.bed -b $gff > ${dirout}/right.overlap

    cat ${dirout}/left.overlap ${dirout}/right.overlap | cut -f4 | sort | uniq -c | awk '{ if($1 == 2) print $2 }' > ${dirout}/control_filtered_junctions.txt

    grep -F -f ${dirout}/control_filtered_junctions.txt ${conbed} > ${dirout}/filtered_junctions.bed

    sed 's/,/\t/g' ${dirout}/filtered_junctions.bed | grep -v description | awk '{OFS="\t"}{print $1,$2+$13,$3-$14,$4,$5,$6}' > ${dirout}/intron.bed

    intersectBed -wao -f 1.0 -a $gff -b ${dirout}/intron.bed | awk 'BEGIN{OFS="\t"}{$16 == 0? s[$9] += 0:s[$9] += $14}END{for (i in s) {print i,s[i]}}' | sort -g -k 1 > ${dirout}/control_exonic_parts.exclusion

    awk 'BEGIN{OFS="\t"}{print $5,$0}' ${dirout}/control_exonic_parts.inclusion | sort -g | cut -f 2- | paste - ${dirout}/control_exonic_parts.exclusion | awk -v "len=$rl" 'BEGIN{OFS="\t"; print "exon_ID" , "length" , "inclusion" , "exclusion" , "PSI"}{NIR=$6/($4+len-1) ; NER=$8/(len-1)}{print $5,$4,$6,$8,(NIR+NER<=0)? "NA":NIR / (NIR + NER)}' > ${dirout}/control_exonic_parts.psi

    sed 1d ${dirout}/control_exonic_parts.psi | sed 's/^e//g' | sort -g | grep -v intron | awk '{print "e"$0}' > ${dirout}/$(basename $i .bam).psi
    sed 1d ${dirout}/control_exonic_parts.psi | sed 's/^e//g' | sort -g | grep intron | awk '{print "i"$0}' | sed "s/;intron//g" >> ${dirout}/$(basename $i .bam).psi

    header="1iexon_ID\tlength\tinclusion\texclusion\t$(basename $i .bam)-psi"
    sed -i $header ${dirout}/$(basename $i .bam).psi

    rm ${dirout}/*exclusion ${dirout}/*inclusion ${dirout}/*exonic_parts.psi ${dirout}/*_filtered_junctions.txt ${dirout}/left* ${dirout}/right* ${dirout}/filtered_junctions.bed ${dirout}/intron.bed

    if [ -n "$rm_bam" ] ; then
      rm $conbam*
    fi

  done


  cut -f1 $(ls ${dirout}/*-control.psi | head -n1) > ${dirout}/control-tot.psi
  for i in ${dirout}/*-control.psi ; do
    cut -f5 $i | paste ${dirout}/control-tot.psi - > ${dirout}/temp ; mv -f ${dirout}/temp ${dirout}/control-tot.psi
  done

  Rscript configure-psi-byIndiv.R ${dirout}/control-tot.psi

  Rscript comp-PSI.R ${dirout}/case-tot.psi ${dirout}/control-tot.psi ${dirout}/psi-comparison.txt

  if [ -n "$plot" ]; then
    Rscript plot-psi.R ${dirout}/psi-comparison.txt ${dirout}/${outname}-plot.pdf
  fi
elif [ -n "$plot" ]; then
    Rscript plot-psi-noComp.R ${dirout}/case-tot.psi ${dirout}/${outname}-plot.pdf
fi
