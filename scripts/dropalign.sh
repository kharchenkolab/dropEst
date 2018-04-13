#!/usr/bin/env bash
if [ "$#" -ne 3 ]; then
    echo "usage: dropalign.sh $tagged.fastq.gz genome.dir out.dir"
    echo "example: dropalign.sh SCG_27.tag.1.fastq.gz /groups/pklab/genomes/mm10 mouse.align"
    exit 1
fi
infile=$1
genome=$2
outdir=$3
nthreads=2
topdir=${outdir}/${infile%.fastq.gz}.tophat

#module load dev/java/jdk1.8
#module load seq/tophat/2.0.9
#export PATH=/opt/java/jdk1.8.0_45/bin:$PATH
echo $topdir
tophat2 -p ${nthreads} --no-coverage-search -g 1 -G ${genome}/genes.gtf -o ${topdir} ${genome}/Bowtie2Index/genome $infile
cat $topdir/align_summary.txt
outf=${outdir}/${infile%.fastq.gz}.aligned.bam
~/drop/Drop-seq_tools-1.0/TagReadWithGeneExon O=${outf} ANNOTATIONS_FILE=${genome}/genes.refFlat TAG=GE CREATE_INDEX=true INPUT=${topdir}/accepted_hits.bam
rm -rf $topdir
