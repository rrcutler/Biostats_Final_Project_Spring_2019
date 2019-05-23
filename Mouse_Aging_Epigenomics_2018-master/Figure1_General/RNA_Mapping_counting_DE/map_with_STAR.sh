#!/bin/bash

# set the name of the job
#$ -N Heart_star
#
# set the maximum memory usage
#$ -l h_vmem=8G
#
#$ -q extended
#
#$ -pe shm 6
#
# set the maximum run time
#$ -l h_rt=150:00:00
#
#$ -l h_stack=15M
#
# send mail when job ends or aborts
#$ -m ea
#
# specify an email address
#$ -M benayoun@stanford.edu
#
# check for errors in the job submission options
#$ -w e
#
#$ -R y


cd /srv/gs1/projects/brunet/BB/Heart/;

export GENIN=/srv/gs1/projects/brunet/BB/MM9_ERCC
export START_EXEC=/srv/gs1/projects/brunet/tools/STAR-STAR_2.4.0j/source

$START_EXEC/STAR --genomeDir $GENIN --readFilesIn 3m4_heart_ATCACG_R1_val_1.fq 3m4_heart_ATCACG_R2_val_2.fq --runThreadN 6 --outFilterMultimapNmax 10 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix 3m4_heart_ATCACG
$START_EXEC/STAR --genomeDir $GENIN --readFilesIn 3m5_heart_CGATGT_R1_val_1.fq 3m5_heart_CGATGT_R2_val_2.fq --runThreadN 6 --outFilterMultimapNmax 10 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix 3m5_heart_CGATGT
$START_EXEC/STAR --genomeDir $GENIN --readFilesIn 3m6_heart_TTAGGC_R1_val_1.fq 3m6_heart_TTAGGC_R2_val_2.fq --runThreadN 6 --outFilterMultimapNmax 10 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix 3m6_heart_TTAGGC
$START_EXEC/STAR --genomeDir $GENIN --readFilesIn 12m4_heart_TGACCA_R1_val_1.fq 12m4_heart_TGACCA_R2_val_2.fq --runThreadN 6 --outFilterMultimapNmax 10 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix 12m4_heart_TGACCA
$START_EXEC/STAR --genomeDir $GENIN --readFilesIn 12m5_heart_ACAGTG_R1_val_1.fq 12m5_heart_ACAGTG_R2_val_2.fq --runThreadN 6 --outFilterMultimapNmax 10 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix 12m5_heart_ACAGTG
$START_EXEC/STAR --genomeDir $GENIN --readFilesIn 12m6_heart_GCCAAT_R1_val_1.fq 12m6_heart_GCCAAT_R2_val_2.fq --runThreadN 6 --outFilterMultimapNmax 10 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix 12m6_heart_GCCAAT
$START_EXEC/STAR --genomeDir $GENIN --readFilesIn 29m4_heart_CAGATC_R1_val_1.fq 29m4_heart_CAGATC_R2_val_2.fq --runThreadN 6 --outFilterMultimapNmax 10 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix 29m4_heart_CAGATC
$START_EXEC/STAR --genomeDir $GENIN --readFilesIn 29m5_heart_ACTTGA_R1_val_1.fq 29m5_heart_ACTTGA_R2_val_2.fq --runThreadN 6 --outFilterMultimapNmax 10 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix 29m5_heart_ACTTGA
$START_EXEC/STAR --genomeDir $GENIN --readFilesIn 29m6_heart_GATCAG_R1_val_1.fq 29m6_heart_GATCAG_R2_val_2.fq --runThreadN 6 --outFilterMultimapNmax 10 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFileNamePrefix 29m6_heart_GATCAG

module load samtools

for f in $(find . -name '*.sam')
do
	fileName=$(basename "${f}" | sed 's/\.sam/\.bam/g');
 	filePath="."
    oFname="${filePath}/${fileName}"
    samtools view -b -S $f > $oFname
done

