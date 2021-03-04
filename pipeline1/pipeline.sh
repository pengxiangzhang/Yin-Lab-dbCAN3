### let's begin with raw reads

## step1: trim raw reads

#input: sample_R1.fq.gz sample_R2.fq.gz
#output: sample_R1_val_1.fq.gz sample_R2_val_2.fq.gz
#here we may consider single end reads and interleaved reads in the feature
trim_galore --illumina -j 12 --paired sample_R1.fq.gz sample_R2.fq.gz

## step 2: align raw reads to genome, assembled by megahit or goldstard genome
## apply minimap2 to align raw reads to genome
## input genome: gold_standard_genome.fa, pair end raw reads: sample_R1_val_1.fq.gz sample_R2_val_2.fq.gz.

minimap2 -t 32 -x sr gold_standard_genome.fa sample_R1_val_1.fq.gz sample_R2_val_2.fq.gz >reads.paf

## step 3: convert fastq format to fasta format

seqtk seq -a sample_R1_val_1.fq.gz >sample_R1.fa
seqtk seq -a sample_R2_val_2.fq.gz >sample_R2.fa

## step 4: apply diamond to align raw reads(paired end) to CAZyme
diamond blastx --strand both --evalue 1e-6 --query sample_R1.fa --db ../dbcan --threads 32 --out sample_R1.fa.dmnd.out --outfmt 6
diamond blastx --strand both --evalue 1e-6 --query sample_R2.fa --db ../dbcan --threads 32 --out sample_R2.fa.dmnd.out --outfmt 6

## alternative step 4.1: apply bowtie2 or minimap to align the raw reads to CAZyme
bowtie2 -p 32 -x ../CAZyDB.version2.seperated sample_R1_val_1.fq.gz -S bowtie.R1.sam
bowtie2 -p 32 -x ../CAZyDB.version2.seperated sample_R2_val_2.fq.gz -S bowtie.R2.sam
# step 4.2  convert sam format to paf, bioconvert may be need to installed
bioconvert sam2paf bowtie.R1.sam bowtie.R1.paf
bioconvert sam2paf bowtie.R2.sam bowtie.R2.paf
cat bowtie.R1.paf bowtie.R2.paf >bowtie.paf

# also, minimap2 can be used to align the raw reads

## step 5:
## input: diamond align result:sample_R1.fa.dmnd.out sample_R2.fa.dmnd.out, reads align result: reads.paf
## output: best.hit.paf
python3 metagenome.py getpafdiamond --diamondR1 sample_R1.fa.dmnd.out --diamondR2 sample_R2.fa.dmnd.out --paf reads.paf >best.hit.paf

## step 6:
python3 metagenome.py compare_megahit_gold --gold_paf gold_standard.GH.fasta.paf --megahit_paf best.hit.paf

## step 7: ## the most import, need to been done, is a part of code in step 6
## in this step, we have not golden standard data,
python3 metagenome.py compare_megahit_gold --megahit_paf best.hit.paf

## get the diamond result
python3 metagenome.py analyzediomand --diamond sample_R1.fa.dmnd.out
## get cds sequence
## input: gold_standard_genome.fa overview.txt uniInput
## output a.fasta
python3 metagenome.py getcds --genomefasta gold_standard_genome.fa --overview overview.txt --uniInput uniInput
