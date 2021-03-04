#0 unpack
trim_galore --illumina -j 12 --paired sample_R1.fq.gz sample_R2.fq.gz
seqtk seq -a sample_R1_val_1.fq.gz >sample_R1.fa
seqtk seq -a sample_R2_val_2.fq.gz >sample_R2.fa

#1 bowtie2
bowtie2 -p 32 -x /mnt/array2/jinfang/dbcan_da/grant_test/CAZyDB.version2.seperated sample_R1_val_1.fq.gz -S bowtie.R1.sam
bowtie2 -p 32 -x /mnt/array2/jinfang/dbcan_da/grant_test/CAZyDB.version2.seperated sample_R2_val_2.fq.gz -S bowtie.R2.sam
# TODO: Change CAZyDB location.
bioconvert sam2paf bowtie.R1.sam output/bowtie.R1.paf
bioconvert sam2paf bowtie.R2.sam output/bowtie.R2.paf
python3 paf_result_analysis.py output/bowtie.R1.paf bowtie.R2.paf output/result/bowtie.sequence_FPKM.json result/bowtie.cazyfamily_FPKM.json

#2 diamond
diamond blastx --strand both --evalue 1e-6 --query sample_R1.fa --db ../dbcan --threads 32 --out output/diamond.R1.paf --outfmt 103
python3 diamond_paf.py output/diamond.R1.paf output/diamond.new.R1.paf
diamond blastx --strand both --evalue 1e-6 --query sample_R2.fa --db ../dbcan --threads 32 --out output/diamond.R2.paf --outfmt 103
python3 diamond_paf.py output/diamond.R2.paf output/diamond.new.R2.paf
python3 paf_result_analysis.py output/diamond.new.R1.paf output/diamond.new.R2.paf result/diamond.sequence_FPKM.json result/diamond.cazyfamily_FPKM.json

#3 minimap2
minimap2 -t 32 -x sr gold_standard_genome.fa sample_R1_val_1.fq.gz sample_R2_val_2.fq.gz >output/minimap2.paf

rm output/*
rm sample_R1_val_1.fq.gz sample_R2_val_2.fq.gz sample_R1.fa sample_R2.fa bowtie.R1.sam bowtie.R2.sam
