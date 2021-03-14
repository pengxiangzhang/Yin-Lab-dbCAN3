import os
import argparse
import time

# TODO: NEW
## TODO: New pipeline: read to contig -> megahit/ Long sequence
## TODO: Long sequence -> dbCAN -> annotation(CAZYFamily) .gff -> stringtie (MPKM TPM)


parser = argparse.ArgumentParser(description='pipeline2')
parser.add_argument('-o', dest='output', type=str, help='Output location?', default="output/")
parser.add_argument('-r1', dest='input1', type=str, help='R1 File', required=True)
parser.add_argument('-r2', dest='input2', type=str, help='R2 File', required=True)
args = parser.parse_args()


def check_return(status):
    if status != 0:
        print("Error")
        exit(1)


def main():
    print("test")


# megahit -1 /mnt/array2/pengxiang/annotation/reads/fq/example0/sample_R1_val_1.fq.gz -2 /mnt/array2/pengxiang/annotation/reads/fq/example0/sample_R2_val_2.fq.gz -o megahit
# run_dbcan.py /mnt/array2/pengxiang/annotation/pipeline2/megahit/final.contigs.fa prok --out_dir dbcan_out --db_dir /mnt/array2/pengxiang/annotation/pipeline2/dbcan/db --dia_cpu 4  --hmm_cpu=4 --hotpep_cpu 4 --tf_cpu 4 --tf_cpu 4
# TODO: Change database
# python3 filter_gff.py /mnt/array2/pengxiang/annotation/pipeline2/dbcan_out/prodigal.gff /mnt/array2/pengxiang/annotation/pipeline2/dbcan_out/overview.txt /mnt/array2/pengxiang/annotation/pipeline2/new.gff
# samtools sort -@ 32 -O bam -o R1.bam R1.sam
# samtools sort -@ 32 -O bam -o R2.bam R2.sam
# stringtie -e -B -p 8 -G /mnt/array2/pengxiang/annotation/pipeline2/new.gff  -o R1.gtf /mnt/array2/pengxiang/annotation/pipeline2/bowtie/R1.bam /mnt/array2/pengxiang/annotation/pipeline2/bowtie/R2.bam


if __name__ == '__main__':
    main()
