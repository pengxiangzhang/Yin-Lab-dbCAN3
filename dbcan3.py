import os
import argparse
import time
import sys

parser = argparse.ArgumentParser(description='dbCAN3')
parser.add_argument('-r1', type=str, help='R1 File', required=True)
parser.add_argument('-r2', type=str, help='R2 File', required=True)
args = parser.parse_args()

def main():
    start = time.time()

    # Prepare:
    sys.stdout.write("Running trim_galore")
    os.system("trim_galore --illumina --no_report_file --basename input -j 12 --paired "+args.r1+" "+args.r2+" -o output/")
    sys.stdout.write("Running seqtk R1")
    os.system("seqtk seq -a output/input_val_1.fq.gz  > output/R1.fa")
    sys.stdout.write("Running seqtk R2")
    os.system("seqtk seq -a output/input_val_2.fq.gz  > output/R2.fa")
    prepare = time.time()
    print('Prepare Run Time: ' + str(prepare - start))

    # Bowtie:
    sys.stdout.write("Running bowtie2 R1")
    os.system("bowtie2 -p 32 -x /home/penxiang/dbCAN3/test/bowtiebase/CAZyDB output/input_val_1.fq.gz -S output/bowtie.R1.sam")
    sys.stdout.write("Running bowtie2 R2")
    os.system("bowtie2 -p 32 -x /home/penxiang/dbCAN3/test/bowtiebase/CAZyDB output/input_val_2.fq.gz -S output/bowtie.R2.sam")
    sys.stdout.write("Running bioconvert R1")
    os.system("bioconvert sam2paf output/bowtie.R1.sam output/bowtie.R1.paf")
    sys.stdout.write("Running bioconvert R2")
    os.system("bioconvert sam2paf output/bowtie.R2.sam output/bowtie.R2.paf")
    sys.stdout.write("Running bowtie2 analysis")
    os.system("python3 paf_result_analysis.py output/bowtie.R1.paf output/bowtie.R2.paf result/bowtie2.sequence_FPKM.csv result/bowtie2.cazyfamily_FPKM.csv")
    bowtie = time.time()
    print('Bowtie Run Time: ' + str(bowtie - prepare))


    # diamond:
    sys.stdout.write("Running diamond R1")
    os.system("diamond blastx --strand both --evalue 1e-6 --query output/R1.fa --db diamondbase.dmnd --threads 32 --out output/diamond.R1.paf --outfmt 103")
    sys.stdout.write("Running diamond R2")
    os.system("diamond blastx --strand both --evalue 1e-6 --query output/R2.fa --db diamondbase.dmnd --threads 32 --out output/diamond.R2.paf --outfmt 103")
    sys.stdout.write("Running diamond_paf R1")
    os.system("python3 diamond_paf.py output/diamond.R1.paf output/diamond.new.R1.paf")
    sys.stdout.write("Running diamond_paf R2")
    os.system("python3 diamond_paf.py output/diamond.R2.paf output/diamond.new.R2.paf")
    sys.stdout.write("Running diamond analysis")
    os.system("python3 paf_result_analysis.py output/diamond.new.R1.paf output/diamond.new.R2.paf result/diamond.sequence_FPKM.csv result/diamond.cazyfamily_FPKM.csv")
    diamond = time.time()
    print('Bowtie Run Time: ' + str(diamond - bowtie))

    # minimap2
    sys.stdout.write("Running minimap2")
    os.system("minimap2 -t 32 -x sr CAZyDB.07312020.fa.cds output/output_val_1.fq.gz output/output_val_2.fq.gz > output/minimap2.paf")

    stop = time.time()
    print('Total Run Time: ' + str(stop - start))
if __name__ == '__main__':
    main()