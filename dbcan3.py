import os
import argparse
import time
import sys

parser = argparse.ArgumentParser(description='dbCAN3')
parser.add_argument('-b', '--bowtie', dest='bowtie', action='store_true', help='Run Bowtie2', default=False)
parser.add_argument('-d', '--diamond', dest='diamond', action='store_true', help='Run diamond', default=False)
parser.add_argument('-m', '--minimap', dest='minimap', action='store_true', help='Run minimap2', default=False)
parser.add_argument('-r1', dest='input1', type=str, help='R1 File', required=True)
parser.add_argument('-r2', dest='input2', type=str, help='R2 File', required=True)
args = parser.parse_args()


def main():
    start = time.time()

    # Prepare:
    print("Running trim_galore")
    os.system(
        "trim_galore --illumina --no_report_file --basename input -j 12 --paired " + args.input1 + " " + args.input2 + " -o output/")
    prepare_end = time.time()
    print('Prepare Run Time(second): ' + str(prepare_end - start))

    if args.bowtie:
        bowtie_start = time.time()
        print("Running bowtie2 R1")
        os.system(
            "bowtie2 -p 32 -x /home/penxiang/dbCAN3/test/bowtiebase/CAZyDB output/input_val_1.fq.gz -S output/bowtie.R1.sam")
        print("Running bowtie2 R2")
        os.system(
            "bowtie2 -p 32 -x /home/penxiang/dbCAN3/test/bowtiebase/CAZyDB output/input_val_2.fq.gz -S output/bowtie.R2.sam")
        print("Running bioconvert R1")
        os.system("bioconvert sam2paf output/bowtie.R1.sam output/bowtie.R1.paf")
        print("Running bioconvert R2")
        os.system("bioconvert sam2paf output/bowtie.R2.sam output/bowtie.R2.paf")
        print("Running bowtie2 analysis")
        os.system(
            "python3 paf_result_analysis.py output/bowtie.R1.paf output/bowtie.R2.paf result/bowtie2.sequence_FPKM.csv result/bowtie2.cazyfamily_FPKM.csv")
        bowtie_end = time.time()
        print('Bowtie Run Time(second): ' + str(bowtie_end - bowtie_start))

    if args.diamond:
        diamond_start = time.time()
        print("Running seqtk R1")
        os.system("seqtk seq -a output/input_val_1.fq.gz  > output/R1.fa")
        print("Running seqtk R2")
        os.system("seqtk seq -a output/input_val_2.fq.gz  > output/R2.fa")
        print("Running diamond R1")
        os.system(
            "diamond blastx --strand both --evalue 1e-6 --query output/R1.fa --db diamondbase.dmnd --threads 32 --out output/diamond.R1.paf --outfmt 103")
        print("Running diamond R2")
        os.system(
            "diamond blastx --strand both --evalue 1e-6 --query output/R2.fa --db diamondbase.dmnd --threads 32 --out output/diamond.R2.paf --outfmt 103")
        print("Running diamond_paf R1")
        os.system("python3 diamond_paf.py output/diamond.R1.paf output/diamond.new.R1.paf")
        print("Running diamond_paf R2")
        os.system("python3 diamond_paf.py output/diamond.R2.paf output/diamond.new.R2.paf")
        print("Running diamond analysis")
        os.system(
            "python3 paf_result_analysis.py output/diamond.new.R1.paf output/diamond.new.R2.paf result/diamond.sequence_FPKM.csv result/diamond.cazyfamily_FPKM.csv")
        diamond_end = time.time()
        print('Bowtie Run Time(second): ' + str(diamond_end - diamond_start))

    if args.minimap:
        minimap_start = time.time()
        print("Running minimap2")
        os.system(
            "minimap2 -t 32 -x sr CAZyDB.07312020.fa.cds output/input_val_1.fq.gz output/input_val_2.fq.gz > output/minimap2.paf")
        print("Running minimap_paf")
        os.system("python3 minimap_paf.py output/minimap2.paf output/minimap2_r1.paf output/minimap2_r2.paf")
        print("Running minimap analysis")
        os.system(
            "python3 paf_result_analysis.py output/minimap2_r1.paf output/minimap2_r2.paf result/minimap.sequence_FPKM.csv result/minimap.cazyfamily_FPKM.csv")
        minimap_end = time.time()
        print('Minimap Run Time(second): ' + str(minimap_end - minimap_start))

    # os.system("rm -rf output/*")
    stop = time.time()
    print('Total Run Time: ' + str(stop - start))


if __name__ == '__main__':
    main()
