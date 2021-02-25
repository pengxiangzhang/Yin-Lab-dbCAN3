import os
import argparse
import time

parser = argparse.ArgumentParser(description='dbCAN3')
parser.add_argument('-b', '--bowtie', dest='bowtie', action='store_true', help='Run Bowtie2', default=False)
parser.add_argument('-d', '--diamond', dest='diamond', action='store_true', help='Run diamond', default=False)
parser.add_argument('-m', '--minimap', dest='minimap', action='store_true', help='Run minimap2', default=False)
parser.add_argument('-r1', dest='input1', type=str, help='R1 File', required=True)
parser.add_argument('-r2', dest='input2', type=str, help='R2 File', required=True)
args = parser.parse_args()


def check_return(status):
    if status != 0:
        print("Error")
        return 1


def main():
    if not any([args.bowtie, args.diamond, args.minimap]):
        print("ERROR: You must select a program to run(--bowtie, --diamond or --minimap).")
        return 1

    if (args.input1.split(".")[-2]) != (args.input2.split(".")[-2]):
        print("You must have the same file type for both of your input.")
        return 1

    start = time.time()
    # Prepare:
    if (args.input1.split(".")[-2]) == "fq":
        print("Running trim_galore")
        result = os.system(
            "trim_galore --illumina --no_report_file --basename input -j 12 --paired " + args.input1 + " " + args.input2 + " -o output/")
        check_return(result)
        prepare_end = time.time()
        print('Prepare Run Time(second): ' + str(prepare_end - start))

    if args.bowtie:
        bowtie_start = time.time()
        print("Running bowtie2 R1")
        result = os.system(
            "bowtie2 -p 32 -x /home/penxiang/dbCAN3/test/bowtiebase/CAZyDB output/input_val_1.fq.gz -S output/bowtie.R1.sam")
        check_return(result)
        print("Running bowtie2 R2")
        result = os.system(
            "bowtie2 -p 32 -x /home/penxiang/dbCAN3/test/bowtiebase/CAZyDB output/input_val_2.fq.gz -S output/bowtie.R2.sam")
        check_return(result)
        print("Running bioconvert R1")
        result = os.system("bioconvert sam2paf output/bowtie.R1.sam output/bowtie.R1.paf")
        check_return(result)
        print("Running bioconvert R2")
        result = os.system("bioconvert sam2paf output/bowtie.R2.sam output/bowtie.R2.paf")
        check_return(result)
        print("Running bowtie2 analysis")
        result = os.system(
            "python3 paf_result_analysis.py output/bowtie.R1.paf output/bowtie.R2.paf result/bowtie2.sequence_FPKM.csv result/bowtie2.cazyfamily_FPKM.csv")
        check_return(result)
        bowtie_end = time.time()
        print('Bowtie Run Time(second): ' + str(bowtie_end - bowtie_start))

    if args.diamond:
        diamond_start = time.time()
        print("Running seqtk R1")
        result = os.system("seqtk seq -a output/input_val_1.fq.gz  > output/R1.fa")
        check_return(result)
        print("Running seqtk R2")
        result = os.system("seqtk seq -a output/input_val_2.fq.gz  > output/R2.fa")
        check_return(result)
        print("Running diamond R1")
        result = os.system(
            "diamond blastx --strand both --evalue 1e-6 --query output/R1.fa --db diamondbase.dmnd --threads 32 --out output/diamond.R1.paf --outfmt 103")
        check_return(result)
        print("Running diamond R2")
        result = os.system(
            "diamond blastx --strand both --evalue 1e-6 --query output/R2.fa --db diamondbase.dmnd --threads 32 --out output/diamond.R2.paf --outfmt 103")
        check_return(result)
        print("Running diamond_paf R1")
        result = os.system("python3 diamond_paf.py output/diamond.R1.paf output/diamond.new.R1.paf")
        check_return(result)
        print("Running diamond_paf R2")
        result = os.system("python3 diamond_paf.py output/diamond.R2.paf output/diamond.new.R2.paf")
        check_return(result)
        print("Running diamond analysis")
        result = os.system(
            "python3 paf_result_analysis.py output/diamond.new.R1.paf output/diamond.new.R2.paf result/diamond.sequence_FPKM.csv result/diamond.cazyfamily_FPKM.csv")
        check_return(result)
        diamond_end = time.time()
        print('Bowtie Run Time(second): ' + str(diamond_end - diamond_start))

    if args.minimap:
        minimap_start = time.time()
        print("Running minimap2")
        result = os.system(
            "minimap2 -t 32 -x sr CAZyDB.07312020.fa.cds output/input_val_1.fq.gz output/input_val_2.fq.gz > output/minimap2.paf")
        check_return(result)
        print("Running minimap_paf")
        result = os.system("python3 minimap_paf.py output/minimap2.paf output/minimap2_r1.paf output/minimap2_r2.paf")
        check_return(result)
        print("Running minimap analysis")
        result = os.system(
            "python3 paf_result_analysis.py output/minimap2_r1.paf output/minimap2_r2.paf result/minimap.sequence_FPKM.csv result/minimap.cazyfamily_FPKM.csv")
        check_return(result)
        minimap_end = time.time()
        print('Minimap Run Time(second): ' + str(minimap_end - minimap_start))

    # os.system("rm -rf output/*")
    stop = time.time()
    print('Total Run Time: ' + str(stop - start))


if __name__ == '__main__':
    main()
