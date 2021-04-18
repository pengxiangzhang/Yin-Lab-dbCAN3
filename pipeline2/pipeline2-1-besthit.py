import os
import argparse
import time

parser = argparse.ArgumentParser(description='pipeline1')
parser.add_argument('-b', '--bowtie2', dest='bowtie', action='store_true', help='Run Bowtie2?', default=False)
parser.add_argument('-d', '--diamond', dest='diamond', action='store_true', help='Run diamond?', default=False)
parser.add_argument('-m', '--minimap2', dest='minimap', action='store_true', help='Run minimap2?', default=False)
parser.add_argument('-a', '--bwa', dest='bwa', action='store_true', help='Run BWA?', default=False)
parser.add_argument('-s', '--source', dest='source', action='store_true', help='Store source file?', default=False)
parser.add_argument('-o', dest='output', type=str, help='Output location?', default="output/")
parser.add_argument('-r1', dest='input1', type=str, help='R1 File', required=True)
parser.add_argument('-r2', dest='input2', type=str, help='R2 File', required=True)
args = parser.parse_args()


def check_return(status):
    if status != 0:
        print("Error")
        exit(1)


def get_count_reads_fq(file):
    r = os.popen("zcat " + file + " | echo $((`wc -l`/4))")
    text = r.read()
    r.close()
    return str(int(text))


def get_count_reads_fa(file):
    r = os.popen("grep '>' " + file + " | wc -l")
    text = r.read()
    r.close()
    return str(int(text))


def main():
    if not any([args.bowtie, args.diamond, args.minimap]):
        print("ERROR: You must select a program to run(--bowtie, --diamond or --minimap).")
        exit(1)

    if not args.output.endswith('/'):
        args.output = args.output + "/"
    print("Writing to: " + args.output)

    start = time.time()
    r2counts = "0"

    # Prepare:
    check_return(os.system("mkdir -p " + args.output))
    check_return(os.system("mkdir -p " + args.output + "data/"))
    check_return(os.system("mkdir -p " + args.output + "result/"))
    print("Start Read Count")
    if (args.input1.split(".")[-2]) == "fq":
        check_return(os.system(
            "cp " + args.input1 + " " + args.output + "data/input1.fq.gz && cp " + args.input2 + " " + args.output + "data/input2.fq.gz"))
        print("Running trim_galore")
        check_return(os.system(
            "trim_galore --illumina --no_report_file -j 12 --paired " + args.output + "data/input1.fq.gz " + args.output + "data/input2.fq.gz -o " + args.output + "data/"))
        check_return(os.system("rm -rf " + args.output + "data/input1.fq.gz " + args.output + "data/input2.fq.gz"))
        r1counts = get_count_reads_fq(args.input1)
        fq = True

    elif (args.input1.split(".")[-1]) == "fa":
        r1counts = get_count_reads_fa(args.input1)
        check_return(os.system(
            "cp " + args.input1 + " " + args.output + "data/R1.fa && cp " + args.input2 + " " + args.output + "data/R2.fa"))
        fq = False

    else:
        print("Error: File type not supported, please provide .fa or .fq.gz datatype.")
        exit(1)

    prepare_end = time.time()
    print("R1 Read Count: " + r1counts)
    print('Prepare Run Time(second): ' + str(prepare_end - start))

    if args.bowtie:
        bowtie_start = time.time()
        print("Running bowtie2")
        if fq:
            check_return(os.system(
                "bowtie2 -p 32 -x " + args.output + "bowtie2-index/cds " + args.output + "data/input1_val_1.fq.gz -S " + args.output + "data/bowtie.R1.sam"))
            check_return(os.system(
                "bowtie2 -p 32 -x " + args.output + "bowtie2-index/cds " + args.output + "data/input2_val_2.fq.gz -S " + args.output + "data/bowtie.R2.sam"))
        else:
            check_return(os.system(
                "bowtie2 -p 32 -x " + args.output + "bowtie2-index/cds -f " + args.output + "data/R1.fa -S " + args.output + "data/bowtie.R1.sam"))
            check_return(os.system(
                "bowtie2 -p 32 -x " + args.output + "bowtie2-index/cds -f " + args.output + "data/R2.fa -S " + args.output + "data/bowtie.R2.sam"))

        print("Running bioconvert")
        check_return(os.system(
            "bioconvert sam2paf " + args.output + "data/bowtie.R1.sam " + args.output + "data/bowtie.R1.paf"))
        check_return(os.system(
            "bioconvert sam2paf " + args.output + "data/bowtie.R2.sam " + args.output + "data/bowtie.R2.paf"))
        print("Running bowtie2 analysis")
        check_return(os.system(
            "python3 paf_result_analysis.py " + args.output + "data/bowtie.R1.paf " + args.output + "data/bowtie.R2.paf " + r1counts + " " + r2counts + " " + args.output + "result/bowtie2.sequence_FPKM.csv " + args.output + "result/bowtie2.cazyfamily_FPKM.csv"))
        bowtie_end = time.time()
        print('Bowtie Run Time(second): ' + str(bowtie_end - bowtie_start))

    if args.diamond:
        diamond_start = time.time()
        if fq:
            print("Running seqtk")
            check_return(os.system(
                "seqtk seq -a " + args.output + "data/input1_val_1.fq.gz  > " + args.output + "data/R1.fa"))
            check_return(os.system(
                "seqtk seq -a " + args.output + "data/input2_val_2.fq.gz  > " + args.output + "data/R2.fa"))
        print("Running diamond")
        check_return(os.system(
            "diamond blastx --strand both -k 1 --evalue 1e-10 --query " + args.output + "data/R1.fa --db " + args.output + "diamond-db/cds.dmnd --threads 32 --out " + args.output + "data/diamond.R1.paf --outfmt 103"))
        check_return(os.system(
            "diamond blastx --strand both -k 1 --evalue 1e-10 --query " + args.output + "data/R2.fa --db " + args.output + "diamond-db/cds.dmnd --threads 32 --out " + args.output + "data/diamond.R2.paf --outfmt 103"))

        print("Running Removing unhit")
        check_return(os.system(
            "python3 paf_remove_unhit.py " + args.output + "data/diamond.R1.paf " + args.output + "data/diamond.new.R1.paf"))
        check_return(os.system(
            "python3 paf_remove_unhit.py " + args.output + "data/diamond.R2.paf " + args.output + "data/diamond.new.R2.paf"))
        print("Running diamond analysis")
        check_return(os.system(
            "python3 paf_result_analysis.py " + args.output + "data/diamond.new.R1.paf " + args.output + "data/diamond.new.R2.paf " + r1counts + " " + r2counts + " " + args.output + "result/diamond.sequence_FPKM.csv " + args.output + "result/diamond.cazyfamily_FPKM.csv"))
        diamond_end = time.time()
        print('Bowtie Run Time(second): ' + str(diamond_end - diamond_start))

    if args.minimap:
        minimap_start = time.time()
        print("Running minimap2")
        if fq:
            check_return(os.system(
                "minimap2 -t 32 -x sr " + args.output + "minimap-db/cds.fa " + args.output + "data/input1_val_1.fq.gz  " + args.output + "data/input2_val_2.fq.gz >  " + args.output + "data/minimap2.paf"))
        else:
            check_return(os.system(
                "minimap2 -t 32 -x sr " + args.output + "minimap-db/cds.fa " + args.output + "data/R1.fa  " + args.output + "data/R2.fa >  " + args.output + "data/minimap2.paf"))
        print("Running split paf")
        check_return(os.system(
            "python3 paf_split.py " + args.output + "data/minimap2.paf " + args.output + "data/minimap2_r1.paf " + args.output + "data/minimap2_r2.paf"))
        print("Running minimap analysis")
        check_return(os.system(
            "python3 paf_result_analysis.py " + args.output + "data/minimap2_r1.paf " + args.output + "data/minimap2_r2.paf " + r1counts + " " + r2counts + " " + args.output + "result/minimap.sequence_FPKM.csv " + args.output + "result/minimap.cazyfamily_FPKM.csv"))
        minimap_end = time.time()
        print('Minimap Run Time(second): ' + str(minimap_end - minimap_start))

    if args.bwa:
        bwa_start = time.time()
        print("Running BWA")
        if fq:
            check_return(os.system(
                "bwa mem -t 32 " + args.output + "bwa-db/cds.fa " + args.output + "data/input1_val_1.fq.gz > " + args.output + "data/bwa_R1.sam"))
            check_return(os.system(
                "bwa mem -t 32 " + args.output + "bwa-db/cds.fa " + args.output + "data/input2_val_2.fq.gz > " + args.output + "data/bwa_R2.sam"))
        else:
            check_return(os.system(
                "bwa mem -t 32 " + args.output + "bwa-db/cds.fa " + args.output + "data/R1.fa > " + args.output + "data/bwa_R1.sam"))
            check_return(os.system(
                "bwa mem -t 32 " + args.output + "bwa-db/cds.fa " + args.output + "data/R2.fa > " + args.output + "data/bwa_R2.sam"))
        print("Running bioconvert")
        check_return(os.system(
            "bioconvert sam2paf " + args.output + "data/bwa_R1.sam " + args.output + "data/bwa_R1.paf"))
        check_return(os.system(
            "bioconvert sam2paf " + args.output + "data/bwa_R2.sam " + args.output + "data/bwa_R2.paf"))
        print("Running diamond_paf")
        check_return(os.system(
            "python3 paf_remove_unhit.py " + args.output + "data/bwa_R1.paf " + args.output + "data/bwa_R1.new.paf"))
        check_return(os.system(
            "python3 paf_remove_unhit.py " + args.output + "data/bwa_R2.paf " + args.output + "data/bwa_R2.new.paf"))

        print("Running bwa analysis")
        check_return(os.system(
            "python3 paf_result_analysis.py " + args.output + "data/bwa_R1.new.paf " + args.output + "data/bwa_R2.new.paf " + r1counts + " " + r2counts + " " + args.output + "result/bwa.sequence_FPKM.csv " + args.output + "result/bwa.cazyfamily_FPKM.csv"))
        bwa_end = time.time()
        print('Minimap Run Time(second): ' + str(bwa_end - bwa_start))

    if not args.source:
        os.system("rm -rf " + args.output + "source")

    stop = time.time()
    print('Total Run Time: ' + str(stop - start))


if __name__ == '__main__':
    main()
