import argparse
import os
import time

parser = argparse.ArgumentParser(description='pipeline2')
parser.add_argument('-b', '--bowtie2', dest='bowtie', action='store_true', help='Run Bowtie2?', default=False)
parser.add_argument('-d', '--diamond', dest='diamond', action='store_true', help='Run diamond?', default=False)
parser.add_argument('-m', '--minimap2', dest='minimap', action='store_true', help='Run minimap2?', default=False)
parser.add_argument('-a', '--bwa', dest='bwa', action='store_true', help='Run BWA?', default=False)
parser.add_argument('-s', '--source', dest='source', action='store_true', help='Store source file?', default=False)
parser.add_argument('-o', dest='output', type=str, help='Output location?', default="output/")
parser.add_argument('-r1', dest='input1', type=str, help='R1 File', required=True)
parser.add_argument('-r2', dest='input2', type=str, help='R2 File', required=True)
args = parser.parse_args()


# python3 pipeline2.py -b -d -m -a -s -o run -r1 /mnt/array2/pengxiang/annotation/reads/fq/example1/sample_R1_val_1.fq.gz -r2 /mnt/array2/pengxiang/annotation/reads/fq/example1/sample_R2_val_2.fq.gz
def check_return(status):
    if status != 0:
        print("Error")
        exit(1)


def main():
    if not args.output.endswith('/'):
        args.output = args.output + "/"
    print("Writing to: " + args.output)

    start = time.time()
    check_return(os.system("mkdir -p " + args.output))
    check_return(os.system("mkdir -p " + args.output + "data/"))
    check_return(os.system("mkdir -p " + args.output + "result/"))

    check_return(os.system("megahit -1 " + args.input1 + " -2 " + args.input2 + " -o " + args.output + "data/megahit"))
    check_return(os.system(
        "run_dbcan.py " + args.output + "data/megahit/final.contigs.fa prok --out_dir " + args.output + "data/dbcan_out --db_dir /mnt/array2/pengxiang/annotation/pipeline2/dbcan/db --dia_cpu 4  --hmm_cpu=4 --hotpep_cpu 4 --tf_cpu 4 --tf_cpu 4"))
    # TODO: Change db location

    check_return(os.system(
        "python3 filter_gff.py " + args.output + "data/dbcan_out/prodigal.gff " + args.output + "data/dbcan_out/overview.txt " + args.output + "data/new.gff"))
    check_return(os.system(
        "gff2bed < " + args.output + "data/new.gff | sed \"s/cds-//\" | awk \'{if ($8==\"CDS\") print $0}\' | bedtools getfasta -fi " + args.output + "data/megahit/final.contigs.fa -bed - -s -name > " + args.output + "data/bedtool.fa"))

    pipeline1 = "python3 pipeline2-1-besthit.py -r1 " + args.input1 + " -r2 " + args.input2 + " -o " + args.output + "result/ "

    if args.bowtie:
        check_return(os.system("mkdir " + args.output + 'result/bowtie2-index'))
        check_return(
            os.system("bowtie2-build " + args.output + "data/bedtool.fa " + args.output + "result/bowtie2-index/cds"))
        pipeline1 = pipeline1 + "-b "

    if args.diamond:
        check_return(os.system("mkdir " + args.output + 'result/diamond-db'))
        check_return(
            os.system("python3 RNA_translate.py " + args.output + "data/bedtool.fa " + args.output + "data/protein.fa"))
        check_return(os.system(
            "diamond makedb --in " + args.output + "data/protein.fa -d" + args.output + "result/diamond-db/cds"))  # cds.dmnd
        pipeline1 = pipeline1 + "-d "

    if args.minimap:
        check_return(os.system("mkdir " + args.output + 'result/minimap-db'))
        check_return(os.system("cp " + args.output + "data/bedtool.fa " + args.output + "result/minimap-db/cds.fa"))
        pipeline1 = pipeline1 + "-m "

    if args.bwa:
        check_return(os.system("mkdir " + args.output + 'result/bwa-db'))
        check_return(os.system(
            "cp " + args.output + "data/bedtool.fa " + args.output + "result/bwa-db/cds.fa;cd " + args.output + "result/bwa-db;bwa index cds.fa"))
        pipeline1 = pipeline1 + "-a "

    if not args.source:
        os.system("rm -rf " + args.output + "source")
    else:
        pipeline1 = pipeline1 + "-s "
    check_return(os.system("pwd"))
    print(pipeline1)
    check_return(os.system(pipeline1))

    if not args.source:
        os.system("rm -rf " + args.output + "source")

    end = time.time()
    print('Total Run Time: ' + str(end - start))


if __name__ == '__main__':
    main()

    # samtools sort -@ 32 -O bam -o R1.bam R1.sam
    # samtools sort -@ 32 -O bam -o R2.bam R2.sam
    # stringtie -e -B -p 8 -G /mnt/array2/pengxiang/annotation/pipeline2/new.gff  -o R1.gtf /mnt/array2/pengxiang/annotation/pipeline2/bowtie/R1.bam /mnt/array2/pengxiang/annotation/pipeline2/bowtie/R2.bam

# megahit -1 /mnt/array2/pengxiang/annotation/reads/fa/example0/sample_R1.fa -2 /mnt/array2/pengxiang/annotation/reads/fa/example0/sample_R2.fa -o megahit
# run_dbcan.py /mnt/array2/pengxiang/annotation/pipeline2/megahit/final.contigs.fa prok --out_dir dbcan_out --db_dir /mnt/array2/pengxiang/annotation/pipeline2/dbcan/db --dia_cpu 4  --hmm_cpu=4 --hotpep_cpu 4 --tf_cpu 4 --tf_cpu 4
# python3 filter_gff.py /mnt/array2/pengxiang/annotation/pipeline2/dbcan_out/prodigal.gff /mnt/array2/pengxiang/annotation/pipeline2/dbcan_out/overview.txt /mnt/array2/pengxiang/annotation/pipeline2/new.gff

# gff2bed < /mnt/array2/pengxiang/annotation/pipeline2/new.gff | sed "s/cds-//" | awk '{if ($8=="CDS") print $0}'| bedtools getfasta -fi /mnt/array2/pengxiang/annotation/pipeline2/megahit/final.contigs.fa -bed - -s -name > test
