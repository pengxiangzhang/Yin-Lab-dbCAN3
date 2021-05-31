#!/usr/bin/python
import sys

from Bio.Seq import Seq


# Input(argv1): a gene file
# Output(argv2): a protein file

def translate(path, new_path):
    file1 = open(path, 'r')
    f = open(new_path, 'w')
    for line in file1.readlines():
        if line.startswith(">"):
            f.write(line)
        else:
            f.write(str(Seq(line).translate()))
            f.write("\n")

    f.close()


def main():
    translate(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()
