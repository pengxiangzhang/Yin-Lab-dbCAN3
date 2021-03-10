#!/usr/bin/python

import sys


# Input(argv1): a paf file
# Output(argv2): a paf file only with R1
# Output(argv3): a paf file only with R2


def spit_minimap(input, output1, output2):
    file1 = open(input, 'r')
    r1 = open(output1, 'w')
    r2 = open(output2, 'w')
    for line in file1.readlines():
        i = line.split("	")
        j = i[0].split("/")
        if j[1] == "1":
            r1.write(line)
        if j[1] == "2":
            r2.write(line)
    r1.close()
    r2.close()


def main():
    spit_minimap(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == '__main__':
    main()
