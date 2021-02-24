#!/usr/bin/python

import time
import sys


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
    start = time.time()
    spit_minimap(sys.argv[1], sys.argv[2], sys.argv[3])
    stop = time.time()
    print('spit_minimap run: ' + str(stop - start))


if __name__ == '__main__':
    main()
