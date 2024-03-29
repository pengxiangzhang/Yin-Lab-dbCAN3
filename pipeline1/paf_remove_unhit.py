#!/usr/bin/python

import sys


# Input(argv1): a paf file
# Output(argv2): a paf file without un-hit

def remove_not_hit(path, new_path):
    file1 = open(path, 'r')
    f = open(new_path, 'w')
    for line in file1.readlines():
        x = line.split("	")
        if x[2] != "*" and x[4] != 255:
            f.write(line)
    f.close()


def main():
    remove_not_hit(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()
