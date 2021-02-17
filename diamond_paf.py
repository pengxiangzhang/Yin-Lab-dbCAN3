#!/usr/bin/python

import time
import sys

# Input(argv1): diamond generated paf file
# Output(argv2): diamond generated paf file without unhit

def remove_not_hit(path,new_path):
	file1 = open(path, 'r')
	f = open(new_path, 'w')
	for line in file1.readlines():
		x = line.split("	")
		if x[2]!="*" and x[4]!=255:
			f.write(line)
	f.close()


def main():
	start = time.time()
	reads_R1 = remove_not_hit(sys.argv[1],sys.argv[2])

	stop = time.time()
	print('run: ' + str(stop - start))

if __name__ == '__main__':
    main()