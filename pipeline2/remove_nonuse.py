import sys


# argv1: file

def main():
    with open(sys.argv[1], "r+") as f:
        d = f.read()
        t = d.replace('(-)', '')
        f.seek(0, 0)
        f.write(t)
    f.close()
    with open(sys.argv[1], "r+") as f:
        d = f.read()
        t = d.replace('(+)', '')
        f.seek(0, 0)
        f.write(t)
    f.close()


if __name__ == '__main__':
    main()
