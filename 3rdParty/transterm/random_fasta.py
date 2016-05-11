#!/usr/bin/env python

# generate a FASTA file with length random characters
import random, sys

def random_base(at):
    if random.random() < at:
        return "AT"[random.randint(0,1)]
    else:
        return "GC"[random.randint(0,1)]

def random_fasta(fasta, length, at):
    fasta = open(fasta, 'w')
    print(">random uniform length=", length, file=fasta)

    for i in range(1,length+1):
        fasta.write(random_base(at))
        if i % 60 == 0:
            fasta.write("\n")
    fasta.close()

def random_genes(cfile, length):
    coords = open(cfile, 'w')
    i = 1
    while i < length:
        i += random.randint(5, 500)
        j = i + random.randint(200, 2000)
        if(j < length):
            print("UNK", i, j, "random", file=coords)
        i = j
    coords.close()

# main program:
def main():
    at, length, file, coordfile = float(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
    random_fasta(file, length, at)
    random_genes(coordfile, length)

if __name__ == '__main__': main()
