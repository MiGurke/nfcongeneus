#!/usr/bin/env python3
from Bio import SeqIO
import argparse
import gzip

def GetArguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--fasta", type=str, required=True)
    parser.add_argument("-s", "--start", type=str, required=True)
    parser.add_argument("-e", "--end", type=str, required=True)
    parser.add_argument("-b", "--bamName", type=str, required=True)
    parser.add_argument("-c", "--chunk", type=str, required=True)
    args = parser.parse_args()
    return args

arg = GetArguments()

with gzip.open(arg.fasta, "rt") as handle:
    for genomes in SeqIO.parse(handle, 'fasta'):
        id = genomes.id
        seq = str(genomes.seq)
        print(">"+arg.bamName+"_"+arg.chunk+"\n"+seq[int(arg.start):int(arg.end)])
