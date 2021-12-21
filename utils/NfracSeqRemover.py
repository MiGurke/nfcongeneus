#!/usr/bin/env python3
from Bio import SeqIO
import argparse
import glob, os
import sys

def GetArguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--fasta_folder", type=str, required=True)
    parser.add_argument("-m", "--min_frac", type=float, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)
    args = parser.parse_args()
    return args

def CalcNfrac(seq, id):
    N = seq.count("N")
    n = seq.count("n")
    l = len(seq)
    try:
        frac = (n + N) / l
    except ZeroDivisionError:
        sys.exit(id + " has length 0. Please have a look at this file. Exiting...")
    return(frac)

def RemoveSeqs(file, minfrac, out):
    outfile = out+"/"+os.path.basename(file)[:-6]+"_"+str(minfrac)+"_NfracRM.fasta"
    out = open(outfile, "a")
    for genomes in SeqIO.parse(file, 'fasta'):
        id = genomes.id
        seq = genomes.seq
        frac = CalcNfrac(seq, id)
        if frac < minfrac:
            out.write(">"+id+"_Nfrac:"+str(frac)+"\n"+str(seq)+"\n")
    out.close

########################
def main():
    arg = GetArguments()
    for fasta in glob.glob(arg.fasta_folder+"/*.fasta"):
        RemoveSeqs(fasta, arg.min_frac, arg.outdir)

if __name__=='__main__':
    main()
