#!/usr/bin/env python3
from Bio import SeqIO
import argparse
import glob, os

def GetArguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--fasta_folder", type=str, required=True, help="The folder that holds the NfracSeqRemover.py outputted fasta files.")
    parser.add_argument("-m", "--min_ind", type=int, required=True, help="The minimum number of sequences a file must hold in order to be outputted into the the output directory.")
    parser.add_argument("-o", "--outdir", type=str, required=True, help="The diretory were filtered output should be written into. Output will be renamed and will not overwrite any previous output.")
    args = parser.parse_args()
    return args

def MinIndivs(file, min, out):
    counter = 0
    for genomes in SeqIO.parse(file, 'fasta'):
        counter = counter + 1
    if counter => min:
        outfile = out+"/"+os.path.basename(file)[:-6]+"_"+str(counter)+"_Indivs.fasta"
        os.system("cp "+file+" "+outfile)

########################
def main():
    arg = GetArguments()
    for fasta in glob.glob(arg.fasta_folder+"/*.fasta"):
        MinIndivs(fasta, arg.min_ind, arg.outdir)

if __name__=='__main__':
    main()
