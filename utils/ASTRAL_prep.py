#!/usr/bin/env python3
import argparse
import re
from Bio import Phylo
import glob, os

def GetArguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--fasta_folder", type=str, required=True)
    parser.add_argument("-o", "--outfile", type=str, required=True)
    args = parser.parse_args()
    return args

def RedTip(tree):
    t = Phylo.read(tree, "newick")
    red = re.compile('_.*$')
    for leaf in t.get_terminals():
        name = leaf.name
        leaf.name = red.sub('',name)
    return(t)


########################
def main():
    arg = GetArguments()
    tlist = []
    for tree in glob.glob(arg.fasta_folder+"/*bestTree*"):
        tlist.append(RedTip(tree))
    Phylo.write(tlist, arg.outfile, "newick")


if __name__=='__main__':
    main()
