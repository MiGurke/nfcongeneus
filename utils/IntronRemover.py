#!/usr/bin/env python3
from Bio import SeqIO
import sqlite3
import argparse
import gffutils
import glob, os
from Bio import AlignIO
import sys

########################################
###USE BEFORE FRAC AND INDIV REMOVER!###
########################################
def GetArguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--anno", type=str, required=True)
    parser.add_argument("-f", "--fasta_folder", type=str, required=True)
    parser.add_argument("-o", "--output_folder", type=str, required=True)
    args = parser.parse_args()
    return args

def CreateAnnoDB(anno):
    DB = anno+".DB"
    gffutils.create_db(anno, DB, force = True, merge_strategy="create_unique")
    return DB

def ExonLister(DB,fasta):
    fname = os.path.basename(fasta)
    db = gffutils.FeatureDB(DB)
    start = fname[fname.index(":") + 1 : fname.index("-")]
    end = fname[fname.index("-") + 1 : -6]
    chrom = fname[0 : fname.index(":")]
    reg = db.region(region=fname[0 : -6], featuretype='gene')
    for i in reg:
        id = i.id
    elist = []
    for i in db.children(id, featuretype='exon'):
        se = [i.start,i.end]
        elist.append(se)
    return elist

def RmIntrons(elist,fasta):
    fname = os.path.basename(fasta)
    start = fname[fname.index(":") + 1 : fname.index("-")]
    end = fname[fname.index("-") + 1 : -6]
    try:
        alg = AlignIO.read(open(fasta), "fasta")
    except ValueError as e:
        sys.exit("Value Error: " + str(e) + "\n" + fname + " is empty. Please remove empty files from " + os.path.dirname(fasta) + "\n")
    exonalg = alg[:,0:0]
    for exon in elist:
        s = exon[0] - int(start)
        e = exon[1] - int(start)
        exonalg += alg[:,s:e]
        #alglist.append(alg[:,s:e])
        #print( start, end, exon, s, e, fname)
        #print(alg[:,s:e])
    return(exonalg)

def WriteAlg(alg, outfile):
    output = open(outfile, 'w')
    AlignIO.write(alg,output, 'fasta')
    output.close()

########################
arg = GetArguments()
db = CreateAnnoDB(arg.anno)
for fasta in glob.glob(arg.fasta_folder+"/*.fasta"):
    elist = ExonLister(db,fasta)
    exonalg = RmIntrons(elist,fasta)
    WriteAlg(exonalg, arg.output_folder + os.path.basename(fasta)[:-6] + "_IntronRM.fasta")
