#!/usr/bin/env python3
from Bio import SeqIO
import sqlite3
import argparse
import gffutils
import glob, os

def GetArguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--anno", type=str, required=True)
    parser.add_argument("-f", "--fasta_folder", type=str, required=True)
    args = parser.parse_args()
    return args

def CreateAnnoDB(anno):
    DB = anno+".DB"
    gffutils.create_db(anno, DB, force = True, merge_strategy="create_unique")
    return DB

def RmIntrons(DB,fasta):
    fname = os.path.basename(fasta)
    db = gffutils.FeatureDB(DB)
    start = fname[fname.index(":") + 1 : fname.index("-")]
    end = fname[fname.index("-") + 1 : -6]
    chrom = fname[0 : fname.index(":")]
    reg = db.region(region=fname[0 : -6], featuretype='gene')
    for i in reg:
        id = i.id
    for i in db.children(id, featuretype='exon'):
        print(i)
########################
arg = GetArguments()
db = CreateAnnoDB(arg.anno)
for fasta in glob.glob(arg.fasta_folder+"/*.fasta"):
    RmIntrons(db,fasta)
