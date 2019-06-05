#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import argparse
import multiprocessing
import dendropy
from collections import defaultdict
from Bio.Seq import translate
from Bio.Seq import reverse_complement
gene_list = [x.strip() for x in open(sys.argv[1], 'r')]

#fasta_list = open(sys.argv[2]).read().splitlines();print (fasta_list)
species_list = open(sys.argv[2]).read().splitlines();print (species_list)

fastadict = dict()

for species in species_list:
    fastafile = "1_assemblies/" + species + ".trinity.NI.fasta"
#    species = fastafile.strip().split('.trinity.NI')[0];print (species)
    with open(fastafile, 'r') as fasta_open:
        for line in fasta_open:
            if '>' in line:
                name = line.strip().split(' ')[0][1:]
                fastadict[name] = next(fasta_open).strip();
    last_file = "2_last/" + species + ".trinity.NI.LASTAL.TAB.Qlocal.plusblast.tophit.txt"
    lastfile_open = open(last_file, 'r')
    for line in lastfile_open:
        info = line.strip().split('\t');
	OHR = int((float(info[3]) / float(info[5]))*100)
        query = info[6];
        strand = info[9]
        start = int(info[7]);
        seq = fastadict[query]
        db = info[1]
        aalocal = info[14]
	E = info[25]
        seqcut = ''
        prev_neg_frameshift = "no"
        if db in gene_list:
            #print (start)
            CDS_filename = "3_homolog/" + db + ".homologs.fasta"
            with open(CDS_filename, "a+") as gen_fh:
                if strand == '-':
                    seq = reverse_complement(seq)
                    for b in aalocal:
                        if b == '-':
                            continue
                        if b == '\\':
                            end = start + 1
                            seqblock = seq[start:end]
                            #print (b,start," ",end,seqblock)
                            seqcut = seqcut+str(seqblock)+"nn"
                            start = start + 1
                        elif b == '/':
                            end = start + 2
                            seqblock = seq[start:end]
                            #print (b,start," ",end,seqblock)
                            seqcut = seqcut[:-3]+str(seqblock)+"nnnn"
                            start = start + 2
                            prev_neg_frameshift = "yes"
                        elif prev_neg_frameshift == "yes":
                            prev_neg_frameshift = "no"
                        else:
                            end = start + 3
                            seqblock = seq[start:end]
                            #print (b,start," ",end,seqblock)
                            seqcut = seqcut+str(seqblock)
                            start = start + 3
                else:
                    for b in aalocal:
                        if b == '-':
                            continue
                        if b == '\\':
                            end = start + 1
                            seqblock = seq[start:end]
                            #print (b,start," ",end,seqblock)
                            seqcut = seqcut+str(seqblock)+"nn"
                            start = start + 1
                        elif b == '/':
                            end = start + 2
                            seqblock = seq[start:end]
                            #print (b,start," ",end,seqblock)
                            seqcut = seqcut[:-3]+"nnnn"+str(seqblock)
                            start = start + 2
                            prev_neg_frameshift = "yes"
                        elif prev_neg_frameshift == "yes":
                            prev_neg_frameshift = "no"
                        else:
                            end = start + 3
                            seqblock = seq[start:end]
                            #print (b,start," ",end,seqblock)
                            seqcut = seqcut+str(seqblock)
                            start = start + 3
                #print (seqcut)
                outname = '>' + species + '~~' + query + '~~' + str(OHR) + '~~' + E
                outseq = seqcut
                print(outname, file=gen_fh)
                print(outseq, file=gen_fh)
                gen_fh.close()
        else:
            continue
       

