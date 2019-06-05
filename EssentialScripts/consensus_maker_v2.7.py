#!/usr/bin/env python3

from __future__ import print_function
import os
import sys
from collections import defaultdict
import re
import screed
from collections import Counter, OrderedDict
import itertools
import argparse
from operator import itemgetter 

__author__ = "lteasnail"
# IUPAC table
iupac_codes = {
    'A': 'A', 'G': 'G','C': 'C', 'T': 'T', 
    'K': 'GT', 'M': 'AC','S': 'CG', 'R': 'AG', 'W': 'AT', 'Y': 'CT', 
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG',
    'X': 'ACGT', 'N': 'ACGT', 'n': 'AGTC' 
}

# Unique combination of bases that can be found at one position sorted alphabetically, and their replacement base
# Note all ambiguous codes in lower case
code_replacements = {
    'A' : 'A', 'C' : 'C', 'G' : 'G', 'T' : 'T', 
    'AC' : 'm', 'AG' : 'r', 'AT' : 'w', 'CG' : 's', 'CT' : 'y', 'GT' : 'k',
    'ACG' : 'n', 'ACT' : 'n', 'AGT' : 'n', 'CGT' : 'n', 
    'ACGT' : 'n'
}

debug = 0

def parse_args():
    parser = argparse.ArgumentParser(description = 'Makes a consensus of contigs')
    parser.add_argument("infile", metavar = "infile",  type = str, help = "Fasta file to process")
    parser.add_argument("outfile", metavar = "outfile", type = str, help = "Output file for nucleotide consensus.")
    parser.add_argument("summaryfile", metavar = "summaryfile", type = str, help = "Output file for consensus summary.")
    parser.add_argument("statsfile", metavar = "statsfile", type = str, help = "Output file for consensus statistics.")
    parser.add_argument("-mismatches", dest = "mismatches", default = "0.02", type = float, help="The maximum p-distance between sequences which will be joined into a consensus. [default 0.02]")
    parser.add_argument("-header_spliter", dest="header_spliter", default = "~", type = str, help = "The character you want to split the header with. [default: ~]")
    parser.add_argument("-minlength", dest="minlength", default = "0.5", type = float, help = "Minimum length of contig as proportion to full sequence to retain [default: 0.5]")
    parser.add_argument("-dw_overlap", dest="dw_overlap", default = "15", type = int, help = "Minimum length of overlap for end downweighting. [default: 15]")
    parser.add_argument("-dw_factor", dest="dw_factor", default = "1", type = float, help = "Downweighting factor for overlapping ends. 1 is full, 0 is no downweight. [default: 1]")
    parser.add_argument("-debug", dest="debug", default = "0", type = int, help = "1 for debug mode, 0 for normal mode [default: 0]")
    parser.add_argument("-taxacount", dest="taxacount", default = "1", type = int, help = "Total number of taxa in full analysis")
    
    return parser.parse_args()

# Returns ungapped length of sequence
def get_seq_length(seq):
    countseq = 0
    gaps = ['-', '~']
    for i, nucleotide in enumerate(seq):
        if nucleotide not in gaps:
            countseq += 1
    return countseq

# Function to sort sequences by ungapped lengths
def sort_seqs_by_length(seqs_names):
    seqs = []
    for n, s in seqs_names:
        seqs.append((s, get_seq_length(s)))
    return list(i[0] for i in sorted(seqs, key=itemgetter(1), reverse=True))

# Create a consensus sequence
def create_consensus(seq1, seq2, dw_length, downweight, tail_dw):
    # Sequences must be strings, have the same length, and be aligned
    out_seq = ""
    countseq1 = 0
    countseq2 = 0
    gaps = ['-', '~']
    
    # Work out the downweight zones 3' end for each sequence (no gaps)
    dw_end1 = sum(1 for i in list(seq1) if i not in gaps) - dw_length
    dw_end2 = sum(1 for i in list(seq2) if i not in gaps) - dw_length
    if debug==1 : print ("End1:", dw_end1, "End2:", dw_end2, ", remove tail:", downweight)
    if debug==1 : print ("Seq1:", seq1, "\nSeq2:", seq2)
    
    for i, nucleotide in enumerate(seq1):
        couple = [nucleotide, seq2[i]]
        if couple[0] not in gaps:
            countseq1 += 1
        if couple[1] not in gaps:
            countseq2 += 1
            
        # If sequences overlap, and we have decided to downweight, remove bases only at the ends
        if downweight and (couple[0] not in gaps and couple[1] not in gaps):
            if (countseq1 <= dw_length and tail_dw[0] == 1) or (countseq1 >= dw_end1 and tail_dw[1] == 1): couple[0] =  '-' 
            if (countseq2 <= dw_length and tail_dw[2] == 1) or (countseq2 >= dw_end2 and tail_dw[3] == 1): couple[1] =  '-' 
 
        # If sequence 1 is a gap, return sequence 2
        if couple[0] in gaps:
            out_seq += couple[1]
        # If sequence 2 is a gap, return sequence 1 
        elif couple[1] in gaps:
            out_seq += couple[0]
        # If they are equal, return sequence 1
        elif couple[0].upper() == couple[1].upper():
            out_seq += couple[0]
        # Bases are not gaps and they are not equal
        else:
            # Convert each to unambiguous nucleotides
            n1 = list(iupac_codes[couple[0].upper()])
            n2 = list(iupac_codes[couple[1].upper()])
            # Create a sorted string of all unique nucleotides
            n3 = sorted(set(n1+n2))
            nucs = ''.join(str(n) for n in n3)
            # Replace with our IUPAC code
            replacement = code_replacements[nucs.upper()]
            if debug == 1 : print("position:", countseq1, couple[0], couple[1], nucs, replacement) 
            out_seq += replacement
    if debug == 1: print("Cons:", out_seq) 
    return out_seq

# Create a list of species and contigs
def make_header_species_dict(fasta_file, header_spliter='~'):
    taxa = OrderedDict()
    for seq in screed.open(fasta_file):
        seq_sample = seq.name.split(header_spliter)[0]
        rest_of_header = header_spliter.join(seq.name.split(header_spliter)[1:])
        # Store contig names and sequences in a structure by taxa
        if seq_sample in taxa:
            taxa[seq_sample].append((rest_of_header, seq.sequence))
        else:
            taxa[seq_sample] = [(rest_of_header, seq.sequence), ]
    return taxa

# Compare two sequences to determine whether we can create a consensus
def compare_seqs(seq_1, seq_2, mismatches, dw_length, factor, tail_dw):
    num_missing = 0
    size_overlap = 0
    gaps = ['-', '~']
    
    # Count all non-gap bases for each sequence, work out where 
    # downweighting zone starts from end of sequence
    seq_1_split = list(seq_1)
    dw_end1 = sum(1 for i in seq_1 if i not in gaps) - dw_length
    seq_2_split = list(seq_2)
    dw_end2 = sum(1 for i in seq_2 if i not in gaps) - dw_length

    # Work out number of differences - NOTE: sequences must be same length
    length = len(seq_1)
    num_diff = 0
    num_dw_diff = 0
    countseq1 = 0
    countseq2 = 0
    for base in range(len(seq_1_split)):
        # Keep count of where we are for each sequence
        if seq_1_split[base] not in gaps:
            countseq1 += 1
        if seq_2_split[base] not in gaps:
            countseq2 += 1
        # Count the number of missing bases from either sequence  
        if seq_1_split[base] in gaps or seq_2_split[base] in gaps:
            num_missing += 1
        # Else if bases do not match, add to differences   
        elif seq_1_split[base].upper() != seq_2_split[base].upper():
            num_diff += 1
            size_overlap += 1
            # Add to total differences in DW zone
            if countseq1 <= dw_length or countseq1 >= dw_end1 or countseq2 <= dw_length or countseq2 >= dw_end2: num_dw_diff += 1 
            if countseq1 <= dw_length: tail_dw[0]=1
            if countseq1 >= dw_end1: tail_dw[1]=1
            if countseq2 <= dw_length: tail_dw[2]=1
            if countseq2 >= dw_end2: tail_dw[3]=1
        else:
            size_overlap += 1
     
    # Calculate and print the p-distance
    num_diff = float(num_diff)
    num_dw_diff = float(num_dw_diff)
    size_overlap = float(size_overlap)
    if debug == 1 : print("size overlap:", size_overlap, ", num diff:", num_diff, ", num_missing:", num_missing, ", num_dw_diff:", num_dw_diff)
    
    # Only make consensus at this stage if there is an overlap and mismatches within allowed range
    if size_overlap > 0:
        pdistance = num_diff/(length - num_missing)
        if debug == 1  : print("Original p_dist:", pdistance, ",", float(mismatches))
        if pdistance <= float(mismatches):
            # Can create consensus, within allowed mismatches
            return 1
        elif num_dw_diff > 0:
            # Try downweighting
            pdistance = (num_diff - (num_dw_diff*factor))/(length - num_missing)
            if debug == 1  : print("Downweighted p_dist:", pdistance, ",", float(mismatches), "remove tail:", True)
            if pdistance <= float(mismatches):
                return 2
            else:
                return -1
        else:
            # Cannot create consensus, more than allowed mismatches
            return -1
    else:
        # No overlap
        return 0

# After merging overlapping contigs, merge non-overlapping contigs
# ONLY if there are NO overlapping contigs for this gene and taxon
def determine_non_overlapping(seqnames, mismatches, overlap, factor):
    seqs = sort_seqs_by_length(seqnames)
    overlapping = 0
    if len(seqs) > 1:
        # Sort sequences in set by ungapped length
        for lhs, rhs in itertools.combinations(seqs, 2):
            if debug == 1 : print("Determine Nonoverlapping Comparing seq1 vs seq2")
            # For each pair of sequences in the current set of sequences
            if compare_seqs(lhs, rhs, mismatches, overlap, factor, [0,0,0,0]) != 0:
                overlapping = 1
                break
            
    # If we found no overlapping contigs, merge all into one sequence
    if len(seqs) > 1 and overlapping == 0: 
        name = '~grp1uc' + str(len(seqs))
        seq = merge_non_overlapping(seqnames)
        return [(name,seq)]
    # Else overlap found, return original sequences    
    else:
        return seqnames

# Merge non-overlapping contigs
def merge_non_overlapping(names_and_seqs):
    seqs = sort_seqs_by_length(names_and_seqs)
    while True:
        # For each pair of sequences in the current set of sequences, sorted by ungapped length
        for lhs, rhs in itertools.combinations(seqs, 2):
            if debug == 1 : ("Merging non overlapping seq1:", lhs , " vs seq2:", rhs)
            new_seq = create_consensus(lhs, rhs, 0, 1, [])
            # Remove the original two sequences from the sequence set, and add the new one. 
            seqs.remove(lhs)
            seqs.remove(rhs)
            seqs.append(new_seq)
            break
        else:
            break
    return new_seq

# Determine whether to make a consensus sequence from pairs of contigs by gene/taxon
def decide_consensus(names_and_seqs, mismatches=0.02, overlap=15, factor=1):
    seqnames = {}
    newseqcount = {}
    seqs = sort_seqs_by_length(names_and_seqs)
    for n, s in names_and_seqs:
        seqnames.update({s:n})

    while True:
        if len(seqs) < 2:
            # We've only got 1 seq, so bail.
            if debug == 1 : print("Only have one sequence left, no consensus necessary")            
            break
        for lhs, rhs in itertools.combinations(seqs, 2):
            if debug == 1 : print("Decide_consensus: Comparing seq1:", seqnames.get(lhs), "seq2:", seqnames.get(rhs))
            # For each pair of sequences in the current set of sequences
            tail_dw=[0,0,0,0]
            result = compare_seqs(lhs, rhs, mismatches, overlap, factor, tail_dw)
            downweight = True if result == 2 else False
            if result > 0:
                # If we can join this pair, do so
                new_seq = create_consensus(lhs, rhs, overlap, downweight, tail_dw)
                # Remove the original two sequences from the sequence set, and
                # add the new one. 
                seqs.remove(lhs)
                if lhs in seqnames:
                    seqnames.pop(lhs)
                seqs.remove(rhs)
                if rhs in seqnames:
                    seqnames.pop(rhs)
                # Add the new long sequence to the top of the list
                seqs = [new_seq]+seqs;
   
                # Keep track of how many contigs have been merged in this sequence
                if lhs in newseqcount:
                    contig_count = newseqcount.get(lhs)+1
                    newseqcount.pop(lhs)
                elif rhs in newseqcount:
                    contig_count = newseqcount.get(rhs)+1
                    newseqcount.pop(rhs)
                else:
                    contig_count = 2 
                newseqcount.update({new_seq:contig_count})
                if debug == 1 : print("Made consensus sequence at contig count:", contig_count)
                
                # Break from the for loop, not the while loop. i.e., we want to
                # restart the "for all pairs" loop every time we join two seqs,
                # so that we compare the joined seq with the others.
                break
        else:
            # If we get the whole way though the "for each pairs" loop without
            # break-ing, then there's nothing that can be joined, so exit the
            # while loop even though len(seqs) > 1.
            break
    new_names_and_seqs = []
    consensus = 1
    for seq in seqs:
        name = seqnames.get(seq, '~grp'+ str(consensus) + "uc" + str(newseqcount.get(seq))) 
        new_names_and_seqs.append((name,seq))
        consensus += 1
    
    # Once we have finished merging overlapping contigs, attempt to merge non-overlapping    
    merged_names_and_seqs = determine_non_overlapping(new_names_and_seqs, mismatches, overlap, factor)
    
    return merged_names_and_seqs

def main():
    args = parse_args()
    while args.infile is None:
        args.infile = raw_input("\nYou have to specify an input file: ")
        if not os.path.exists(args.infile):
            args.infile = None
    while args.outfile is None:
        args.outfile = raw_input("\nYou have to specify an output file: ")
    while args.summaryfile is None:
        args.summaryfile = raw_input("\nYou have to specify a summary output file: ")
    while args.statsfile is None:
        args.statsfile = raw_input("\nYou have to specify a statistics output file: ")
    while args.taxacount is None or not isinstance(args.taxacount, int):
        args.taxacount = raw_input("\nYou have to specify the total number of taxa in the analysis: ")
            
    taxa_seqs = make_header_species_dict(args.infile, args.header_spliter)
    gene_idd = os.path.basename(args.infile).strip().split('.')[0]
    global debug
    debug = args.debug

    # Lists to keep track of all contigs for all taxa
    taxa_count = 0
    taxa_multicontig = 0
    cds_list = []

    stats_fh = open(args.statsfile, "a+")
    # Process contigs for this gene
    with open(args.outfile, "w") as gen_fh:
        for taxa in taxa_seqs:
            taxa_count += 1
            contig_list = []
            print("Processing taxa:", taxa, "for gene:", gene_idd)
            # Decide on consensus sequences
            consensus_seqs = decide_consensus(taxa_seqs[taxa], args.mismatches, args.dw_overlap, args.dw_factor)
            # Print out consensus sequences
            for name, seq in consensus_seqs:
                # Collect statistics on consensus sequences
                seq_length = get_seq_length(seq)
                cds_list.append(seq_length)
                contig_list.append(seq_length)
                # Print consensus contig to gene sequences file
                print(">", taxa, '~', name, "\n", seq, sep="", end='\n', file=gen_fh)
            # If this taxon has more than one contig for this gene, add to count
            if len(contig_list) > 1:
                taxa_multicontig += 1
            # Print total number of sequences for this taxa to statistics file
            print(gene_idd, taxa, len(contig_list), max(contig_list), sep="\t", end="\n", file=stats_fh)
    # Print total number of taxa and contigs for this gene to summary file
    with open(args.summaryfile, "a+") as summ_fh:
        print(gene_idd, taxa_count, len(cds_list), taxa_count/len(cds_list), taxa_multicontig, max(cds_list), (sum(cds_list)/max(cds_list))/args.taxacount, sep="\t", end='\n', file=summ_fh)
        
main()
