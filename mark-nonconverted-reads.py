#!/usr/bin/env python3
"""
Script to take a Sam/Bam file (streamed though stdin, by default) of a bisulfite-treated
library, and filter out 'unconverted' reads if they have 3 or more unconverted Cs in a
non-CpG context.

Unconverted reads are flagged as failing platform / vendor quality checks (0x200) and have
the tag "XX:Z:UC" added at the end.

Note: One probably doesn't want to do this for non-human organisms, as they may have
significant levels of non-CpG methylation!
"""


import sys
import argparse
import pysam
import re

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required = True, help = "Reference fasta file")
    parser.add_argument("--bam", required = False, help = "Bam or sam file (must end in .bam or .sam)")
    parser.add_argument("--out", required = False, help = "Name for output bam file")

    args = parser.parse_args()
    return args


def parse_fasta(reference):
	"""
	Parses the reference fasta file, and reads it into a dictionary. For the header, will
	strip anything after the first space. For example:
	
	>chr1 AC:XXX gi:XXX LN:XXX
	
	will be saved to the dictionary as:
	chr1
	"""

	fasta_dict = {}
	fasta = open(reference, "r")
	for line in fasta:

		if line.startswith(">"): # Header line
			header = line.strip().split(" ")[0][1:] # Remove '>' and anything after ' '
			fasta_dict[header] = []

		else: # Sequence line
			fasta_dict[header].append(line.strip())
	fasta.close()	

	# If it's a multiline fasta, join the individual lines to one sequence
	for header in fasta_dict:
		fasta_dict[header] = "".join(fasta_dict[header])

	return fasta_dict
	

def parse_bam(bam_file, fasta_dict):
    """
    Reads through each alignment in the bam file, if it's part of a proper pair then checks
    for non-CpG Cs. If there are 3 or more, calls filter_snps() to check if they are
    unconverted Cs, or some snp / misalignment
    """

    unconverted_count = 0
    for read in bam_file:
    
        # Ignore reads not part of a proper pair
        if not read.is_proper_pair or read.is_qcfail or read.is_duplicate or read.is_secondary \
        or read.is_supplementary:
            out.write(read)
        
        elif (read.is_reverse and read.is_read2) or (read.mate_is_reverse and read.is_read1): # 'Top' strand
            # If there are 3 or more non-CpG Cs in the read, call filter_snps
            if read.query_alignment_sequence.count("C") - read.query_alignment_sequence.count("CG") >= 3:
                unconverted_count += filter_snps(read, fasta_dict)
            else:
                out.write(read)
        
        elif (read.is_reverse and read.is_read1) or (read.mate_is_reverse and read.is_read2): # 'Bottom' strand
            # If there are 3 or more non-CpG Gs in the read, call filter_snps
            if read.query_alignment_sequence.count("G") - read.query_alignment_sequence.count("CG") >= 3:
                unconverted_count += filter_snps(read, fasta_dict)
            else:
                out.write(read)

    return unconverted_count

def filter_snps(read, fasta_dict):
    """
    For reads with 3 or more unconverted Cs, check to be sure they actually align to a C
    on the reference genome (so are not mutations, snps, misalignments, etc.)
    """

    # For all bases in read, gets genomic coordinate of alignment
    coords = read.get_reference_positions(full_length = True)
    chrom = read.reference_name

    # If the alignment was softclipped at all, call softclip() to remove those indices
    # from the coords list
    if coords[0] is None or coords[-1] is None:
        coords = softclip(coords) # Replacing coords with a subset of itself

    my_unconverted = 0
    if (read.is_reverse and read.is_read2) or (read.mate_is_reverse and read.is_read1): # 'Top' strand
        for coord in coords:

            # If the base on the read is a C and not part of an indel
            if read.query_sequence[coords.index(coord)] == "C" and coord != None:
                
                # If the C is not part of a CpG
                if coords.index(coord) < len(read.query_sequence) - 1 and \
                   read.query_sequence[coords.index(coord) + 1] != "G":
                    
                    # If the C maps to a C on the genome (it's unconverted)
                    if fasta_dict[chrom][coord] == "C":
                        my_unconverted += 1
               
                # If the C is the last base of the read
                elif coords.index(coord) == len(read.query_sequence) - 1:
                    
                    # If the C maps to a C on the genome (it's unconverted)
                    if fasta_dict[chrom][coord] == "C":
                        my_unconverted += 1

    elif (read.is_reverse and read.is_read1) or (read.mate_is_reverse and read.is_read2): # 'Bottom' strand
        for coord in coords:
            
            # If the base on the read is a G and not part of an indel
            if read.query_sequence[coords.index(coord)] == "G" and coord != None:
                
                # If the G is not part of a CpG
                if coords.index(coord) > 0 and read.query_sequence[coords.index(coord) - 1] != "C":
                    
                    # If the G maps to a G on the genome (it's unconverted)
                    if fasta_dict[chrom][coord] == "G":
                        my_unconverted += 1
                
                # If the G is the first base of the read
                elif coords.index(coord) == 0:
                    
                    # If the G maps to a G on the genome (it's unconverted)
                    if fasta_dict[chrom][coord] == "G":
                        my_unconverted += 1

    # If there are 3 or more unconverted Cs on the read, set some flags and return
    if my_unconverted >= 3:
        read.set_tags(read.get_tags() + [("XX", "UC")])
        read.flag += 512
        out.write(read)
        return 1

    out.write(read)
    return 0

def softclip(coord_list):
    """
    If the alignment coords list includes softclipped bases ('None' in the list), cuts them
    from the beginning / end of the list
    """

    start = None
    end = None
    
    x = 0
    # Find first non-'None' index in the list
    while coord_list[x] is None:
        x += 1
    start = x
    
    x = len(coord_list) - 1
    # Find the last non-'None' index in the list
    while coord_list[x] is None:
        x -= 1
    end = x

    # The new coord list goes from the first - last aligned base
    new_coords = coord_list[start:end + 1]
    return new_coords



if __name__ == "__main__":
    args = argparser()

    fasta_dict = parse_fasta(args.reference)

    if args.bam:
        if args.bam.endswith(".bam"):
            mysam = pysam.AlignmentFile(args.bam, "rb")
        elif args.bam.endswith(".sam"):
            mysam = pysam.AlignmentFile(args.bam, "r")
    else:
        mysam = pysam.AlignmentFile("-", "r")

    if args.out:
        out = pysam.AlignmentFile(args.out, "wb", template = mysam)
    else:
        out = pysam.AlignmentFile("-", "wb", template = mysam)

    header = mysam.header.to_dict()

    unconverted = parse_bam(mysam, fasta_dict)

    header['PG'].append({"ID": "Find unconversions", "DS": unconverted})
