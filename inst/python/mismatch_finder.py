#!/usr/bin/env python

###############################################################################
# Find mismatches - SNVs, insertions, deletions
#   - Input: BAM with index, reference fasta, gene's chr, start, end, out_prefix
#   - Output: a textfile containing SNVs & indels info
#
# Command: ./mismatch_finder.py --BAM your.bam --fasta mm10.fa --chrom chr15 \
#                               --start 100 --end 2000 --outPrefix ./geneA_
# 
# by Chi-Lam Poon
################################################################################

import pysam
import gzip, argparse


class VariantFinder:
    def __init__(self, ref_fasta, ref_fasta_start):
        self.ref_fa = ref_fasta
        self.ref_fa_start = ref_fasta_start

    def get_mm_positions(self, cigartuples, reference_start):
        # NOTE: pysam are 0-base
        # concating coding sequence blocks
        seq_pos = []
        deletion_pos = []
        current_ref_pos = reference_start
        for i in range(len(cigartuples)):
            cigar_event = cigartuples[i][0]
            if cigar_event == 4:
                # soft clipping
                seq_pos.extend(['S'] * cigartuples[i][1])
            elif cigar_event == 5:
                # hard clipping
                continue
            elif cigar_event == 0:
                # mapped
                pos_list = list(range(current_ref_pos, current_ref_pos + cigartuples[i][1]))
                seq_pos.extend(pos_list)
                current_ref_pos += cigartuples[i][1]
            elif cigar_event == 3:
                # intron
                current_ref_pos += cigartuples[i][1]
            elif cigar_event == 1:
                # insertion to reference
                seq_pos.extend(['I'] * cigartuples[i][1])
            elif cigar_event == 2:
                # deletion from reference
                pos_list = list(range(current_ref_pos, current_ref_pos + cigartuples[i][1]))
                deletion_pos.extend(pos_list)
                current_ref_pos += cigartuples[i][1]
            else:
                raise('Found new cigar string: ' + str(cigar_event))
        return seq_pos, deletion_pos

    def find_mismatches(self, alignment, gene_start, gene_end):
        # gene_start, end: 1-base
        aln_ref_start = alignment.reference_start
        aln_ref_end = alignment.reference_end
        if aln_ref_start >= gene_end-500 or aln_ref_end <= gene_start+500:
            # no overlaps
            mismatches = '\t'.join(['NA'] * 3)
        else:
            snv = []
            deletion = []
            insertion = []
            seq = alignment.seq
            aln_ref_start = alignment.reference_start

            seq_pos, del_pos = self.get_mm_positions(alignment.cigartuples, aln_ref_start)
            # collect deletion nucleotides on reference
            for d in del_pos:
                if d < self.ref_fa_start or d > self.ref_fa_start + len(self.ref_fa):
                    continue
                nt = self.ref_fa[d - self.ref_fa_start].upper()
                deletion.append(str(d+1) + '_' + nt)

            for idx, s in enumerate(seq_pos):
                if s == 'S':
                    # soft clips skipped
                    continue
                elif s == 'I':
                    # collect insertion nucleotides on read & ref positions
                    read_nt = seq[idx]

                    pos_before = seq_pos[max(0, idx-1)]
                    pos_after = seq_pos[min(len(seq_pos)-1, idx+1)]
                    # if >1 insertions in a row
                    if pos_before != 'I':
                        last_inser_pos = pos_before
                    else:
                        pos_before = last_inser_pos if last_inser_pos else aln_ref_start

                    add_idx = 1
                    while pos_after == 'I':
                        if idx + 1 > len(seq_pos)-1:
                            pos_after = aln_ref_end
                            break
                        pos_after = seq_pos[min(len(seq_pos)-1, idx+1+add_idx)]
                        add_idx += 1

                    ref_pos = str(pos_before + 1) + '_' + str(pos_after + 1)
                    insertion.append(ref_pos + '_' + read_nt)
                elif gene_start - 1 <= s < gene_end:
                    # mapped seq - compare nucleotides on reference and read
                    read_nt = seq[idx].upper()
                    ref_nt = self.ref_fa[s - self.ref_fa_start].upper()
                    if read_nt != ref_nt:
                        variant = str(s+1) + '_' + ref_nt + '|' + read_nt
                        snv.append(variant)

            # organize info
            snv_on_read = ';'.join(snv) if snv else 'NA'
            insertion_on_read = ';'.join(insertion) if insertion else "NA"
            deletion_on_read = ';'.join(deletion) if deletion else "NA"
            mismatches = '\t'.join([snv_on_read, insertion_on_read, deletion_on_read])
        return mismatches

def main(args):
    bam_file = args['BAM'][0]
    fasta_file = args['fasta'][0]
    gene_chr = args['chrom'][0]
    gene_start = int(args['start'][0])
    gene_end = int(args['end'][0])
    gene=args['gene'][0]
    out_prefix = args['outPrefix'][0]

    # Open Sesame!
    bam = pysam.AlignmentFile(bam_file, 'rb')
    fasta = pysam.FastaFile(fasta_file)
    out_file = out_prefix + gene + '.mismatches.txt.gz'
    out = gzip.open(out_file, 'wb')

    # fetch
    sub_bam = bam.fetch(gene_chr, gene_start-10000, gene_end+10000) # assume longest read len is 10k
    ref_fasta = fasta.fetch(gene_chr, gene_start-1, gene_end)

    # find mismatches & write
    header = 'chrom\treadId\tSNV\tinsertion\tdeletion\n'
    out.write(header.encode())

    findr = VariantFinder(ref_fasta, gene_start-1)
    for aln in sub_bam:
        mm_info = findr.find_mismatches(aln, gene_start, gene_end)
        row = gene_chr + '\t' + aln.query_name + '\t' + str(mm_info) + '\n'
        out.write(row.encode())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mapping reads to snps')
    parser.add_argument('--BAM', required=True, nargs=1)
    parser.add_argument('--fasta', required=True, nargs=1)
    parser.add_argument('--chrom', required=True, nargs=1)
    parser.add_argument('--start', required=True, nargs=1)
    parser.add_argument('--end', required=True, nargs=1)
    parser.add_argument('--gene', required=True, nargs=1)
    parser.add_argument('--outPrefix', required=True, nargs=1)
    args = vars(parser.parse_args())
    main(args)
