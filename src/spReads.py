#!/usr/bin/env python
"""
File: spReads.py
Author: Xu-Kai Ma
Date: 2017-04-23
Description: This file can be used to get linear junction reads
Version: 0.1.0
"""
import argparse
import pysam
from collections import defaultdict

#######################################################################
# Classes
######################################################################

#######################################################################
# Functions
######################################################################
def getCigar(read):
    '''
    get cigars
    '''
    num2op = {0: 'M', 1: 'I', 2: 'D',
              3: 'N', 4: 'S', 5: 'H',
              6: 'P', 7: '=', 8: 'X'}
    cigar = [ (num2op[num], length) for num, length in read.cigartuples]

    return cigar

def isSpRead(cigar):
    spFlag = False
    for op, length in cigar:
        if op == 'N':
            spFlag = True
    return spFlag

def parseJunction(cigar):
    '''
    Use cigar to get block count, block sizes, block offsets
    '''
    bCount = 0
    bSizes = []
    bOffsets= []
    introns = [] # for filter intron length

    bSize = 0
    bOffset = 0

    newCigar = cigar + [('N', 0)] # for the last block, and not change cigar
    for op, length in newCigar:
        if op in 'SHI':
            continue
        elif op in '=MXD':
            bSize += length
        elif op in 'N':
            bCount += 1
            bSizes.append(bSize)
            bOffsets.append(bOffset)

            # For last one
            if length != 0:
                introns.append(length)
            else:
                # this is the last one, mannually added
                break

            # prepare for next block
            bOffset += bSize
            bOffset += length
            bSize = 0

    return bCount, bSizes, bOffsets, introns

def extract(args):
    bam = pysam.AlignmentFile(args.bam)

    with open(args.output, 'w') as f:
        for read in bam.fetch():
            cigar = getCigar(read)

            # whether is spReads
            if not isSpRead(cigar):
                continue

            chrom = read.reference_name
            start = read.reference_start # it is 0-based
            end = read.reference_start + read.reference_length
            name = read.query_name
            strand = '+' if read.is_reverse else '-'
            bCount, bSizes, bOffsets, introns = parseJunction(cigar)

            # filters
            filterFlag = False
            for bSize in bSizes:
                if bSize < args.anchor:
                    filterFlag = True
            for intron in introns:
                if intron < args.intron:
                    filterFlag = True
            if filterFlag:
                continue

            # for output
            line = [chrom, str(start), str(end), name, '0', strand,
                    str(start), str(end), '0,0,0', str(bCount),
                    ','.join(str(i) for i in bSizes) + ',',
                    ','.join(str(i) for i in bOffsets) + ',']
            line = '\t'.join(line) + '\n'
            
            f.write(line)

def annotate(args):
    with open(args.ss, 'r') as ss, \
         open(args.spRead, 'r') as spReads, \
         open(args.output, 'w') as output:

        sp2reads = defaultdict(list)
        for line in spReads:
            line = line.strip().split('\t')
            chrom = line[0]
            start = int(line[1])
            strand = line[5]
            name = line[3]

            bCount = int(line[9])
            bSizes = [int(i) for i in line[10].split(',')[:-1]]
            bOffsets = [int(i) for i in line[11].split(',')[:-1]]

            for i in range(bCount - 1):
                sp_s = start + bOffsets[i] + bSizes[i]
                sp_e = start + bOffsets[i+1]

                if not args.ignoreStrand:
                    readSp = '{}\t{}\t{}\t{}'.format(chrom, sp_s, sp_e, strand)
                else:
                    readSp = '{}\t{}\t{}'.format(chrom, sp_s, sp_e)

                sp2reads[readSp].append(name)

        for line in ss:
            lineList = line.strip().split('\t')

            fourCols = '\t'.join(lineList[0:4])

            if not args.ignoreStrand:
                sp = '\t'.join(lineList[0:4])
            else:
                sp = '\t'.join(lineList[0:3])

            if sp in sp2reads:
                if args.frag:
                    num = len(set(sp2reads[sp]))
                else:
                    num = len(sp2reads[sp])

                o_line = '{}\t{}\t{}\t{}\n'.format(fourCols, lineList[4], num, 
                                                ','.join(sp2reads[sp]))
            else:
                num = 0
                o_line = '{}\t{}\t{}\n'.format(fourCols, lineList[4], num)

            output.write(o_line)

            # if not args.ignoreStrand:
                # sp = '\t'.join(lineList[0:4])

                # if sp in sp2reads:
                    # if args.frag:
                        # num = len(set(sp2reads[sp]))
                    # else:
                        # num = len(sp2reads[sp])

                    # o_line = '{}\t{}\t{}\t{}\n'.format(sp, lineList[4], num, 
                                                    # ','.join(sp2reads[sp]))
                # else:
                    # num = 0
                    # o_line = '{}\t{}\t{}\n'.format(sp, lineList[4], num)

            # else:
                # sp = '\t'.join(lineList[0:3])

                # if sp in sp2reads:
                    # if args.frag:
                        # num = len(set(sp2reads[sp]))
                    # else:
                        # num = len(sp2reads[sp])
                    # o_line = '{}\t{}\t{}\t{}\t{}\n'.format(sp, lineList[3], 
                                                           # lineList[4], num, 
                                                           # ','.join(sp2reads[sp]))
                # else:
                    # num = 0
                    # o_line = '{}\t{}\t{}\t{}\n'.format(sp, lineList[3],
                                                       # lineList[4], num)


#######################################################################
# Main
######################################################################
def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    p_extract = subparsers.add_parser('extract', 
                                    help='extract splicing reads')
    p_extract.add_argument('-b', '--bam', dest='bam', 
                           required=True,
                           help='min anchor length in the every anchor')
    p_extract.add_argument('-a', '--anchor', dest='anchor', 
                           type=int, default=10,
                           help='min anchor length in the every anchor')
    p_extract.add_argument('-i', '--intron', dest='intron', 
                           type=int, default=50,
                           help='min intron length')
    p_extract.add_argument('-o', '--output', dest='output',
                           default='sp_reads.bed',
                           help='output file')
    p_extract.set_defaults(func=extract)

    p_annotate = subparsers.add_parser('annotate',
                                       help='annotate the splicing reads')
    p_annotate.add_argument('--is', dest='ignoreStrand',
                            action='store_true',
                            help='whether ignore the strand information.')
    p_annotate.add_argument('-r', '--spRead', dest='spRead',
                            help='spRead file, the output of spReads.py extract')
    p_annotate.add_argument('-s', '--ss', dest='ss',
                            help='splice site file, format should be :\
                            chrom, sp start(intron), sp end(intron), name')
    p_annotate.add_argument('-o', '--output', dest='output',
                           default='sp_with_reads.txt',
                           help='output file')
    p_annotate.add_argument('-f', '--frag', dest='frag',
                           action='store_true',
                           help='For one sp site, if two paired reads are\
                           both overlapped, then just count one')
    p_annotate.set_defaults(func=annotate)

    args = parser.parse_args()
    print(args)
    args.func(args)


if __name__ == '__main__':
    main()

