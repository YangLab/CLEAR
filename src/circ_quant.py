#!/usr/bin/env python
"""
File: circ_quant.py
Author: Xu-Kai Ma
Date: 2018-05-31
Description: This file will calculate HPBcirc, HPBlinear and CIRCscore
Version: 0.1.1 # use stringTie
"""
import argparse
import pysam
import tempfile
import os.path
import shutil
from collections import defaultdict
import subprocess
import sys

import pybedtools
import spReads

#######################################################################
# Classes
######################################################################


class Extract_args():
    def __init__(self, bam, output,
                 anchor=2, intron=10):
        self.bam = bam
        self.output = output
        self.anchor = anchor
        self.intron = intron


class Annotate_args():
    def __init__(self, spRead, ss, output,
                 ignoreStrand=True, frag=True):
        self.spRead = spRead
        self.ss = ss
        self.output = output
        self.ignoreStrand = ignoreStrand
        self.frag = frag


#######################################################################
# Functions
######################################################################
def which(program):
    '''
    Check the path of external programs, and source codes are modified from
    https://github.com/infphilo/tophat/blob/master/src/tophat.py.
    '''
    def is_executable(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    for path in os.environ["PATH"].split(os.pathsep):
        progpath = os.path.join(path, program)
        if is_executable(progpath):
            return progpath
    return None

def get_read_length(AlignSeg):
    '''
        AlignSeg should be pysam.AlignedSegment object.
    '''
    length = AlignSeg.query_length
    if length > 0:
        return length

    length = AlignSeg.infer_read_length()

    return length


def count_mapped(bam_f, length_flag):
    # bam = pysam.AlignmentFile(bam_f)
    samfile = pysam.AlignmentFile(bam_f, 'rb')

    if not length_flag:

        idxstats = pysam.idxstats(bam_f)
        read_num = sum([int(line.split('\t')[2]) for
                        line in idxstats.split('\n')[:-1]])

        for read in samfile.fetch():
            read_length = get_read_length(read)

        total_bases = read_num * read_length
    else:
        total_bases = 0
        for read in samfile.fetch():
            # line = line.strip().split()
            total_bases += get_read_length(read)

    return total_bases


def cal_circ_hpb(circ_file, mapped_base_num, output):
    with open(circ_file, 'r') as circ_f,\
            open(output, 'w') as out:
        for line in circ_f:
            line = line.strip().split('\t')

            hpb = float(line[12]) * 1000000000 / mapped_base_num
            line.append(str(hpb))

            line = '\t'.join(line) + '\n'
            out.write(line)


def get_max_iso(ref_f, bam_f, output_d, output_f):

    if which('genePredToGtf') is None:
        sys.exit('genePredToGtf is required for maximal isoform selection!')
    if which('stringtie') is None:
        sys.exit('stringtie is required for maximal isoform selection!')

    rename_ref = open('{}/rename_ref.gp'.format(output_d), 'w')
    with open(ref_f, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            line = line[1:] + ['0', '{}:{}'.format(line[0], line[1])]
            line = '\t'.join(line) + '\n'

            rename_ref.write(line)
    rename_ref.close()

    rename_ref = '{}/rename_ref.gp'.format(output_d)
    rename_gtf = '{}/rename_ref.gtf'.format(output_d)
    gp_gtf_cmd = ['genePredToGtf', 'file', rename_ref, rename_gtf]
    subprocess.call(gp_gtf_cmd)

    gene_abund = '{}/gene_abund.txt'.format(output_d)
    assem_gtf = '{}/assem.gtf'.format(output_d)
    stringtie_cmd = ['stringtie', '-p', '5', '-e', '-G', rename_gtf,
                     '-A', gene_abund, '-o', assem_gtf, bam_f]
    subprocess.call(stringtie_cmd)

    gene_maxiso = defaultdict(lambda: ('guard', -1.0))
    with open(gene_abund, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            # remove first line
            if line[0] == 'Gene ID':
                continue

            gene = line[0].split(':')
            gene, iso = gene[0], gene[1]
            tpm = float(line[-2])

            if tpm > gene_maxiso[gene][1]:
                gene_maxiso[gene] = (iso, tpm)

    output = '{}/{}'.format(output_d, output_f)
    isos = set( gene_maxiso[i][0] for i in gene_maxiso )
    have_output = set()
    with open(ref_f, 'r') as refs,\
            open(output, 'w') as out:
        for line in refs:
            sp_line = line.strip().split('\t')
            # print(sp_line[1])
            if sp_line[1] in isos and sp_line[1] not in have_output:
                out.write(line)
                have_output.add(sp_line[1])


def get_sps(input, output):
    # most copy from getSPSite.py
    sp2name = defaultdict(list)

    with open(input, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            eCount = int(line[8])
            eStarts = [int(i) for i in line[9].split(',')[:-1]]
            eEnds = [int(i) for i in line[10].split(',')[:-1]]

            if eCount == 1:
                continue

            strand = line[3]
            if strand != '-':
                # not '-', then treat as +
                s_i, e_i, step = 0, eCount-1, 1
            else:
                s_i, e_i, step = eCount-2, -1, -1

            for i in range(s_i, e_i, step):
                sp = "{}\t{}\t{}\t{}".format(line[2], eEnds[i], eStarts[i+1], line[3])

                index = (i - s_i) * step
                name = "{}:{}:{}".format(line[0], line[1], index)

                sp2name[sp].append(name)

    with open(output, 'w') as f:
        for k, v in sp2name.items():
            line = '{}\t{}'.format(k, ','.join(v)) + '\n'
            f.write(line)


def rm_circ_sp(circ, in_sp, out_sp, selected_circ, hpb_thred=1):
    with open(circ, 'r') as circs,\
            open(selected_circ, 'w') as s_circs:
        for line in circs:
            sp_line = line.strip().split('\t')
            if float(sp_line[-1]) >= hpb_thred:
                s_circs.write(line)

    in_sp = pybedtools.BedTool(in_sp)

    in_sp.intersect(b=selected_circ, v=True).saveas(out_sp)


def cal_gene_hpb(csj_f, mapped_base_num, output):
    gene_hpb = defaultdict(list)
    with open(csj_f, 'r') as in_f:
        for line in in_f:
            line = line.strip().split('\t')

            gene = line[4].split(':')[0]
            hpb = float(line[5]) * 1000000000 / mapped_base_num
            gene_hpb[gene].append(hpb)

    with open(output, 'w') as out:
        for gene in gene_hpb:
            hpb = sum(gene_hpb[gene]) / len(gene_hpb[gene])

            line = '{}\t{}\n'.format(gene, hpb)
            out.write(line)


def cal_score(circ_hpb, gene_hpb, output, ratio, addition=0.01):
    gene_hpbs = defaultdict(float)
    with open(gene_hpb, 'r') as genes:
        for line in genes:
            line = line.strip().split('\t')
            gene, hpb = line[0], float(line[1])
            gene_hpbs[gene] = hpb

    with open(circ_hpb, 'r') as circs,\
            open(output, 'w') as out:
        for line in circs:
            line = line.strip()

            sp_line = line.split('\t')
            gene, c_hpb = sp_line[14], float(sp_line[-1])

            circ_score = ratio * (c_hpb + addition) /\
                (gene_hpbs[gene] + addition)

            line = '{}\t{}\t{}\n'.format(line,
                                         gene_hpbs[gene],
                                         circ_score)
            out.write(line)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--circ', dest='circ',
                        required=True,
                        help='Input circular RNA file from CIRCexplorer2.')
    parser.add_argument('-b', '--bam', dest='bam',
                        required=True,
                        help='Input mapped reads from HISAT2 in BAM format.')
    # for get max
    # parser.add_argument('-S', '--sl', dest='sl',
                        # required=True,
                        # help='stringTie mapping file')
    parser.add_argument('-r', '--ref', dest='ref',
                        required=True,
                        help='The refFlat format gene annotation file.')
    parser.add_argument('--threshold', dest='threshold',
                        default=1, type=float,
                        help='Threshold of HPB for choosing circRNAs to \
                        filter linear SJ.[default: 1]')
    parser.add_argument('--ratio', dest='ratio',
                        default=1, type=float,
                        help='The ratio is used for adjusting comparison between \
                        circ and linear.[default: 93.0/85]')

    parser.add_argument('-l', '--length', dest='length',
                        default=False, action='store_true',
                        help='Wether to consider all reads\' length? [default: False]')

    # parser.add_argument('-f', '--frag', dest='frag',
                        # default=False, action='store_true',
                        # help='Is mapped as paired-end? [default: False]')

    parser.add_argument('-t', '--tmp', dest='tmp',
                        default=False, action='store_true',
                        help='Keep tmp dir? [default: False]')

    parser.add_argument('-o', '--output', dest='output',
                        default='circRNA_quant.txt',
                        help='Output file. [default: circRNA_quant.txt]')

    args = parser.parse_args()
    print('###Parameters:')
    print(args)
    print('###Parameters')

    o_prefix = os.path.dirname(args.output)
    if o_prefix == '':
        o_prefix == '.'
    tmp_dir = tempfile.mkdtemp(dir=o_prefix,
                               prefix=os.path.basename(args.output))

    mapped_base_num = count_mapped(args.bam, args.length)

    cal_circ_hpb(args.circ, mapped_base_num,
                 output='{}/circ_hpb.txt'.format(tmp_dir))

    # print('to get max iso')
    # get maximum isoform splice sites
    get_max_iso(args.ref, args.bam,
                output_d=tmp_dir,
                output_f='max_iso.ref')

    # sys.exit('stop at get max iso')

    get_sps(input='{}/max_iso.ref'.format(tmp_dir),
            output='{}/max_iso_sps.txt'.format(tmp_dir))
    rm_circ_sp(circ='{}/circ_hpb.txt'.format(tmp_dir),
               in_sp='{}/max_iso_sps.txt'.format(tmp_dir),
               out_sp='{}/max_iso_sps_noCirc.txt'.format(tmp_dir),
               selected_circ='{}/selected_circ_hpb.txt'.format(tmp_dir),
               hpb_thred=args.threshold)

    # to calculate linear hpb
    extract_args = Extract_args(bam=args.bam,
                                output='{}/spReads.bed'.format(tmp_dir))
    spReads.extract(extract_args)
    annotate_args = Annotate_args(spRead='{}/spReads.bed'.format(tmp_dir),
                                  ss='{}/max_iso_sps_noCirc.txt'.format(tmp_dir),
                                  output='{}/csj_count.txt'.format(tmp_dir))
    spReads.annotate(annotate_args)

    cal_gene_hpb('{}/csj_count.txt'.format(tmp_dir),
                 mapped_base_num,
                 output='{}/gene_hpb.txt'.format(tmp_dir))

    # for CIRCscore
    cal_score(circ_hpb='{}/circ_hpb.txt'.format(tmp_dir),
              gene_hpb='{}/gene_hpb.txt'.format(tmp_dir),
              output=args.output,
              ratio=args.ratio)

    # print(mapped_num)
    if not args.tmp:
        shutil.rmtree(tmp_dir)


if __name__ == '__main__':
    main()
