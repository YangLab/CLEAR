#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import subprocess
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

import pysam

import circ_quant

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

def hisat_align(m1, m2, index, gtf, hisat_dir, thread):
    # setup hisat mapping environment
    sp_flag = False
    pair_flag=False
    if m2:
        pair_flag=True

    if not which('hisat2'):
        sys.exit('Please install hisat2 at first.')

    if not which('hisat2_extract_splice_sites.py'):
        print('Without hisat2_extract_splice_sites.py, hisat2 will not use --known-splicesite-infile option.')
    else:
        sp_flag = True
        get_sp_sites = ['hisat2_extract_splice_sites.py',
                        gtf]
        with open('{}/sp.txt'.format(hisat_dir), 'wb') as f:
            print('# start to get sp sites for hisat mapping')
            sps = subprocess.check_output(get_sp_sites)
            f.write(sps)

    # for mapping steps
    if sp_flag and not pair_flag:
        hisat_map = ['hisat2',
                    '--no-softclip',
                    '--score-min', 'L,-16,0',
                    '--mp', '7,7',
                    '--rfg', '0,7',
                    '--rdg', '0,7',
                    '--dta',
                    '-k', '1',
                    '--max-seeds', '20',
                    '-p', thread,
                    '-x', index,
                    '--known-splicesite-infile', '{}/sp.txt'.format(hisat_dir),
                    '-U', m1,
                    '-S', '{}/align.sam'.format(hisat_dir)]
    elif not sp_flag and not pair_flag:
        hisat_map = ['hisat2',
                    '--no-softclip',
                    '--score-min', 'L,-16,0',
                    '--mp', '7,7',
                    '--rfg', '0,7',
                    '--rdg', '0,7',
                    '--dta',
                    '-k', '1',
                    '--max-seeds', '20',
                    '-p', thread,
                    '-x', index,
                    '-U', m1,
                    '-S', '{}/align.sam'.format(hisat_dir)]
    elif sp_flag and pair_flag:
        hisat_map = ['hisat2',
                    '--no-softclip',
                    '--score-min', 'L,-16,0',
                    '--mp', '7,7',
                    '--rfg', '0,7',
                    '--rdg', '0,7',
                    '--dta',
                    '-k', '1',
                    '--max-seeds', '20',
                    '-p', thread,
                    '-x', index,
                    '--known-splicesite-infile', '{}/sp.txt'.format(hisat_dir),
                    '-1', m1,
                    '-2', m2,
                    '-S', '{}/align.sam'.format(hisat_dir)]
    else:
        hisat_map = ['hisat2',
                    '--no-softclip',
                    '--score-min', 'L,-16,0',
                    '--mp', '7,7',
                    '--rfg', '0,7',
                    '--rdg', '0,7',
                    '--dta',
                    '-k', '1',
                    '--max-seeds', '20',
                    '-p', thread,
                    '-x', index,
                    '-1', m1,
                    '-2', m2,
                    '-S', '{}/align.sam'.format(hisat_dir)]

    with open('{}/hisat_align.log'.format(hisat_dir), 'wb') as f:
        print('# start to align to genome by hisat')
        f.write(subprocess.check_output(hisat_map,
                                        stderr=subprocess.STDOUT))

    # get align.bam and unmapped.fq 
    samfile = pysam.AlignmentFile('{}/align.sam'.format(hisat_dir))
    bamfile = pysam.AlignmentFile('{}/tmp.bam'.format(hisat_dir),
                                  "wb", template=samfile)
    print('# get mapped and unmapped reads')
    with open('{}/unmapped.fq'.format(hisat_dir), 'w') as f:
        for read in samfile.fetch():
            if read.is_unmapped:
                f.write("@{}\n{}\n+\n{}\n".format(read.query_name,
                                                read.query_sequence,
                                                "".join([chr(x + 33) for x in read.query_qualities])
                                                )
                        )
            else:
                bamfile.write(read)

    bamfile.close()
    samfile.close()

    print('# sort bam file')
    pysam.sort("-o",
               '{}/align.bam'.format(hisat_dir),
               '{}/tmp.bam'.format(hisat_dir))
    print('# index bam file')
    pysam.index('{}/align.bam'.format(hisat_dir) )

def fusion_align(fq, index, fusion_dir, thread):
    if sys.version_info[0] > 2:
        print('TopHat2 is not compable to python3.')
        print('You can change to a python2 environment or \
change tophat shebang "#!/usr/bin/env python" to "#!/usr/bin/env python2"')

    if not which('tophat2'):
        sys.exit('Please install tophat2 at first.')

    fusion_map = ['tophat2',
                  '-o', fusion_dir,
                  '-p', thread,
                  '--fusion-search',
                  '--keep-fasta-order',
                  '--bowtie1',
                  '--no-coverage-search',
                  index,
                  fq]
    subprocess.check_output(fusion_map,
                            stderr=subprocess.STDOUT)


def circ_annot(bam, genome_fa, gtf, circ_dir):
    # setup work environment
    if not which('CIRCexplorer2'):
        sys.exit('Please install CIRCexplorer2 at first.')
    if not which('gtfToGenePred'):
        sys.exit('Please install gtfToGenePred at first.')

    gtf_cmd = ['gtfToGenePred',
               '-genePredExt',
               '-allErrors',
               gtf,
               '{}/genePred.tmp'.format(circ_dir)]
    try:
        from subprocess import DEVNULL # py3k
    except ImportError:
        DEVNULL = open(os.devnull, 'wb')
    subprocess.call(gtf_cmd, stdout=DEVNULL, stderr=subprocess.STDOUT)

    with open('{}/genePred.tmp'.format(circ_dir), 'r') as inf,\
            open('{}/annotation.txt'.format(circ_dir), 'w') as out:
        for line in inf:
            line = line.strip().split('\t')
            outline = '\t'.join([line[11]]+line[0:10]) + '\n'
            out.write(outline)
    # os.remove('{}/genePred.tmp'.format(circ_dir))

    circ_parse = ['CIRCexplorer2',
                  'parse',
                  '-f',
                  '-t', 'TopHat-Fusion',
                  bam,
                  '-b', '{}/bsj.bed'.format(circ_dir)]
    subprocess.check_output(circ_parse)

    circ_annotate = ['CIRCexplorer2',
                     'annotate',
                     '-r', '{}/annotation.txt'.format(circ_dir),
                     '-g', genome_fa,
                     '-b', '{}/bsj.bed'.format(circ_dir),
                     '-o', '{}/circular.txt'.format(circ_dir)]
    subprocess.check_output(circ_annotate)


def main():
    parser = argparse.ArgumentParser()

    # input reads
    parser.add_argument('-1', dest='m1',
                        required=True,
                        help='Comma-separated list of read sequence files in FASTQ format. When running with pair-end read, this should contain #1 mates.')
    parser.add_argument('-2', dest='m2',
                        required=False,
                        help='Comma-separated list of read sequence files in FASTQ format. -2 is only used when running with pair-end read. This should contain #2 mates.')

    # reference and index
    parser.add_argument('-g', '--genome', dest='genome',
                        required=True,
                        help='Genome FASTA file')
    parser.add_argument('-i', '--hisat', dest='hisat',
                        required=True,
                        help='Index files for HISAT2')
    parser.add_argument('-j', '--bowtie1', dest='bowtie1',
                        required=True,
                        help='Index files for TopHat-Fusion')

    # gene annotation
    parser.add_argument('-G', '--gtf', dest='gtf',
                        required=True,
                        help='Annotation GTF file.')

    # output
    parser.add_argument('-o', '--output', dest='output',
                        default='clear_output',
                        help='The output directory')

    # others
    parser.add_argument('-p', '--thread', dest='thread',
                        default='5',
                        help='Running threads. [default: 5]')

    args = parser.parse_args()
    print('###Parameters:')
    print(args)
    print('###Parameters')

    # for setup work environment
    hisat_dir = '{}/{}'.format(args.output, 'hisat')
    fusion_dir = '{}/{}'.format(args.output, 'fusion')
    circ_dir = '{}/{}'.format(args.output, 'circ')
    quant_dir = '{}/{}'.format(args.output, 'quant')
    # test: comment out below four commands
    os.mkdir(args.output)

    print('\n###Start hisat2 mapping')
    os.mkdir(hisat_dir)
    hisat_align(args.m1, args.m2, args.hisat, args.gtf, hisat_dir, args.thread)
    print('###End hisat2 mapping\n')

    print('\n###Start tophat-fusion mapping')
    fusion_align('{}/unmapped.fq'.format(hisat_dir),
                 args.bowtie1, fusion_dir, args.thread)
    print('###End tophat-fusion mapping\n')

    print('\n###Start circRNA annotation')
    os.mkdir(circ_dir)
    circ_annot('{}/accepted_hits.bam'.format(fusion_dir),
                 args.genome, args.gtf, circ_dir)
    print('###End circRNA annotation\n')

    print('\n###Start circRNA quantification')
    os.mkdir(quant_dir)
    tmp_argv = sys.argv
    sys.argv = ['circ_quant.py',
                '-c', '{}/circular.txt'.format(circ_dir),
                '-b', '{}/align.bam'.format(hisat_dir),
                '-t',
                '-r', '{}/annotation.txt'.format(circ_dir),
                '-o', '{}/quant.txt'.format(quant_dir)]
    circ_quant.main()
    print('###End circRNA quantification\n')

if __name__ == '__main__':
    main()
