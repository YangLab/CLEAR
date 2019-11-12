# CLEAR (CIRCexplorer 3)
A computational pipeline for **C**ircular and **L**inear RNA **E**xpression **A**nalysis from **R**ibosomal-RNA depleted (**R**ibo–) RNA-seq (CLEAR; CIRCexplorer 3)

## Schema
![pipeline](/docs/pipeline.png)

## Installation requirements
* Software
    - [CIRCexplorer2](https://github.com/YangLab/CIRCexplorer2/tree/master/circ2) (>=2.3.6)
    - [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) (>=2.0.5)
    - [StringTie](https://ccb.jhu.edu/software/stringtie) (>1.3.6)
* Package (python 2.7 +)
    - [pysam](http://pysam.readthedocs.org/en/latest/) (>=0.8.4)
    - [pybedtools](http://daler.github.io/pybedtools/)

## Installation
```bash
git clone https://github.com/YangLab/CLEAR
cd CLEAR
python ./setup.py install
```

## Usage
Start from fastq file:
```
usage: clear_quant [-h] -1 M1 [-2 M2] -g GENOME -i HISAT -j BOWTIE1 -G GTF
              [-o OUTPUT] [-p THREAD]

optional arguments:
  -h, --help            Show this help message and exit.
  -1 M1                 Comma-separated list of read sequence files in FASTQ
                        format. When running with pair-end read, this should
                        contain #1 mates.
  -2 M2                 Comma-separated list of read sequence files in FASTQ
                        format. -2 is only used when running with pair-end
                        read. This should contain #2 mates.
  -g GENOME, --genome GENOME
                        Genome FASTA file.
  -i HISAT, --hisat HISAT
                        Index files for HISAT2.
  -j BOWTIE1, --bowtie1 BOWTIE1
                        Index files for TopHat-Fusion.
  -G GTF, --gtf GTF     Annotation GTF file.
  -o OUTPUT, --output OUTPUT
                        The output directory.
  -p THREAD, --thread THREAD
                        Running threads. [default: 5]
```
Start from CIRCexplorer2 output file:
```
usage: circ_quant [-h] -c CIRC -b BAM -r REF [--threshold THRESHOLD]
                     [--ratio RATIO] [-l] [-t] [-o OUTPUT]

optional arguments:
  -h, --help            Show this help message and exit.
  -c CIRC, --circ CIRC  Input circular RNA file from CIRCexplorer2.
  -b BAM, --bam BAM     Input mapped reads from HISAT2 in BAM format.
  -r REF, --ref REF     The refFlat format gene annotation file.
  --threshold THRESHOLD
                        Threshold of FPB for choose circRNAs to filter linear
                        SJ.[default: 1]
  --ratio RATIO         The ratio is used for adjust comparison between circ
                        and linear.[default: 1]
  -l, --length          Whether to consider all reads' length? [default: False]
  -t, --tmp             Keep tmp dir? [default: False]
  -o OUTPUT, --output OUTPUT
                        Output file. [default: circRNA_quant.txt]
```

### Example
Start from fastq file:
```bash
clear_quant -1 mate_1.fastq -2 mate_2.fastq -g hg38.fa -i hg38.hisat_index -j hg38.bowtie_index -G annotation.gtf -o output_dir
```
Start from CIRCexplorer2 output file:
```bash
circ_quant -c CIRCexplorer2_output.txt -b hisat_aligned.bam -t -r annotation.refFlat -o quant.txt
```

*hisat_aligned.bam should not contain unmapped reads.*

### Output
* output_dir/quant/quant.txt

| Field       | Description                           |
| :---------- | :------------------------------------ |
| chrom       | Chromosome                            |
| start       | Start of circular RNA                 |
| end         | End of circular RNA                   |
| name        | Circular RNA/Junction reads           |
| score       | Flag of fusion junction realignment   |
| strand      | + or - for strand                     |
| thickStart  | No meaning                            |
| thickEnd    | No meaning                            |
| itemRgb     | 0,0,0                                 |
| exonCount   | Number of exons                       |
| exonSizes   | Exon sizes                            |
| exonOffsets | Exon offsets                          |
| readNumber  | Number of junction reads              |
| circType    | Type of circular RNA                  |
| geneName    | Name of gene                          |
| isoformName | Name of isoform                       |
| index       | Index of exon or intron               |
| flankIntron | Left intron/Right intron              |
| FPBcirc     | Expression of circRNA                 |
| FPBlinear   | Expression of cognate linear RNA      |
| CIRCscore   | Relative expression of circRNA        |


## Citation
**Ma XK\*, Wang MR, Liu CX, Dong R, Carmichael GG, Chen LL and Yang L#. A CLEAR pipeline for direct comparison of circular and linear RNA expression. 2019, bioRxiv doi: 10.1101/668657**


## License
Copyright (C) 2019 YangLab. Licensed GPLv3 for open source use or contact YangLab (yanglab@@picb.ac.cn) for commercial use.
