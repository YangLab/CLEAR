# CLEAR
A computational pipeline for circular and linear RNA expression analysis from ribosomal-RNA depleted (ribo–) RNA-seq (CLEAR

## Schema
![pipeline](/docs/pipeline.png)

## Installation requirements
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
  -h, --help            show this help message and exit
  -1 M1                 Comma-separated list of read sequence files in FASTQ
                        format. When running with pair-end read, this should
                        contain #1 mates.
  -2 M2                 Comma-separated list of read sequence files in FASTQ
                        format. -2 is only used when running with pair-end
                        read. This should contain #2 mates.
  -g GENOME, --genome GENOME
                        Genome FASTA file
  -i HISAT, --hisat HISAT
                        Index files for HISAT2
  -j BOWTIE1, --bowtie1 BOWTIE1
                        Index files for TopHat-Fusion
  -G GTF, --gtf GTF     Annotation GTF file.
  -o OUTPUT, --output OUTPUT
                        The output directory
  -p THREAD, --thread THREAD
                        Running threads. [default: 5]
```
Start from CIRCexplorer2 output file:
```
usage: circ_quant [-h] -c CIRC -b BAM -r REF [--threshold THRESHOLD]
                     [--ratio RATIO] [-l] [-t] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -c CIRC, --circ CIRC  input circular RNA file
  -b BAM, --bam BAM     to get mapped reads
  -r REF, --ref REF     refFlat format gene annotation file.
  --threshold THRESHOLD
                        threshold of HPB for choose circRNAs to filter linear
                        sp.[default: 1]
  --ratio RATIO         the ratio is used for adjust comparison between circ
                        and linear.[default: 93.0/85]
  -l, --length          Wether to consider all reads' length? [default: False]
  -t, --tmp             Keep tmp dir? [default: False]
  -o OUTPUT, --output OUTPUT
                        output file. [default: circRNA_quant.txt]
```

### Example
Start from fastq file:
```bash
clear_quant -1 mate_1.fastq -2 mate_2.fastq -g hg38.fa -i hg38.hisat_index -j hg38.bowtie_index -G annotation.gtf -o output_dir
```
Start from fastq file:
```bash
circ_quant.py -c CIRCexplorer2_output.txt -b hisat_aligned.bam -t -r annotation.refFlat -o quant.txt
```

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
**Ma XK\*, Liu CX, Wang MR, Dong R, Chen LL# and Yang L#. A CLEAR pipeline to directly compare circular and linear RNA expression. 2019 (Submitted)**


## License
Copyright (C) 2019 YangLab.
Licensed GPLv3 for open source use
or contact YangLab for commercial use
