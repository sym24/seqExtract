# seqExtract
Extract sequence from reference geonome by giving gene names

rescue.py and seq_extract_rescue.py are together to rescue missing genes from cosmic file.

seq_extract.py only need file names as input.
##seq_extract_rescue.py
####Prerequisite files:
* GTF file(gene transform format)
* Human reference gene
* Cosmic gene list or gene name list

####Prerequisite system:
* Pybedtools
* Pysam
* Python 2.7

####Usage:
```
usage: seq_extract_rescue.py [-h] -g file -t gtf_file -r ref_genome [-o STR]
                             [-p STR] [-b INT] [-f] [-c]

Extract target genes and generate bloom filter

optional arguments:
  -h, --help            show this help message and exit
  -g file, --genes file
                        absolute path of text file containing target gene
                        names. [Required]
  -t gtf_file, --gtf gtf_file
                        absolute path of gtf file. [Required]
  -r ref_genome, --ref ref_genome
                        absolute path of the reference sequence file.
                        [Requried]
  -o STR, --outdir STR  output file path. Default:
                        [/home/ymingsun/work/seqExtract]
  -p STR, --prefix STR  output file prefix. Default: [target_genes]
  -b INT, --buffer INT  size of extended buffer on both end of extracted
                        sequence. Defaule: [0]
  -f, --force           overwrite existing output files
  -c, --cosmic          using cosmic gene file containing chromsome
                        information as input gene list
```
**Note:**  `-c` option specify whether using cosmic gene list or gene name list.

Missing genes are rescued by extracting new names of old genes.
