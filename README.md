# seqExtract
Extract sequence from reference geonome by giving gene names
rescue.py and seq_extract_rescue.py are together to rescue missing genes from cosmic file.
seq_extract.py only need file names as input.

```
~$ ./seq_extract.py --help
usage: seq_extract.py [-h] -g file -t gtf_file -r ref_genome [-o STR] [-p STR]
                      [-b INT] [-f]

Extract target genes and generate bloom filter

optional arguments:
  -h, --help            show this help message and exit
  -g file, --genes file
                        Absolute path of text file containing target gene
                        names. [Required]
  -t gtf_file, --gtf gtf_file
                        Absolute path of gtf file.[Required]
  -r ref_genome, --ref ref_genome
                        Absolute path of the reference sequence file.
                        [Requried]
  -o STR, --outdir STR  Output file path. Default:
                        [/home/ymingsun/work/seqExtract]
  -p STR, --prefix STR  Output file prefix. Default: [target_genes]
  -b INT, --buffer INT  Size of extended buffer on both end of extracted
                        sequence. Defaule: [0]
  -f, --force           Overwrite existing output files

```
