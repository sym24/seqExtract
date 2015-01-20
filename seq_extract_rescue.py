#!/usr/bin/env python
"""
Required files:
hgnc_table.tsv in directory contains this script
"""
import os
import sys
from pybedtools import BedTool
import pysam
import subprocess
import argparse
import csv
import rescue

def get_gene_name(target_file):
    """Generate set contains gene name without duplication."""
    with open(target_file, 'r') as fi:
        gene_set = set(line.rstrip() for line in fi)
    print 'Total number of genes: ', len(gene_set)
    return (gene_set, len(gene_set))

def get_bounds(gene_set, buffer, gtf_file):
    """Go through gtf file extract the longest sequence for genes in gene_set"""
    print '...Extract gene boundaries...'
    bounds = {}
    for feature in BedTool(gtf_file):
        if feature[2] == 'exon' and feature.attrs.has_key('gene_id'):
             gene = feature.attrs['gene_id']
             if gene in gene_set:
                 startp = int(feature.start) + 1 - buffer
                 endp = int(feature.stop) + buffer
                 try:
                     # update the starting and end point
                     if startp < bounds[gene][1]:
                         bounds[gene][1] = startp
                     if endp > bounds[gene][2]:
                         bounds[gene][2] = endp
                 except KeyError:
                     # add new key
                     bounds[gene] = [feature[0], startp, endp]

    find_num = 0
    miss_genes = []
    for gene in gene_set:
        if gene in bounds:
            find_num += 1
        else:
            miss_genes.append(gene)
    print 'Number of genes successfully captured: %d'%(find_num)
    if miss_genes:                
        print 'Number of genes failed to captured: %d'%len(miss_genes)
        print 'Missing genes: %s'%', '.join(miss_genes)
    return bounds, miss_genes 

def extract_seq(target_fa, ref_file, bounds):
    """Extract genes from reference genome and write it to output file."""
    print '...Extract gene sequence...'
    with open(target_fa, 'a') as fi:    
        ref_fasta = pysam.Fastafile(ref_file)
        for gene, coord in bounds.iteritems():
            seq = ref_fasta.fetch(coord[0], coord[1], coord[2])
            fi.write("%s%s %s:%d-%d\n%s\n" % ('>',gene,coord[0],coord[1],coord[2],seq))

def makeBloomFilter(input_file, outdir, prefix):
    """Make bloom filter out of targeted genes"""
    target_bf = os.path.join(outdir, prefix+'.bf')
    try:
        print '>>>Run samtools'
        subprocess.check_output(['samtools', 'faidx', input_file])
        cmd_pass = ['biobloommaker', '-p', prefix, '-o', outdir, input_file]
        print '>>>Run biobloommaker'
        subprocess.check_output(cmd_pass)
    except subprocess.CalledProcessError as e:
        sys.exit(e)    
    return target_bf

def main():
    
    parser = argparse.ArgumentParser(description='Extract target genes and '
                                     'generate bloom filter')
    # Optional variables
    parser.add_argument('-o', '--outdir', dest='outdir', metavar='STR', type=str,
                        help='Output file path. Default: [%(default)s]',
                        default=os.getcwd())
    parser.add_argument('-p', '--prefix', dest='prefix', metavar='STR', type=str, 
                        help='Output file prefix. Default: [%(default)s]',
                        default='target_genes')
    parser.add_argument('-b', '--buffer', dest='buffer', metavar='INT', type=int,
                        help='Size of extended buffer on both end of '
                        'extracted sequence. Defaule: [%(default)s]', default='0')
    parser.add_argument('-f', '--force', dest='overwrite', action='store_true',
                        help='Overwrite existing output files')
    parser.add_argument('-c', '--cosmic', dest='cosmic', action = 'store_true',
                        help='Using cosmic gene file containing '
                        'chromsome information as input gene list')

    # Required arguments
    
    parser.add_argument('-g', '--genes', dest='genes', metavar='file',
                        type=argparse.FileType('r'), required=True,
                        help='Absolute path of text file containing '
                        'target gene names. Required.')
    parser.add_argument('-t', '--gtf', dest='gtf', metavar='gtf_file',
                        type=file, required=True,
                        help='Absolute path of gtf file. Required')
    parser.add_argument('-r', '--ref', dest='ref', metavar='ref_genome',
                        type=file, required=True,
                        help='Absolute path of the reference sequence file. Requried.')

    args = parser.parse_args()

    # Make output file name
    target_fa = os.path.join(args.outdir, args.prefix+'.fa')

    # Make output directory if not exist
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    elif os.path.isfile(target_fa):
        if args.overwrite:
            print 'File {} exists, will overwrite'.format(target_fa)
        else:
            sys.exit('File {} exist\nPlease use "-f" to overwrite or '
                     'change outdir or prefix'.format(target_fa))

    print 'Genome file: ', os.path.abspath(args.ref.name)
    print 'Gene list: ', os.path.abspath(args.genes.name)
    print 'GTF file: ', os.path.abspath(args.gtf.name)
    print 'Output FASTA file: ', target_fa
    
    if(args.cosmic):
        print 'Input genes list from cosmic gene info sheet.'
        # Extract gene name from input file 
        gene_set, chrom = rescue.get_cosmic_info(args.genes.name)
        # Get gene coordinates from GTF file  
        bounds, miss_genes = get_bounds(gene_set, args.buffer, args.gtf.name)
        extract_seq(target_fa, args.ref.name, bounds)

        print 'Rescure genes failed to extract'
        new_set, still_miss_set, newdict, replace_num = \
                 rescue.rescue_miss_gene(miss_genes, 'hgnc_table.tsv')
        new_bounds, failed_rescue = rescue.check_rescue_bounds\
                    (miss_genes, chrom, newdict, new_set, args.buffer, args.gtf.name)

        print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Results:'
        extracted = len(bounds) + len(new_bounds)
        print 'Initially number of cosmic gene extracted C ==> ', len(bounds) 
        print 'Initially number of cosmic gene unable to extract: ', len(miss_genes)
        print 'Cosmic genes unable to be replaced with new names:\n%s'%still_miss_set
        print '==>Number A ', len(still_miss_set)
        print 'Number of cosmic genes successfully get new name: ', replace_num
        print ('Cosmic genes with new names unable to extract sequence from gtf: '
               '\n%s \nnumber B ==> %d'%(failed_rescue, len(failed_rescue)))
        print 'Number of cosmic genes successfully rescued with new names D ==> %d/%d'%\
              (len(new_bounds), len(miss_genes))
        print 'Totall gene number : %d'%len(gene_set)

        rescue.extract_rescue_seq(target_fa, args.ref.name, new_bounds)
    else:
        # Extract gene name from input file 
        gene_set, gene_num = get_gene_name(args.genes.name)
        # Get gene coordinates from GTF file  
        bounds, miss_genes = get_bounds(gene_set, args.buffer, args.gtf.name)
        extract_seq(target_fa, args.ref.name, bounds)
        
    # Generate bloom filter
    target_bf = makeBloomFilter(target_fa, args.outdir, args.prefix)
    
if __name__ == '__main__':
    main()
