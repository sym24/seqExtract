#!/usr/bin/env python

import os
import sys
from pybedtools import BedTool
import pysam
import csv

def get_cosmic_info(target_file):
    """Extract cosmic genes approved symbol into a set and
    a dictionary contains chromsome location
    """
    chrom = {}
    gene_set = set()
    with open(target_file, 'r') as csvfile:
        incsv = csv.reader(csvfile)
        hasheader = csv.Sniffer().has_header(csvfile.read(1024))
        csvfile.seek(0)
        if hasheader:
            next(incsv)
        for line in incsv:
            gene_set.add(line[0])
            if line[0] not in chrom:
                chrom[line[0]] = line[3]
        if 'Gene Symbol' in gene_set:
            print 'Failed to skip first row of csv file'
    print 'Total number of genes: ', len(gene_set)
    return gene_set, chrom

def rescue_miss_gene(miss_genes, hugo_table):
    """Rescue missing genes using hgnc_table
    by replacing old genes with new gene name            
    """
    print '...Rescue missing genes...'
    newdict = {}
    # Store old gene names as dic key, new gene names and 
    with open(hugo_table, 'r') as fi1:
        for line in csv.reader(fi1, dialect='excel-tab'):
            old = line[2].split(', ') + line[3].split(', ')
            for ele in old:
                if ele in newdict:
                    newdict[ele].append(line[0])
                else:
                    newdict[ele] = [line[0]]
    new_list = []
    still_miss_list = []
    new_gene_num = 0
    for oldgene in miss_genes:
        if oldgene in newdict:
            new_gene_num += 1
            for ele in newdict[oldgene]:
                new_list.append(ele)
        else:
            still_miss_list.append(oldgene)
    still_miss_set = set(still_miss_list)
    new_set = set(new_list)
    print 'New genes set: ', ', '.join(new_set)
    print 'Still missing old genes: ', ', '.join(still_miss_set)
    print '%d old genes found new names: '%new_gene_num
    return (new_set, still_miss_set, newdict, new_gene_num)

def check_rescue_bounds(miss_gene, chrom, newdict, new_set, buffer, gtf_file):
    """
    Extract rescued genes' bounds and check chromsome match to initial cosmic genes
    """
    print '...Get rescued gene boundaries...'
    bounds = {}
    for feature in BedTool(gtf_file):
        if feature[2] == 'exon' and feature.attrs.has_key('gene_id'):
             gene = feature.attrs['gene_id']
             if gene in new_set:
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

    failed_rescue = []
    new2old = {}
    for oldgene in miss_gene:
        if oldgene in newdict:
            print '{oldgene: newgenes}  {%s: %s}'%(oldgene, newdict[oldgene])
            flag1 = flag2 = 0
            for newgene in newdict[oldgene]:
                if newgene in bounds:
                    coord = bounds[newgene]
                    if coord[0][3:] == chrom[oldgene]:
                        print '...chromsome matches...old gene %s rescued: '%oldgene
                        print '...new gene %s is extracted t'%newgene
                        # if chromsome match, add old gene to new gene bounds
                        bounds[newgene].append(oldgene)
                        if newgene not in new2old:
                            new2old[newgene]=oldgene
                        else:
                            print '......!!!New gene %s already in bounds!!!...... '%newgene
                            print 'New gene %s ==> Existing old gene %s'%\
                            (newgene, new2old[newgene])
                    else:
                        flag2 += 1
                        print 'Sequence chromsome does not match. Delete key', newgene
                        del bounds[newgene]
                        if flag2 == len(newdict[oldgene]):
                            failed_rescue.append(oldgene)
                else:
                    flag1 += 1
                    print 'Unable to extract bounds for ', newgene
                    if flag1 == len(newdict[oldgene]):
                        failed_rescue.append(oldgene)
        else:
            print 'Old gene %s is not found in new gene dictionary'%oldgene

    return (bounds, failed_rescue)

def extract_rescue_seq(target_fa, ref_file, bounds):
    """
    Extract genes from reference genome and write it to output file.
    """
    print '...Extract rescued gene sequence...'
    with open(target_fa, 'a') as fi:    
        ref_fasta = pysam.Fastafile(ref_file)
        for gene, coord in bounds.iteritems():
            seq = ref_fasta.fetch(coord[0], coord[1], coord[2])
            fi.write("%s%s(%s) %s:%d-%d\n%s\n"\
            %('>',gene,', '.join(coord[3:]),coord[0],coord[1],coord[2],seq))

if __name__ == '__main__':
    print 'rescue funtions for sequence extraction'
