#!/usr/bin/env python 

import sys 

from utils import GFFParser

try:
    rquant_file = sys.argv[1]
except:
    exit("error")

gtf_file = "gencode.v14.gtf"
genes_gtf = GFFParser.Parse(gtf_file)

gene_name_map = dict() 
for gene in genes_gtf:
    gene_name_map[gene["gene_info"]["ID"]]=gene["gene_info"]["Name"]
    
rpkm = GFFParser.Parse(rquant_file)

for gene in rpkm:

    trans_str = ''
    for idx, trans in enumerate(gene["transcripts"]):
        trans_str += "%s\t%s\t" % (trans[0], gene["transcript_score"][idx])
 
    gid = gene["name"]
    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (gene["chr"], gene["source"], gene["start"], gene["stop"], gene["strand"], gid, gene_name_map[gid], trans_str) 
