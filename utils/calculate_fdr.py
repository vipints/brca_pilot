"""

- helper functions for estimating False Disocvery Rates from Pvalues 

"""
import scipy as SP

from collections import defaultdict 

from utils import GFFParser


def gene_id_name_map(gtf_file):
    """
    read the gtf file to get a mapping between gene id and gene name 
    """

    gtf_db = GFFParser.Parse(gtf_file) 

    gene_id2name = defaultdict() 
    gene_name2id = defaultdict() 
    for gene in gtf_db:
        gene_id2name[gene["gene_info"]["ID"]]=gene["gene_info"]["Name"]
        gene_name2id[gene["gene_info"]["Name"]]=gene["gene_info"]["ID"]

    return gene_id2name, gene_name2id


def get_rdiff_p_value(rdiff_res):
    """
    parse the rdiff outfile 

    @args rdiff_res: rdiff result file genes with p-values  
    @type rdiff_res: str 
    """

    rfh = open(rdiff_res, "rU") 
    genes_p_val = defaultdict(list) 
    for rec in rfh: 
        rec = rec.strip("\n\r").split("\t")
        if rec[0] == "gene": ## remove header 
            continue 
        if float(rec[1]) < 1.0: ## consider p-value less than 1.0 to calculate FDR 
            genes_p_val[rec[0]].append(float(rec[1])) 
    rfh.close() 
    print("considering %d genes for FDR calculation" % len(genes_p_val))

    ## FIXME merge genome annotation to fix this problem
    """
    genes_corr_p_val = defaultdict() 
    for gid, p_val in genes_p_val.items():
        if len(p_val) > 1: 
            p_val.sort() 
            genes_corr_p_val[gid] = p_val[-1] 
            continue 
        genes_corr_p_val[gid] = p_val[0]
    print len(genes_corr_p_val) 
    """

    p_val_to_genes = defaultdict(list) 
    for gid, val in genes_p_val.items():
        p_val_to_genes[val].append(gid) 

    PVAL = SP.zeros_like(p_val_to_genes.keys())
    for idx, pval in enumerate(p_val_to_genes.keys()):
        PVAL[idx] = pval 

    return PVAL, p_val_to_genes 


def q_value_to_gene(q_val, p_val, pval_genes, gene_id_name_map):
    """
    adjusted p value to the genes
    
    """
    candidates = dict(RAD51=0,
    PALB2=0,
    FANCD2=0,
    BRCA1=0,
    SHFM1=0,
    FANCG=0,
    CDK2=0,
    ATM=0,
    RAD51C=0,
    XRCC3=0,
    DMC1=0,
    TP53=0,
    FANCI=0,
    MRE11A=0,
    RAD50=0,
    HMG20B=0,
    MLH1=0,
    C11orf30=0,
    RPA1=0,
    H2AFX=0)

    for idx, qv in enumerate(q_val):
        corrsp_p_val = p_val[idx] 

        for gene in pval_genes[corrsp_p_val]:
            gname = gene_id_name_map[gene] 
            
            if gname in candidates:
                print gene, gname, corrsp_p_val, qv



def estimate_q_values(PV,m=None,pi=1):
    """
    estimate q vlaues from a list of Pvalues
    this algorihm is taken from Storey, significance testing for genomic ...
    m: number of tests, (if not len(PV)), pi: fraction of expected true null (1 is a conservative estimate)

    @author: Oliver Stegle
    """

    if m is None:
        m = len(PV)
    lPV = len(PV)
    
    #1. sort pvalues
    PV = PV.squeeze()
    IPV = PV.argsort()
    PV  = PV[IPV]

    #2. estimate lambda
    if pi is None:
        lrange = SP.linspace(0.05,0.95,max(lPV/100,10))
        pil    = SP.double((PV[:,SP.newaxis]>lrange).sum(axis=0))/lPV
        pilr   = pil/(1-lrange)
        #ok, I think for SNPs this is pretty useless, pi is close to 1!
        pi =1
        #if there is something useful in there use the something close to 1
        if pilr[-1]<1:
            pi = pilr[-1]
            
    #3. initialise q values
    QV_ = pi * m/lPV* PV
    #4. update estimate
    for i in xrange(lPV-2,0,-1):
        QV_[i] = min(pi*m*PV[i]/(i+1),QV_[i+1])
    #5. inverst sorting
    QV = SP.zeros_like(PV)
    QV[IPV] = QV_
    return QV


if __name__ == '__main__':
    a=SP.rand(5)
    print "p-value ", a 
    c=estimate_q_values(a)
    print "q-value ", c 
    
