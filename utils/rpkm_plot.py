#!/usr/bin/env python 

import sys 
import numpy
import pandas 
import matplotlib.pyplot as plt
from collections import defaultdict 

from utils import GFFParser

try:
    rquant_file = sys.argv[1]
    fout_name = sys.argv[2]
except:
    exit("error")

rpkm = GFFParser.Parse(rquant_file)

score_freq = defaultdict(int) 
for gene in rpkm:
    for idx, trans in enumerate(gene["transcripts"]):
        val = gene["transcript_score"][idx]
        val = float("%.1f" % round(float(val), 1))
        score_freq[val] += 1 
          
score_freq_tab = numpy.zeros((len(score_freq), 2))
for idx, score in enumerate(sorted(score_freq.keys())):
    score_freq_tab[idx] = numpy.array([score, score_freq[score]])

df_score_freq = pandas.DataFrame(score_freq_tab) 
#import ipdb 
#ipdb.set_trace() 

## plotting settings  
width = 0.1
fig = plt.figure()

ind = numpy.arange(len(score_freq))
plt.bar(ind, numpy.log10(df_score_freq[1]), color="#A9CCE3", edgecolor='#A9CCE3')
    
max_index = len(df_score_freq[0])-1
#xlocations = [df_score_freq[0][int(max_index*i)] for i in [0.25, 0.5, 0.75, 1]] 

plt.yticks(fontsize=9)
plt.xticks(fontsize=9) 
plt.xlabel('rpkm values range: %d - %d' % (df_score_freq[0][0], df_score_freq[0][max_index]), fontsize=9)
plt.ylabel("frequency (log10)", fontsize=9)

plt.savefig(fout_name)

