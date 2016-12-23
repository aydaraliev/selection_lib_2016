#!/home/aydar/anaconda3/bin/python

import sys
import pandas as pd
import os

chr_to_concat = []

for filename in os.listdir("./"):
    if filename.endswith(".norm"):
        chrom_n = filename.split('.')
        chromosome=pd.read_table(filename, header = None)
        chromosome = chromosome[[0,1,6]]
        chrom_n_vector = pd.Series([chrom_n[0] for i in range(len(chromosome[0]))])
        chromosome[7] = chrom_n_vector
        chr_to_concat.append(chromosome)

all_chromosomes = pd.concat(chr_to_concat)
all_chromosomes.index = range(len(all_chromosomes))
all_chromosomes.columns = ['SNP', 'POS', 'NORM_STAT', 'CHR']
all_chromosomes.POS = all_chromosomes.POS.astype('int')
all_chromosomes.CHR = all_chromosomes.CHR.astype('int')
all_chromosomes = all_chromosomes.sort(['CHR', 'POS'], ascending = [1,1])      
all_chromosomes.NORM_STAT = all_chromosomes.NORM_STAT.abs()
all_chromosomes = all_chromosomes[all_chromosomes.NORM_STAT > 3]

       
all_chromosomes.to_csv("%s.csv"%os.path.relpath('.', ".."), index = None)
