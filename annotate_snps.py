__author__ = 'aydar'

import vcf
import pandas as pd
import matplotlib.pyplot as plt
import os
import peakdetect as pk
import scipy.stats.mstats as mst
import numpy as np


# parse stats
def parse_stats(directory_with_stat, head_presence):
    chromosomes = {}
    directory_with_stat_list = os.listdir(directory_with_stat)
    for i in directory_with_stat_list:
        chrom_num = i.split(".")[0]
        chrom = pd.read_table("%s/%s" % (directory_with_stat, i), header=head_presence)
        if head_presence == None:
            chrom.columns = ["rs", "pos", "c", "d", "e", "f", "stat", "giovanni"]
            chrom.giovanni = -np.log10(mst.zscore(chrom.stat.abs()))
        elif head_presence == 0:
            chrom.rename(columns = {"xpehh":"stat"}, inplace = True)
        chromosomes[chrom_num] = chrom

    return chromosomes

# function to extract gene annotation
def parse_vcf_annot(vcf_annot):
    genes_annot = []
    position = []
    snps_names = []
    vcf_reader = vcf.Reader(open(vcf_annot, 'r'))
    for record in vcf_reader:
        gene = record.INFO['ANN'][0].split('|')[3]
        pos = record.POS
        rs = record.ID
        position.append(pos)
        genes_annot.append(gene)
        snps_names.append(rs)
    parsed_df = pd.DataFrame({'RS': snps_names, "POS": position, 'GENE': genes_annot})
    return parsed_df


#calculate heterozygosity
def calculate_heterozygosity(vcfile):
    heterozygosity_results = pd.DataFrame(columns = ('rs', 'pos', 'chr', 'het_exp', 'het_obs'))
    vcf_reader = vcf.Reader(open(vcfile, 'r'))
    i = 0
    for record in vcf_reader:
        i += 1
        if i%10000 == 0:
            print(i)
        N = len(record.samples) 
        het_obs = record.num_het/N
        p = record.num_hom_ref/N + 1/2*het_obs
        q = record.num_hom_alt/N + 1/2*het_obs
        het_exp = 2*p*q
        heterozygosity_results.loc[heterozygosity_results.shape[0]] = [record.ID, record.POS, record.CHROM, het_exp, het_obs]
    return heterozygosity_results


# get peaks vertices using peakdetect, maybe write another function for xpehh, because in that case negative values matter
def detect_peaks(stats_file):
    peaks = pk.peakdetect(stats_file.stat.abs(), stats_file.pos, lookahead = 200)
    peaks = peaks[0]
    return peaks

# annotating function, should be called after plotting (scatter in this case)
def annotate_peaks(stats_table, peaks_pos, genes_to_annot, axis, cvet):
    peaks_coord = [i[0] for i in peaks_pos]
    print(len(peaks_coord))
    genes_to_annot.index = genes_to_annot.RS
    for pos, stat, rs in zip(stats_table.pos, stats_table.stat, stats_table.rs):
        if pos in peaks_coord and rs in genes_to_annot.RS.tolist():
            print('yay')
            gene = genes_to_annot.loc[rs].GENE
            axis.annotate(gene, xy=(pos, stat), xytext=(pos+0.3, stat+0.3), color = cvet, fontsize = 9,
                        rotation=-75, horizontalalignment='right', verticalalignment='top',
                        arrowprops=dict(arrowstyle="->", connectionstyle="arc, angleA=0, armA=0, rad=0"))


#plotting stats
def plot_stats(chrom, color_p, label_legend, **kwargs):
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    leg = plt.scatter(chrom.pos[chrom.stat.abs() > 2], chrom.stat.abs()[chrom.stat.abs() > 2], s=10, c=color_p, alpha=0.5, label = label_legend)
    annotate_peaks(chrom, axis = ax, cvet = color_p, **kwargs)
    plt.autoscale(enable=True, axis='both', tight=True)
    return leg

#plot heterozygocity

def plot_heterozygosity_2_pop(het1, het2, label_leg_name_1, label_leg_name_2):
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    leg1 = plt.scatter(het1.pos, het2.het_exp, s=10, c = 'black', alpha = 0.5, label = label_leg_name_1)
    leg2 = plt.scatter(het2.pos, -het2.het_exp, s = 10, c = 'red', alpha = 0.5, label = label_leg_name_2)



#annotate_axes
def annotate_axes(pop1, pop2, chrom, a, b, stat):
    plt.title(pop1+" vs "+pop2+" chrom %s"%chrom)
    plt.legend(handles=[a,b], loc = 'center left', bbox_to_anchor = (1,0.5))
    plt.ylabel('|%s|'%stat)
    plt.xlabel('Mb')

"""for i in directory_ki:
    one = pd.read_table("./KI/%s" % i, header=None)
    two = pd.read_table("./KO/%s" % i, header=None)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.scatter(one[1][one[6].abs() > 2], one[6].abs()[one[6].abs() > 2], s=10, c="black", alpha=0.5)
    plt.scatter(two[1][two[6].abs() > 2], two[6].abs()[two[6].abs() > 2], s=10, c="red", alpha=0.5)

    plt.title("KI and KO populations, Chrom " + i, fontsize=12, horizontalalignment="center")

    plt.autoscale(enable=True, axis='both', tight=True)

    plt.xlabel('Mb')
    plt.ylabel('|iHS|')
    # plt.legend((l2, l1), ('oscillatory', 'damped'), loc='upper right', shadow=False)
    plt.savefig(i + '_graph.png')"""
