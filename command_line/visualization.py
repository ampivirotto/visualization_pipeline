import allel; print('scikit-allel', allel.__version__)
import random, subprocess
random.seed(14)
import time
import numpy as np
np.random.seed(14)
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
import bcolz
#import ipyrad.analysis as ipa
import pandas as pd
import os
from collections import defaultdict
import file_parser as fp
import sys


def randColor(num, pops):
    number_of_colors = num

    colors = []
    r = int(random.random() * 256)
    g = int(random.random() * 256)
    b = int(random.random() * 256)
    step = 256 / number_of_colors
    for i in range(number_of_colors):
        r += step
        g += step
        b += step
        r = int(r) % 256
        g = int(g) % 256
        b = int(b) % 256
        colors.append((r,g,b))

    colordict = {}
    for x in range(len(pops)):
        colordict[pops[x]] = colors[x]
    print(colordict)
    return colordict


def makeVCF(directory, samples, outfn):
    """
    make new vcf file from list of samples
    """
    ## write list of gtc files from list of samples
    gtcfile = open(outfn + ".gtc", "w")
    for sample in samples:
        for file in os.listdir(directory + "/" + sample):
            if file.endswith(".gtc"):
                gtc.write(os.path.join("/mydir", file))
    gtc.close()

    ## run the gtc2vcf conversion

def runConversion(outfn, vcffile, bs):
    """
    convert from vcf to hdf5 format
    """
    #init converter  https://ipyrad.readthedocs.io/en/latest/API-analysis/cookbook-vcf2hdf5.html
    converter = ipa.vcf_to_hdf5(name=outfn, data=vcffile, ld_block_size=bs)

    #run converter
    converter.run()

def plot_ld(gn, title):
    m = allel.rogers_huff_r(gn) ** 2
    ax = allel.plot_pairwise_ld(m)
    ax.set_title(title)

def ld_prune(gn, size, step, threshold=.1, n_iter=1):
    for i in range(n_iter):
        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
        gn = gn.compress(loc_unlinked, axis=0)
    return gn

def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population, populations, pop_colours):

    #import pdb; pdb.set_trace()
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    for i, pop in enumerate(populations):
        #pdb.set_trace()
        flt = (sample_population == pop)
        #pdb.set_trace()
        if pop == "Unknown":
            ax.plot(x[flt], y[flt], marker='o', color = 'black', label=pop, linestyle='None', markersize=6, mec='k', mew=.5)
        else:
            ax.plot(x[flt], y[flt], marker='o', label=pop, linestyle='None', markersize=6, mec='k', mew=.5)
    #pdb.set_trace()
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))

def fig_pca(directory, outfn, coords, model, title, sample_population=None):
    import pylab

    if sample_population is None:
        sample_population = df_samples.population.values
    populations = list(sample_population.unique())
    colordict = randColor(len(populations), populations)

    # plot coords for PCs 1 vs 2, 3 vs 4
    figData = pylab.figure(figsize= (10, 6))
    #fig = pylab.figure()
    #figlegend = pylab.figure(figsize=(3,2))
    #ax = fig.add_subplot(121)
    #lines = ax.plot(range(10), pylab.randn(10), range(10), pylab.randn(10))

    ax = figData.add_subplot(1, 2, 1)
    plot_pca_coords(coords, model, 0, 1, ax, sample_population, populations, colordict)
    ax = figData.add_subplot(1, 2, 2)
    plot_pca_coords(coords, model, 2, 3, ax, sample_population, populations, colordict)

    # create a second figure for the legend
    figLegend = pylab.figure(figsize = (5, 10))

    # produce a legend for the objects in the other figure
    pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper left')

    # save the two figures to files
    #figData.savefig("plot.png")
    figLegend.savefig(directory + "graphics/" + outfn + "_pca_legend.jpg")

    #ax.legend(bbox_to_anchor=(1, 1), loc='center')
    figData.suptitle(title, y=1.02)
    #fig.tight_layout()

    figData.savefig(directory + "graphics/" + outfn + "_pca.jpg")


def transform(g):

    #get allele counts
    ac = g.count_alleles()[:]

    #count mutliallelic and singletons
    print("Multiallelic " + str(np.count_nonzero(ac.max_allele() > 1)))
    print("Singletons " + str(np.count_nonzero((ac.max_allele == 1) & ac.is_singleton(1))))

    #remove singletons and multiallelic
    flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
    gf = g.compress(flt, axis=0)

    ##transform into 2d matrix
    gn = gf.to_n_alt()

    return gn

def prepData(directory, outfn, newVCF, samples, bs):
    ## can either use an existing vcf file or make a new vcf file from a list of samples
    if newVCF == True:
        makeVCF(directory, samples, outfn)

    vcffile = directory + outfn + ".vcf"

    runConversion(outfn, vcffile, bs)

    callsetfn = directory + 'analysis-vcf2hdf5/' + outfn + ".snps.hdf5"
    callset = h5py.File(callsetfn, mode= 'r')

    #callset = allel.read_vcf(vcffile)

    #get genotype data
    g = allel.GenotypeChunkedArray(callset['genos'])

    ## transform data
    gn = transform(g)

    return gn, callset


def LD(directory, outfn, newVCF=False, samples = None, bs = 20000):
    """
    graph the linkage disequilibrium graphs
    """
    gn, callset = prepData(directory, outfn, newVCF, samples, bs)

    #plot ld
    plot_ld(gn[:1000], 'Pairwise LD.')
    plt.show()

    #subsabmple snps
    n = 100000  # number of SNPs to choose randomly
    vidx = np.random.choice(gn.shape[0], n, replace=False)
    vidx.sort()
    gnr = gn.take(vidx, axis=0)

    plot_ld(gnr[:1000], 'Pairwise LD after random downsampling.')
    plt.show()

    #prune linkage
    gnu = ld_prune(gnr, size=500, step=200, threshold=.1, n_iter=5)

    #linkage after pruning
    plot_ld(gnu[:1000], 'Pairwise LD after LD pruning.')
    plt.show()


'''def pca(directory, outfn, column, newVCF=False, samples = None, bs = 20000):
    """
    main function to run pca visualization
    """
    import pdb

    #gn, callset = prepData(directory, outfn, newVCF, samples, bs)
    callset = allel.read_vcf(directory + outfn + ".vcf")

    g = allel.GenotypeChunkedArray(callset['calldata/GT'])
    gn = transform(g)

    ## get metadata
    df = fp.retrieveMetaData(samples, directory, outfn)

    coords1, model1 = allel.pca(gn, n_components=10, scaler='patterson')

    fig_pca(directory, outfn, coords1, model1, 'Conventional PCA.', sample_population = df[column])
'''
def pca(directory, vcffile, outfn):
    gds = vcffile.strip('vcf').strip('.') + '.gds'
    subprocess.run(['Rscript', '--vanilla', 'pca.R', directory, vcffile, outfn + '.html', gds])

def fst(g, directory, outfn, samplelist):
    ## get subpops - dataframe
    df = fp.retrieveMetaData(None, directory, outfn)

    #test
    #print(df)

    ## get list of subpops
    subdict = {}
    for pop in df['ethnic group'].unique():
        subpop = df[df['ethnic group'] == pop]
        subdict[pop] = list(subpop['id'])

    finalpopdict = defaultdict(list)
    for num in range(len(samplelist)):
        idname = samplelist[num]
        for key in subdict.keys():
            if idname in subdict[key]:
                finalpopdict[key].append(num)

    subpoplist = []
    for key in finalpopdict.keys():
        subpoplist.append(finalpopdict[key])

    ## calculate variance components
    a, b, c = allel.weir_cockerham_fst(g, subpoplist)

    ## average fst per variant
    fst = (np.sum(a, axis=1) / (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))

    return fst

def circos(directory, outfn, vcffile):
    ## make config file

    ## turn vcf file into .dat file
    ## chr - start - finish - snp density
    ## bash script

    ## create data file for heterozygosity
    snpdensity = pd.read_csv(directory + outfn + ".dat", sep='\t', header=None)

    ## using user input name read in vcf (either old vcf or new vcf)
    callset = allel.read_vcf(vcffile)
    pos = callset['variants/POS']
    chrm = callset['variants/CHROM']

    gt = allel.GenotypeArray(callset['calldata/GT'])

    ##test fst
    samplelist = callset['samples']
    fstlist = fst(gt, directory, outfn, samplelist)
    fp.makeDATFile(pos, gt, chrm, fstlist, snpdensity, directory, outfn, 'fst')

    ## count the heterozygoes for each pos
    hetcount = gt.count_het(axis=1)
    fp.makeDATFile(pos, gt, chrm, hetcount, snpdensity, directory, outfn, 'het')



def sfs(directory, vcffile):
    ### create a Site Frequency Spectrum Figure
    callset = allel.read_vcf(vcffile)

    gt = allel.GenotypeArray(callset['calldata/GT'])
    ac = gt.count_alleles()[:]

    derived = ac[:, 1]

    sfslist = allel.sfs(derived)

    xlabel = [x for x in range(1, len(sfslist)+1)]
    plt.style.use('seaborn-darkgrid')
    plt.plot(xlabel, list(sfslist), marker = 'o')
    for i in range(len(xlabel)):
        plt.text(xlabel[i], list(sfslist)[i],  str(list(sfslist)[i]), ha = 'center')

    plt.xlabel("K value")
    plt.ylabel("Number of variants")
    

    plt.savefig(vcffile +  "_sfs.jpg")


def singleChr():

    """
    To do: add code to allow for a single chromosome to output heterogeneity
    """
    ## plot the data and save
    plt.bar(x,bined_het)
    plt.xlabel(" Genetic Position ")
    plt.ylabel(" Number of Heterozygotes ")
    plt.savefig(directory + "/" + outfn + "_het.png")

def heatmap(directory, outfn, vcffile, logfile):

    ## calculate relatedness via vcftools
    #subprocess.Popen(['vcftools', '--vcf', vcffile, '--relatedness2', '--out', directory + "/" + outfn], stdout = logfile)
    subprocess.Popen(['vcftools', '--vcf', vcffile, '--relatedness2', '--out', directory + "/" + outfn])
    time.sleep(4)
    relFile = directory + "/" + outfn + ".relatedness2"
    time.sleep(4)
    ## call file parser to turn into csv
    fp.transposeRel(directory, relFile)
    matrixfile = directory + "/" + outfn + ".csv"

    ## call R script
    subprocess.Popen(['Rscript', '--vanilla', '/content/visualization_pipeline/command_line/heatmap.R', directory, matrixfile, directory + "/" + outfn + "_heatmap.jpg"])

def main(viz_options, directory, outfn, vcffile, colname):
    if 'circos' in viz_options:
        circos(directory, outfn, vcffile)
    if 'sfs' in viz_options:
        sfs(directory, vcffile)
    if 'pca' in viz_options:
        ## needs to be changed to allow user to select column
        #pca(directory, outfn, colname)
        pca(directory, vcffile, outfn)
    if 'heatmap' in viz_options:
        ## call R script
        heatmap(directory, outfn, vcffile, '')

if __name__ == '__main__':
    directory = sys.argv[1]
    vcffile = sys.argv[2]
    outfn = sys.argv[3]
    viz_options = [sys.argv[4]]

    main(viz_options, directory, outfn, vcffile, '')
