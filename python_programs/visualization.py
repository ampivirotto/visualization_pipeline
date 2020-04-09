import allel; print('scikit-allel', allel.__version__)
import random
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
#import toyplot
import os
import matplotlib.pyplot as plt
import seaborn as sns


logfile = open("testing.log", "a")

def randColor(num, pops):
    number_of_colors = num

    color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]

    colordict = {}
    for x in range(len(pops)):
        colordict[pops[x]] = color[x]

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
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    for pop in populations:
        flt = (sample_population == pop)
        ax.plot(x[flt], y[flt], marker='o', color=pop_colours[pop],
                label=pop, markersize=6, mec='k', mew=.5)
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))

def fig_pca(coords, model, title, sample_population=None):
    if sample_population is None:
        sample_population = df_samples.population.values
    populations = sample_population.unique()
    colordict = randColor(len(populations), populations)
    # plot coords for PCs 1 vs 2, 3 vs 4
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 2, 1)
    plot_pca_coords(coords, model, 0, 1, ax, sample_population, populations, colordict)
    ax = fig.add_subplot(1, 2, 2)
    plot_pca_coords(coords, model, 2, 3, ax, sample_population, populations, colordict)
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle(title, y=1.02)
    fig.tight_layout()

def summaryInfo(g, logfile):

    #get allele counts
    ac = g.count_alleles()[:]

    #count mutliallelic and singletons
    logfile.write("Multiallelic " + str(np.count_nonzero(ac.max_allele() > 1)))
    logfile.write("Singletons " + str(np.count_nonzero((ac.max_allele == 1) & ac.is_singleton(1))))

    return ac

def transform(g, ac):
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

    vcffile = directory + "/" + outfn + ".vcf"

    #runConversion(outfn, vcffile, bs)

    callsetfn = './analysis-vcf2hdf5/' + outfn + "_data.snps.hdf5"
    callset = h5py.File(callsetfn, mode= 'r')

    #get genotype data
    g = allel.GenotypeChunkedArray(callset['genos'])

    ac = summaryInfo(g, logfile)

    ## transform data
    gn = transform(g, ac)

    return gn

def retrieveMetaData(samples, directory, outfn):
    ## make empty dictionary
    metadata = {}

    if samples is None:
        samples = [ f.name for f in os.scandir(directory) if f.is_dir() ]
    for x in range(len(samples)-1):
        if not samples[x].startswith("GS"):
            samples.remove(samples[x])

    ## loop through all samples included in vcf

    for sample in samples:
        temp = []
        for file in os.listdir(directory + "/" + sample):
            if file.endswith(".txt"):
                with open(directory + "/" + sample + "/" + file) as f:
                    for line in f:
                        if line.startswith(" - character"):
                            temp = line.split(":")
                            break
                empty = []
                for item in temp:
                    if item.startswith(" - "):
                        continue
                    else:
                        xx = item.split(",")
                        for x in xx:
                            empty.append(x)
                counter = 1
                metalist = []
                for val in empty:
                    if counter % 2 == 0 :
                        metalist.append(val.strip())
                    counter +=1
                metadata[sample] = metalist
                break

    df = pd.DataFrame.from_dict(metadata, orient='index')
    df.to_csv(outfn + ".csv")
    return df

def LD(directory, outfn, newVCF=False, samples = None, bs = 20000):
    """
    graph the linkage disequilibrium graphs
    """
    gn = prepData(directory, outfn, newVCF, samples, bs)

    #plot ld
    plot_ld(gn[:1000], 'Figure 1. Pairwise LD.')
    plt.show()

    #subsabmple snps
    n = 100000  # number of SNPs to choose randomly
    vidx = np.random.choice(gn.shape[0], n, replace=False)
    vidx.sort()
    gnr = gn.take(vidx, axis=0)

    plot_ld(gnr[:1000], 'Figure 2. Pairwise LD after random downsampling.')
    plt.show()

    #prune linkage
    gnu = ld_prune(gnr, size=500, step=200, threshold=.1, n_iter=5)

    #linkage after pruning
    plot_ld(gnu[:1000], 'Figure 3. Pairwise LD after LD pruning.')
    plt.show()


def pca(directory, outfn, newVCF=False, samples = None, bs = 20000):
    """
    main function to run pca visualization
    """

    gn = prepData(directory, outfn, newVCF, samples, bs)

    ## get metadata
    df = retrieveMetaData(samples, directory, outfn)

    coords1, model1 = allel.pca(gn, n_components=10, scaler='patterson')

    fig_pca(coords1, model1, 'Figure 4. Conventional PCA.', sample_population = df[1])
    plt.show()

    allel.pca(gn)

    ## init PCA  https://ipyrad.readthedocs.io/en/latest/API-analysis/cookbook-pca.html
    ##pca = ipa.pca(
    ##    data="./analysis-vcf2hdf5/human_popgen.snps.hdf5",
    ##    mincov=1.0,
    ##    impute_method='None'
    ##)
    ##
    #### run pCA
    ##pca.run(nreplicates=25, seed=14)
    ##pca.draw(0, 1)

def heterozygosity(directory, outfn, newVCF=False, samples=None, bs=20000):
    """
    graph the heterozygosity across genome in bar graph
    """
    ## check if user wants to make a new vcf from new sample subset
    if newVCF == True:
        makeVCF(directory, samples, outfn)

    ## using user input name read in vcf (either old vcf or new vcf)
    callset = allel.read_vcf(directory + "/" + outfn + '.vcf')
    #print(sorted(callset.keys()))
    #print(callset['samples'])

    ## get gentoypes
    gt = allel.GenotypeArray(callset['calldata/GT'])

    ## count the heterozygoes for each pos
    hetcount = gt.count_het(axis=1)

    ## put values in bins
    ## TO DO: update this to reflect actual positions on chr
    bined_het = []
    counter = 1
    total = 0
    for val in hetcount:
        total+=val
        if counter == bs:
            bined_het.append(total/bs)
            counter = 1
            total = 0
        else:
            counter+=1

    ## get x values
    x = [val*bs for val in range(1, len(bined_het)+1)]

    ## plot the data and save
    plt.bar(x,bined_het)
    plt.xlabel(" Genetic Position ")
    plt.ylabel(" Number of Heterozygotes ")
    plt.savefig(directory + "/" + outfn + "_het.png")
