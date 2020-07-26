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
import ipyrad.analysis as ipa
import pandas as pd
import os
from collections import defaultdict


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

def retrieveMetaData(samples, directory, outfn):
    ## make empty dictionary
    metadata = {}

    if samples is None:
        samples = [ f.name for f in os.scandir(directory) if f.is_dir() ]

    remove = []
    for x in range(len(samples)-1):
        try:
            if not samples[x].startswith("GS"):
                remove.append(samples[x])
        except:
            print("ERROR")
    for removal in remove:
        samples.remove(removal)

    ## loop through all samples included in vcf

    counter = 0
    cols = []
    for sample in samples:
        metalist = []
        for file in os.listdir(directory + "/" + sample):
            if file.endswith(".txt"):
                with open(directory + "/" + sample + "/" + file) as f:
                    for line in f:
                        if line.startswith(" - character"):
                            temp = line.strip('\n').split(":")
                            #print(line)

                            empty = []
                            for x in range(0, len(temp)-1):
                                temp[x] = temp[x] + ":"
                            for item in temp:
                                if item.startswith(" - "):
                                    continue
                                else:
                                    xx = item.split(",")
                                    for x in xx:
                                        empty.append(x.lstrip())
                            for val in range(0, len(empty)):
                                if counter == 0:
                                    if empty[val].endswith(":"):
                                        cols.append(empty[val].strip(":"))
                                        location = cols.index(empty[val].strip(":"))
                                    else:
                                        try:
                                            metalist.insert(location, metalist[location] + "," + empty[val])
                                            metalist.pop()
                                        except IndexError:
                                            metalist.insert(location, empty[val])
                                else:
                                    if empty[val].endswith(":"):
                                        if empty[val].strip(":") in cols:
                                            location = cols.index(empty[val].strip(":"))
                                        else:
                                            cols.append(empty[val].strip(":"))
                                            location = cols.index(empty[val].strip(":"))
                                    else:
                                        try:
                                            metalist.insert(location, metalist[location] + "," + empty[val])
                                            metalist.pop()
                                        except IndexError:
                                            metalist.insert(location, empty[val])
                        elif line.startswith(" - source"):
                            if counter == 0:
                                cols.append('source')
                            line = line.split(" : ")
                            location = cols.index('source')
                            metalist.insert(location, line[1].strip())
                        elif line.startswith(" - supp"):
                            if counter == 0:
                                cols.append('id')
                            line = line.split(" : ")
                            files = line[1].split(",")
                            parts = files[1].split("/")
                            name = parts[8][:-13]
                            location = cols.index('id')
                            metalist.insert(location, name)



                metadata[sample] = metalist
                counter += 1
                break


    df = pd.DataFrame.from_dict(metadata, orient='index', columns=cols)
    df.to_csv(outfn + ".csv")
    return df

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


def pca(directory, outfn, column, newVCF=False, samples = None, bs = 20000):
    """
    main function to run pca visualization
    """
    import pdb

    #gn, callset = prepData(directory, outfn, newVCF, samples, bs)
    callset = allel.read_vcf(directory + outfn + ".vcf")

    g = allel.GenotypeChunkedArray(callset['calldata/GT'])
    gn = transform(g)

    ## get metadata
    #df = retrieveMetaData(samples, directory, outfn)
    df = pd.read_csv(directory + "canine_subset_7_8.csv")

    coords1, model1 = allel.pca(gn, n_components=10, scaler='patterson')

    fig_pca(directory, outfn, coords1, model1, 'Conventional PCA.', sample_population = df[column])



def makeDATFile(pos, gt, chrm, typecount, snpdensity, directory, outfn, typel):

    ## creat empty dictionary to store the varaints per chromosome
    typedict = {}

    ## add all chromosomes to the dictionary
    for chrnum in range(1, 23):
        typedict[chrnum] = []
    typedict['X'] = []
    typedict['Y'] = []
    typedict['MT'] = []

    ## add variants with with item from count to temp dictionary
    for x in range(len(pos)):
        temp ={pos[x]:typecount[x]}
        try:
            typedict[int(chrm[x])].append(temp)
        except:
            typedict[chrm[x]].append(temp)

    ## bin them based on original dat file
    ## To change allow user to differ these bin sizes for different tracts
    chrnm = snpdensity[0]
    start = snpdensity[1]
    end = snpdensity[2]
    typelist = []

    for x in range(len(start)):
        total = 0
        counter = 0
        try:
            snps = typedict[int(chrnm[x][2:])]
        except:
            snps = typedict[chrnm[x][2:]]
        for var in snps:
            try:
                pos = list(var.keys())[0]
                if pos < end[x]:
                    if pos > start[x]:
                        total+= list(var.values())[0]
                        counter+=1
                else:
                    break
            except:
                pos = list(var.keys())[0]
        if total ==0 and counter == 0:
            typelist.append(0)
        else:
            try:
                typelist.append(total/counter)
            except:
                print("ERROR")

    ## take the snpdensity file and drop snp density and add heterozygosity
    snpdensity.drop([3], axis=1)
    snpdensity[3] = typelist
    snpdensity.to_csv(directory + outfn + "_" + typel + ".dat", index=False, header=False, sep="\t")

def fst(g, directory, outfn, samplelist):
    ## get subpops - dataframe
    df = retrieveMetaData(None, directory, outfn)

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
    callset = allel.read_vcf(directory + vcffile + ".vcf")
    pos = callset['variants/POS']
    chrm = callset['variants/CHROM']

    gt = allel.GenotypeArray(callset['calldata/GT'])

    ##test fst
    samplelist = callset['samples']
    fstlist = fst(gt, directory, outfn, samplelist)
    makeDATFile(pos, gt, chrm, fstlist, snpdensity, directory, outfn, 'fst')

    ## count the heterozygoes for each pos
    hetcount = gt.count_het(axis=1)
    makeDATFile(pos, gt, chrm, hetcount, snpdensity, directory, outfn, 'het')



def sfs(directory, vcffile):
    ### create a Site Frequency Spectrum Figure
    callset = allel.read_vcf(directory + vcffile + ".vcf")

    gt = allel.GenotypeArray(callset['calldata/GT'])
    ac = gt.count_alleles()[:]

    derived = ac[:, 1]

    sfslist = allel.sfs(derived)

    xlabel = [x for x in range(1, len(sfslist)+1)]

    plt.plot(xlabel, list(sfslist))
    plt.xlabel("K value")
    plt.ylabel("Number of variants")
    plt.savefig(directory + "/" + vcffile +  "_sfs.jpg")


def singleChr():

    """
    To do: add code to allow for a single chromosome to output heterogeneity
    """
    ## plot the data and save
    plt.bar(x,bined_het)
    plt.xlabel(" Genetic Position ")
    plt.ylabel(" Number of Heterozygotes ")
    plt.savefig(directory + "/" + outfn + "_het.png")


def main(viz_options, directory, outfn, vcffile):
    if 'circos' in viz_options:
        circos(directory, outfn, vcffile)
    if 'sfs' in viz_options:
        sfs(directory, vcffile)
    if 'pca' in viz_options:
        ## needs to be changed to allow user to select column
        pca(directory, outfn, 'ethnic group')
    if 'heatmap' in viz_options:
        ## call R script
        print('error')


## TESTING ###
#pca("../../data/GSE74100/", 'myanmar_population_structure_GSE74100', 'ethnic group')
#circos("../../data/GSE74100/", 'GSE74100', "myanmar_population_structure_GSE74100")
#retrieveMetaData(None, '../../data/GSE74100/', 'GSE74100')
#sfs("../../data/GSE74100/", "myanmar_population_structure_GSE74100")

##listdf = []
##directs = ['GSE96736', 'GSE90441', 'GSE83160', 'GSE70454']
##for x in directs:
##    df = retrieveMetaData(None, '/mnt/d/visualization_pipeline/data/' + x + "/", x)
##    listdf.append(df)
##
##df_all = pd.concat(listdf, sort=False)
##
##df_all.to_csv("../data/sardinia_canine_allSamples.csv", index=False)

pca("../data/", 'subset_7_8.recode', 'breed')

#runConversion('sardinia_canine_allSamples', '../data/sardinia_canine_allsamples.vcf', 10000)