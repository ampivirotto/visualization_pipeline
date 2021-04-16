import allel; print('scikit-allel', allel.__version__)
import random, subprocess
random.seed(14)
import time
import numpy as np
np.random.seed(14)
import h5py
import matplotlib.pyplot as plt
#import seaborn as sns  ## update causing errors?
#sns.set_style('white')
#sns.set_style('ticks')
import bcolz
import pandas as pd
from collections import defaultdict
import file_parser as fp
import sys
import pickle
import subprocess
import os


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

    callset = allel.read_vcf(vcffile)

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


def old_pca(directory, outfn, column, newVCF=False, samples = None, bs = 20000):
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

def pca(directory, vcffile, outfn):
    gds = vcffile.strip('vcf').strip('.') + '.gds'
    subprocess.run(['Rscript', '--vanilla', 'pca.R', directory, vcffile, outfn + '.html', gds])

def fst(g, directory, outfn, samplelist):
    ## get subpops - dataframe
    df = fp.retrieveMetaData(None, directory, outfn)

    ## get list of subpops
    subdict = {}
    for pop in df[df.columns[1]].unique():
        subpop = df[df[df.columns[1]] == pop]
        subdict[pop] = list(subpop['id'])

    finalpopdict = defaultdict(list)
    for num in range(len(samplelist)):
        idname = samplelist[num]
        for key in subdict.keys():
            for item in subdict[key]:
                if idname in item:
                    finalpopdict[key].append(num)

    subpoplist = []
    for key in finalpopdict.keys():
        subpoplist.append(finalpopdict[key])

    if len(subpoplist) > 1:
        ## calculate variance components
        a, b, c = allel.weir_cockerham_fst(g, subpoplist)

        ## average fst per variant
        fst = np.sum(a, axis=1) / (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1))
    else:
        fst = []

    return fst

def circos(directory, outfn, vcffile, chipType):
     
    ## get just vcffile name 
    vcffile = os.path.basename(vcffile)

    ## check directory format
    if not directory.endswith("/"):
        directory = directory+'/'
    cdir = directory + "circos/"

    if not os.path.exists(cdir):
        os.mkdir(cdir)

    ## identify karyotype and k id from karyotype from pickle file 
    with open('file_location.pickle', "rb") as f:
        chipDict = pickle.load(f)

    values = chipDict[chipType]
    ktype = values[3]
    kid = values[4]


    ## check vcf format and create an edited vcffile if not in correct format 
    newvcfname = fp.checkVCF(directory, vcffile)

    ## turn vcf file into .dat file
    ## chr - start - finish - snp density
    ## bash script
    command = "awk '/^#/ {next} {printf(\"%s\\t%d\\n\",$1,$2-$2%1000000);}' " + directory + newvcfname + " | sort | uniq -c | awk '{printf(\"" + kid +"%s\\t%s\\t%d\\t%s\\n\",$2,$3,$3+1000000,$1);}' > " + cdir + outfn + ".dat"
    with open(cdir + 'vcf2dat.sh', 'w') as o:
        o.write(command)
    subprocess.run(['bash', cdir + 'vcf2dat.sh'])


    ## read in snp density file created above 
    snpdensity = pd.read_csv(cdir + outfn + ".dat", sep='\t', header=None)

    snpdensity = fp.compareKaryotype(ktype, snpdensity, cdir + outfn + ".dat")

    ## using user input name read in vcf (either old vcf or new vcf)
    callset = allel.read_vcf(directory + newvcfname)
    pos = callset['variants/POS']
    chrm = callset['variants/CHROM']

    gt = allel.GenotypeArray(callset['calldata/GT'])

    ##test fst
    samplelist = callset['samples']
    fstlist = fst(gt, directory, outfn, samplelist)
    if len(fstlist) == 0:
        ## if there's only one population, then it will output an empty fst file.  
        f = open(cdir + outfn + "_fst.dat", 'w')
        f.close()
    else:
        fp.makeDATFile(pos, chrm, fstlist, snpdensity, cdir, outfn, 'fst')

    ## count the heterozygoes for each pos
    hetcount = gt.count_het(axis=1)
    #import pdb; pdb.set_trace()
    fp.makeDATFile(pos, chrm, hetcount, snpdensity, cdir, outfn, 'het')

    ## make config file
    fp.makeCircos(cdir, outfn, ktype)

    ## run circos 
    subprocess.run(["../software/circos/circos-0.69-9/bin/circos", "-outputdir", cdir, "-outputfile", outfn, "-conf", cdir + "circos.conf"])



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
    '''
    for i in range(len(xlabel)):
        plt.text(xlabel[i], list(sfslist)[i],  str(list(sfslist)[i]), ha = 'center')
    '''

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

def heatmap(directory, outfn, vcffile):

    ## calculate relatedness via vcftools
    subprocess.Popen(['vcftools', '--vcf', vcffile, '--relatedness2', '--out', directory + "/" + outfn])
    time.sleep(4)
    relFile = directory + "/" + outfn + ".relatedness2"
    time.sleep(4)
    ## call file parser to turn into csv
    fp.transposeRel(directory, relFile)
    matrixfile = directory + "/" + outfn + ".csv"

    ## call R script
    subprocess.Popen(['Rscript', '--vanilla', '/content/visualization_pipeline/command_line/heatmap.R', directory, matrixfile, directory + "/" + outfn + "_heatmap.jpg"])

def interactive_heatmap(directory, outfn, vcffile):
    ## calculate relatedness via vcftools
    subprocess.Popen(['vcftools', '--vcf', vcffile, '--relatedness2', '--out', directory + "/" + outfn])
    time.sleep(4)
    relFile = directory + "/" + outfn + ".relatedness2"
    time.sleep(4)
    ## call file parser to turn into csv
    fp.transposeRel(directory, relFile)
    matrixfile = directory + "/" + outfn + ".csv"

    subprocess.Popen(['Rscript', '--vanilla', '/content/visualization_pipeline/command_line/interactive_heatmap.R', directory, matrixfile])

def tree(directory, vcffile, outfn):
    f = open(directory +'/tree.nwk', 'w')
    subprocess.run(['vk','phylo','tree','upgma', vcffile], stdout = f)
    subprocess.run(['Rscript','--vanilla','tree_maker.R', directory+'/tree.nwk', directory, outfn+'_tree.png'])

def base_changes(directory, vcffile, outfn):
    process = subprocess.run(['bcftools', 'stats', vcffile], stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
    output = process.stdout
    output = output.decode("utf-8").split('\n')
    subs = output[47:59]
    subs = [i.split('\t') for i in subs]
    subs = [i[2:] for i in subs]
    A, C, G, T = subs[:3], subs[3:6], subs[6:9], subs[9:]
    A = [int(i[1]) for i in A]
    A.insert(0, 0)
    C = [int(i[1]) for i in C]
    C.insert(1,0) 
    G = [int(i[1]) for i in G] 
    G.insert(2,0)
    T = [int(i[1]) for i in T]
    T.insert(3,0)

    bars1, bars2, bars3, bars4 = [A[0],C[0],G[0],T[0]], [A[1],C[1],G[1],T[1]], [A[2],C[2],G[2],T[2]], [A[3],C[3],G[3],T[3]]
    barWidth = 0.15
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    r4 = [x + barWidth for x in r3]

    # Make the plot
    plt.bar(r1, bars1, color='green', width=barWidth, edgecolor='white', label='A')
    plt.bar(r2, bars2, color='orange', width=barWidth, edgecolor='white', label='C')
    plt.bar(r3, bars3, color='blue', width=barWidth, edgecolor='white', label='G')
    plt.bar(r4, bars4, color='red', width=barWidth, edgecolor='white', label='T')
    
    # Add xticks on the middle of the group bars
    plt.xlabel('Reference Base', fontweight='bold')
    plt.ylabel('Number of Changes', fontweight='bold')
    plt.title('Base Substitutions', fontweight='bold')
    plt.xticks([r + barWidth for r in range(len(bars1))], ['A', 'C', 'G', 'T'])
    
    # Create legend & Show graphic
    plt.legend()
    plt.style.use('seaborn-dark')
    plt.savefig(directory + '/' + outfn + '_baseChanges.png')

def tstv(directory, vcffile, outfn):
    process = subprocess.run(['bcftools', 'stats', vcffile], stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
    output = process.stdout
    output = output.decode("utf-8").split('\n')
    line = output[33]
    ts, tv, tstv = line.split('\t')[2], line.split('\t')[3], line.split('\t')[4]

    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    labels = 'Ts', 'Tv'
    sizes = [ts, tv]
    explode = (0.1, 0)

    fig1, ax1 = plt.subplots()
    ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
            shadow=True, startangle=90)
    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    ax1.set_title('Ts/Tv = ' + str(tstv))
    plt.savefig(directory + '/' + outfn + '_tstv.png')

def main(viz_options, directory, outfn, vcffile, chipName):
    if 'circos' in viz_options:
        circos(directory, outfn, vcffile, chipName)
    if 'sfs' in viz_options:
        sfs(directory, vcffile)
    if 'pca' in viz_options:
        pca(directory, vcffile, outfn)
    if 'interactive_heatmap' in viz_options:
        interactive_heatmap(directory, outfn, vcffile)
    if 'heatmap' in viz_options:
        ## call R script
        heatmap(directory, outfn, vcffile)
    if 'tree' in viz_options:
        tree(directory, vcffile, outfn)
    if 'base_changes' in viz_options:
        base_changes(directory,vcffile,outfn)
    if 'Ts/Tv' in viz_options:
        tstv(directory,vcffile,outfn)

        
if __name__ == '__main__':

    directory = sys.argv[1]
    vcffile = sys.argv[2]
    outfn = sys.argv[3]
    viz_options = [sys.argv[4]]
    chip = sys.argv[5]

    main(viz_options, directory, outfn, vcffile, chip)
