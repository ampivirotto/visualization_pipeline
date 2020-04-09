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
import pandas

##def plot_ld(gn, title):
##    m = allel.rogers_huff_r(gn) ** 2
##    ax = allel.plot_pairwise_ld(m)
##    ax.set_title(title)
##
##def ld_prune(gn, size, step, threshold=.1, n_iter=1):
##    for i in range(n_iter):
##        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
##        n = np.count_nonzero(loc_unlinked)
##        n_remove = gn.shape[0] - n
##        print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
##        gn = gn.compress(loc_unlinked, axis=0)
##    return gn
##
##def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population):
##    sns.despine(ax=ax, offset=5)
##    x = coords[:, pc1]
##    y = coords[:, pc2]
##    for pop in populations:
##        flt = (sample_population == pop)
##        ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[pop],
##                label=pop, markersize=6, mec='k', mew=.5)
##    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
##    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))
##
##
##def fig_pca(coords, model, title, sample_population=None):
##    if sample_population is None:
##        sample_population = df_samples.population.values
##    # plot coords for PCs 1 vs 2, 3 vs 4
##    fig = plt.figure(figsize=(10, 5))
##    ax = fig.add_subplot(1, 2, 1)
##    plot_pca_coords(coords, model, 0, 1, ax, sample_population)
##    ax = fig.add_subplot(1, 2, 2)
##    plot_pca_coords(coords, model, 2, 3, ax, sample_population)
##    ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
##    fig.suptitle(title, y=1.02)
##    fig.tight_layout()
##
##callset_fn = './analysis-vcf2hdf5/popgen_data.snps.hdf5'
##
##callset = h5py.File(callset_fn, mode='r')
##
##print(callset.keys())
##
####get genotype data
##g = allel.GenotypeChunkedArray(callset['genos'])
##
#### get allele counts
##ac = g.count_alleles()[:]
##
##print(g)
##print(ac)
##
#### count mutliallelic and singletons
#### add print statements here to a log file
##np.count_nonzero(ac.max_allele() > 1)
##np.count_nonzero((ac.max_allele == 1) & ac.is_singleton(1))
##
#### remove singletons and multiallelic
##flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
##gf = g.compress(flt, axis=0)
##print(gf)
##
#### transform into 2d matrix
##gn = gf.to_n_alt()
##print(gn)
##
#### plot ld
##plot_ld(gn[:1000], 'Figure 1. Pairwise LD.')
##plt.show()
##
#### take fewer snps
##n = 100000  # number of SNPs to choose randomly
##vidx = np.random.choice(gn.shape[0], n, replace=False)
##vidx.sort()
##gnr = gn.take(vidx, axis=0)
##
##plot_ld(gnr[:1000], 'Figure 2. Pairwise LD after random downsampling.')
##plt.show()
##
#### prune linkage
##gnu = ld_prune(gnr, size=500, step=200, threshold=.1, n_iter=5)
##
#### linkage after pruning
##plot_ld(gnu[:1000], 'Figure 3. Pairwise LD after LD pruning.')
##plt.show()
##
##coords1, model1 = allel.pca(gnu, n_components=10, scaler='patterson')



### insert meta data - read in csv, find unique populations, set colors

##fig_pca(coords1, model1, 'Figure 4. Conventional PCA.')
##plt.show()

callset = allel.read_vcf('../data/GSE93106/hmanomni_test.vcf')

print(sorted(callset.keys()))

print(callset['samples'])

gt = allel.GenotypeArray(callset['calldata/GT'])

print(gt)

hetcount = gt.count_het(axis=1)

##for item in hetcount:
##    if item > 0:
##        print(item)

import matplotlib.pyplot as plt
import seaborn as sns

sns.distplot(hetcount)
plt.show()
