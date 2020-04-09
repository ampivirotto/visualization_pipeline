import ipyrad.analysis as ipa
import pandas as pd
import toyplot
import os


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

def pca(directory, outfn, newVCF=False, samples = None, bs = 20000):
    """
    main function to run pca visualization
    """
    ## can either use an existing vcf file or make a new vcf file from a list of samples
    if newVCF == True:
        makeVCF(directory, samples, outfn)
        vcffile = directory + "/" + outfn + ".vcf"
    else:
        vcffile = directory + "/" + outfn + ".vcf"

    runConversion(outfn, vcffile, bs)
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