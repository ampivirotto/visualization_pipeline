from collections import defaultdict
import pandas as pd
import numpy as np
import math
import os
import subprocess

def findmax(directory):
    """
    get max values from snp density, heterozygote, and fst file 
    """
    for filename in os.listdir(directory):
        if filename.endswith("_het.dat"):
            df = pd.read_csv(directory + filename, sep = "\t", header = None)
            maxhet = df[3].max() * 1.1
        elif filename.endswith("_fst.dat"):
            try:
                df = pd.read_csv(directory + filename, sep = "\t", header = None)
                maxfst = df[3].max() * 1.1
            except:
                maxfst = 0
        elif filename.endswith(".dat"):
            df = pd.read_csv(directory + filename, sep = "\t", header = None)
            maxden = df[3].max() * 1.1 
    return maxden, maxhet, maxfst


def makeCircos(directory, gseid, ktype):
    """
    make .conf file 
    """
    ## get max values from .dat files 
    maxden, maxhet, maxfst = findmax(directory)

    ## write to .conf file 
    with open(directory + '/circos.conf', 'w') as o:
        ## karyotype file 
        o.write('karyotype = data/karyotype/{}\n\n'.format(ktype))
        ## start set up of plot
        o.write("chromosomes_units = 1000000\n\n<plots>\n\n")
        ## snp density plot 
        o.write("<plot>\ntype = histogram\nmin=0\nmax={}\nfile= {}.dat\nr0  = 0.76r\nr1  = 0.95r\ncolor = black_a4\nfill_color = blue\nthickness = 2\n</plot>\n".format(maxden, gseid))
        ## heterozygosity plot
        o.write("<plot>\ntype = histogram\nmin=0\nmax={}\nfile= {}_het.dat\nr0  = 0.56r\nr1  = 0.75r\ncolor= black_a4\nfill_color = green\nthickness = 1\n</plot>\n".format(maxhet, gseid))
        ## fst plot
        o.write("<plot>\ntype = histogram\nmin=0\nmax={}\nfile = {}_fst.dat\nr0  = 0.35r\nr1  = 0.55r\ncolor = black_a4\nfill_color = red\nthickness = 1\n</plot>\n".format(maxfst, gseid))
        ## finish set up 
        o.write("</plots>\n<<include etc/ideogram.conf>>\n<<include etc/ticks.conf>>\n<image>\n<<include etc/image.conf>>\n</image>\n<<include etc/colors_fonts_patterns.conf>>\n<<include etc/housekeeping.conf>>")

def compareKaryotype(kfname, sddf, sdfpath):
    """
    takes the karyotype file and SNP density DAT file and compares the naming scheme to ensure they're the same
    for circos plotting  
    """
    ## read in karyotype file 
    df = pd.read_csv("../software/circos/circos-0.69-9/data/karyotype/" + kfname, sep = " ", header = None)  ## have to update where the file is located
    ## get the chromosome ids used in karyotype file 
    kt_ids = list(df[2])

    ids = sddf[0].unique()
    
    found = []
    differs = []

    for idnum in ids:
        if idnum in kt_ids:
            found.append(idnum)
        else:
            differs.append(idnum)
    if len(differs) > 1:
        raise Exception('Check Karyotype and VCF chromosome naming schemes')
    elif len(differs) == 1:
        for item in kt_ids:
            if not item in found:
                temp = sddf[sddf[0] == differs[0]]
                eelse = sddf.loc[sddf[0] != differs[0]]
                temp.loc[:,0] = item

                final = pd.concat([temp, eelse])
                final.to_csv(sdfpath, index = False, header = False, sep = "\t")
                #print(len(sddf), len(eelse) + len(temp))
                return final
    else:
        return sddf


def checkVCF(directory, vcffile):
    """
    check to see if VCF has chr# or # in chr column.  if has chr# then remove the string 'chr' and output to edited vcf file. 
    """
    with open(directory + vcffile) as f:
        for line in f:
            if not line.startswith("#"):
                temp = line.split("\t")
                chrn = temp[0]

                if 'chr' in chrn:
                    newvcffile = vcffile[:-4] + "_edit.vcf"
                    command = "awk \'{gsub(/^chr/,\"\"); print}\' " + directory + vcffile + " > " + directory + newvcffile
                    with open(directory + 'fixVCF.sh', 'w') as o:
                        o.write(command)
                    subprocess.run(['bash', directory + 'fixVCF.sh'])
                    break
                else:
                    newvcffile = vcffile
    return newvcffile


def transposeRel(directory, relFile):
    ## set up dictionary
    cols = []
    maindict = defaultdict(list)
    with open(relFile) as f:
        total_samples = int(math.sqrt(len(f.read().split('\n'))-1))
        
    ## read through output file to parse information
    with open(relFile) as f:
        
        for line in f:
            #print(line)
            if line.startswith("GSM"):
                #print(line)
                splitline = line.split("\t")
                try:
                    indexnum = cols.index(splitline[1])
                except:
                    cols.append(splitline[1])
                    indexnum = cols.index(splitline[1])

                if splitline[0] in maindict.keys():
                    maindict[splitline[0]].insert(indexnum, float(splitline[6].strip("\n")))
                    maindict[splitline[0]].pop(indexnum+1)
                else:
                    maindict[splitline[0]] = [math.nan for x in range(1, total_samples+1)]
                    maindict[splitline[0]].insert(indexnum, float(splitline[6].strip("\n")))
                    maindict[splitline[0]].pop(indexnum+1)

    ## turn into df and csv
    df = pd.DataFrame.from_dict(maindict, orient='index', columns=cols)
    ## remove the .relatedness2 from filename 
    filename = relFile[:-13]
    df.to_csv(filename + ".csv")



def makeDATFile(pos, chrm, typecount, snpdensity, directory, outfn, typel):

    ## creat empty dictionary to store the varaints per chromosome
    typedict = defaultdict(list)

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
    newdict = {}

    row = 0
    for x in range(len(start)):
        total = 0
        counter = 0
        try:
            snps = typedict[int(chrnm[x][2:])]
        except:
            snps = typedict[chrnm[x][2:]]
        for var in snps:
            pos = list(var.keys())[0]
            if pos < end[x]:
                if pos > start[x]:
                    total += list(var.values())[0]
                    counter +=1
            else:
                break 
        if total == 0 and counter == 0:
            newdict[row] = [chrnm[x], start[x], end[x], 0]
        else:
            newdict[row] = [chrnm[x], start[x], end[x], total/counter]
        row += 1

    ## take the snpdensity file and drop snp density and add heterozygosity
    #pdb.set_trace()
    dat = pd.DataFrame.from_dict(newdict, orient = 'index')
    dat.to_csv(directory + outfn + "_" + typel + ".dat", index=False, header=False, sep="\t")

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
    ## might need directory - check in debugging 
    #df.to_csv(directory + outfn + ".csv")
    df.to_csv(directory + outfn + ".csv")
    return df
