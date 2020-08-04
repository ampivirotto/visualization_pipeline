from collections import defaultdict
import pandas as pd
import numpy as np
import math

def transposeRel(directory, relFile):
    ## set up dictionary
    cols = []
    maindict = defaultdict(list)

    ## read through output file to parse information
    with open(relFile) as f:
        for line in f:
            if line.startswith("GSM"):
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
                    maindict[splitline[0]] = [math.nan for x in range(1, 241)]
                    maindict[splitline[0]].insert(indexnum, float(splitline[6].strip("\n")))
                    maindict[splitline[0]].pop(indexnum+1)

    ## turn into df and csv
    df = pd.DataFrame.from_dict(maindict, orient='index', columns=cols)
    filename = relFile.strip(".relatedness2")
    df.to_csv(directory + "/" + filename + ".csv")



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