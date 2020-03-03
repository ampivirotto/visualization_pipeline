#-------------------------------------------------------------------------------
# Name:        vpp - make Sample Sheet
# Author:      tuk32868
#-------------------------------------------------------------------------------
import os, glob
import pandas as pd
import numpy as np

def getColumns(filename):
    """
    gets the id, sentrix barcode and sentrix position from file name
    """
    filename = os.path.basename(filename).strip(".idat")
    fileparts = filename.split("_")
    return fileparts

def makeSampleSheet(listoffiles):
    dict1 = {}
    for i in range(len(listoffiles)):
        newlist = getColumns(listoffiles[i])
        #if newlist[3] == "Re":
            #dict1[i] = [newlist[0], newlist[1], newlist[2]]
        dict1[i] = [newlist[0], newlist[1], newlist[2], newlist[3]]
    df = pd.DataFrame.from_dict(dict1, orient="index", columns=['sample_id', 'sentrix_barcode', 'sentrix_position', 'probe'])
    return df

#def pullGeoData(geonum):

def retrieveData(dirpath):
    idatfiles = [i for i in os.listdir(dirpath) if i.endswith(".idat")]
    return idatfiles

def main(dirpath):
    files = retrieveData(dirpath)
    df = makeSampleSheet(files)
    df.to_csv(dirpath + "\\samplesheet.csv", index=False)

if __name__ == '__main__':
    #main(dirpath)
    main("D:\\visualization_pipeline\\data\\beadarray_data")