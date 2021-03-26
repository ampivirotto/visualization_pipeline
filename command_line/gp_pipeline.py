import GEOparse, os, wget, subprocess, sys, time, pickle

#import visualization as viz_py

def identifyChip(chipType):
    """
    used to identify the associated chip files based on user input
    """
    '''
    #this part didn't work; illumina_dict.pickle not found
    with open('../illumina_files/illumina_dict.pickle', "rb") as f:
        chipDict = pickle.load(f)
    values = chipDict[chipType]
    print('BPM: ' + values[0] + '\n')
    print('EGT: ' + values[1] + '\n')
    print('CSV: ' + values[2] + '\n')
    return values[0], values[1], values[2]
    '''
    return 'CanineHD_B.bpm', 'CanineHD_A.egt', 'CanineHD_B.csv'
def checkDir(directory):
    """
    identify whether the dictionary exists or not - if it doens't make one
    create the LOG file
    """
    ## test if directory is there
    if not os.path.exists(directory):
        os.mkdir(directory)
        orig_stout = sys.stdout
        f = open(directory + '/' + str(time.time()) + '.log', 'w')
        sys.stdout = f
        print("Making new directory: " + directory + "\n")
    else:
        orig_stout = sys.stdout
        f = open(directory + '/' + str(time.time()) + '.log', 'w')
        sys.stdout = f
        print("Found directory: " + directory + "\n")
    return orig_stout, f

def retrieveGEOFiles(geonum, directory):
    """
    using GEOparse module - pull down data using GEO number supplied by user
    make files of the metadata along with platform.
    donwload supplment files
    TO DO: make a csv of metadata for use in PCA
    """
    samplelist = []

    ##download data  https://geoparse.readthedocs.io/en/latest/GEOparse.html
    print("###############STARTING DOWNLOAD################ \n\n\n")
    print(geonum + '\n\n')
    gse = GEOparse.get_GEO(geo=geonum, destdir=directory)

    for gsm_name, gsm in gse.gsms.items():

        samplelist.append(gsm_name)

        filename = directory + "/" + gsm_name + ".txt"
        o = open(filename, "w")
        o.write("Name: " + gsm_name)
        o.write("\nMetadata:")
        for key, value in gsm.metadata.items():
            o.write("\n - %s : %s" % (key, ", ".join(value)))
            if key == 'supplementary_file':
                for item in value:
                    wget.download(item, directory)
        o.close()

    for gpl_name, gpl in gse.gpls.items():
        filename = directory + "/" + gpl_name + ".platform"
        o = open(filename, "w")
        o.write("Name: " + gpl_name)
        o.write("\nMetadata:")
        for key, value in gpl.metadata.items():
            o.write("\n - %s : %s" % (key, ", ".join(value)))
        o.close()

    print("################### FINISHED DOWNLOAD ###################### \n\n")

    return samplelist

def makeBashFile(directory, bpm, csv, egt, output):
    """
    write the bash file to run bcftools
    TO DO: adjust to where things are saved and what things are saved following new project management
    """
    ## write bash file
    print("Making Bash File ... \n\n")
    bash = open(directory + '/run1.sh', "w")
    bash.write("direct=\'" + directory + "\'\n")
    bash.write("bpm=\'" + bpm + "\'\n")
    bash.write("egt=\'" + egt + "\'\n")
    bash.write("csv=\'" + csv + "\'\n")
    bash.write("output=\'" + output + "\'\n\n")
    bash.close()

    ## mash bash files
    filenames = [directory + '/run1.sh', 'main.sh']
    with open(directory + '/final.sh', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                outfile.write(infile.read())
    print("Finished making Bash File... \n\n")

def runBash(directory):
    '''
    run bash file just created
    TO DO: make bash read.me to describe what it does
    '''
    
    ## run bash files
    file = directory + "/final.sh"
    #command = "bash " + file
    subprocess.call(['dos2unix', file]) #does not work in linux
    subprocess.call(['bash', file])

def main(geonum, chipType, allSamples, output, viz=None):

    start = time.time()
    
    directory = "/content/data/" + geonum

    ogstdout, f = checkDir(directory)
    
    bpm, egt, csv = identifyChip(chipType)
    
    print('retrieveGEOFiles')
    samplelist = retrieveGEOFiles(geonum, directory)
    
    
    
    makeBashFile(directory, bpm, csv, egt, output)

    runBash(directory)

##    if allSamples == False:
##        files = [(x,x) for x in samplelist]
##        #else introduce a way for user to select files
##        return files
##    else:
##        viz_py.main(viz, directory, output + "_" + geonum, output)

    end = time.time()
    print("Time elapsed: " + str(end-start))

    sys.stdout = ogstdout
    f.close()



geonum = sys.argv[1]
chipType = sys.argv[2]
output = sys.argv[3]
allSamples = sys.argv[4]

main(geonum, chipType, allSamples, output)



