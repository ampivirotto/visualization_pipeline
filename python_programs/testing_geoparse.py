import GEOparse
import os
import wget
import subprocess
import sys
import time

def identifyChip(chipType):
    """
    used to identify the associated chip files based on user input
    TO DO: have information saved in a dictionary and open dictionary when this is called
    """
    chipDict = {'human_omni_express':['HumanOmniExpress-24-v1-0-B.bpm', 'HumanOmniExpress_24v1-0_A.egt', 'HumanOmniExpress-24-v1-0-B.csv'],
                'inf_omni_zhonghua':['InfiniumOmniZhongHua-8v1-4_A1.bpm', 'InfiniumOmniZhongHua-8v1-4_A1.egt', 'InfiniumOmniZhongHua-8v1-4_A1.csv'],
                'human_omni': ['HumanOmni5-4v1-1_A.bpm', 'HumanOmni5-4v1-1.egt', 'HumanOmni5-4v1-1_A.csv']}

    values = chipDict[chipType]

    print('BPM: ' + values[0] + '\n')
    print('EGT: ' + values[1] + '\n')
    print('CSV: ' + values[2] + '\n')

    return values[0], values[1], values[2]

def checkDir(directory):
    """
    identify whether the dictionary exists or not - if it doens't make one
    create the LOG file

    """
    ## test if directory is there
    if not os.path.exists(directory):
        os.mkdir(directory)
        sys.out = open(directory + '/' + str(time.time()) + '.log', 'w')
        print("Making new directory: " + directory + "\n")
    else:
        sys.out = open(directory + '/' + str(time.time()) + '.log', 'w')
        print("Found directory: " + directory + "\n")

def retrieveGEOFiles(geonum, directory):
    """
    using GEOparse module - pull down data using GEO number supplied by user
    make files of the metadata along with platform.
    donwload supplment files
    TO DO: make a csv of metadata for use in PCA
    """
    ##download data  https://geoparse.readthedocs.io/en/latest/GEOparse.html
    print("###############STARTING DOWNLOAD################ \n\n\n")
    print(geonum + '\n\n')
    gse = GEOparse.get_GEO(geo=geonum, destdir=directory)

    for gsm_name, gsm in gse.gsms.items():
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

    print(" ################### FINISHED DOWNLOAD ###################### \n\n")

def makeBashFile(directory, bpm, csv, egt):
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
    subprocess.call(['dos2unix', file])
    subprocess.call(['bash', file])

def main(geonum, directory, chipType):

    checkDir(directory)

    bpm, egt, csv = identifyChip(chipType)

    retrieveGEOFiles(geonum, directory)

    makeBashFile(directory, bpm, csv, egt)

    runBash(directory)

    """
    insert here the graphing portion - allow user to select graphing options
    """

    sys.out.close()


if __name__ == '__main__':
##    geonum = sys.argv[1]
##    userdirect = sys.argv[2]
##    chipType = sys.argv[3]
##    output = sys.argv[4]

    geonum = 'GSE93106'
    output='popgen.vcf'
    chipType = 'human_omni_express'
    userdirect = "/mnt/d/visualization_pipeline/data/"

    directory = userdirect+geonum

    main(geonum, directory, chipType)


