import re, os, sys
def change_idat_names(accession):    
    '''iaap-cli only accepts IDAT files in the form of XXXX_XXXX_Grn.idat or XXXX_XXXX_Red.idat
        This function renames them to work with iaap-cli. accession is the GEO accession number.'''
    
    print('changing idat names')
    directory = '/content/data/' + accession + '/'
    #os.chdir(directory)    
    sample_folders = [i for i in os.listdir(directory) if os.path.isdir(directory + i)]
    for folder in sample_folders:
        idats = [i for i in os.listdir(directory + folder) if i.endswith('.idat')]
        for idat in idats:
            if '_Red' in idat:
                new_name = re.sub('_R.*', '_Red.idat', idat)
                os.rename(directory + folder + '/' + idat, directory + folder + '/' + new_name)
            else:
                new_name = re.sub('_R.*', '_Grn.idat', idat)
                os.rename(directory + folder + '/' + idat, directory + folder + '/' + new_name)
            #print(new_name)
change_idat_names(sys.argv[1])
