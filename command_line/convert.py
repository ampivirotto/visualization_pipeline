import os, sys, subprocess
from multiprocessing import Pool

def convert(folders):
    global accession
    for folder in folders:
        cmd = '/content/iaap-cli/./iaap-cli gencall /content/illumina_files/CanineHD_B.bpm /content/illumina_files/CanineHD_A.egt'
        subprocess.call(cmd)

        
        print('done')

accession = sys.argv[1]
print( os.listdir('/content/data/' + accession))
processes = 4
p = Pool(processes)
directories = [i for i in os.listdir('/content/data/' + accession) if os.path.isdir('/content/data/' + accession+'/'+i)]
split_directories = []

maxx = len(directories)
while directories:
    split_directories.append(directories[:maxx//processes])
    directories = directories[maxx//processes:]
print(split_directories)
p.map(convert, split_directories)


