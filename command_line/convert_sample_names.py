import os, subprocess, re

acc = 'GSE83225'
path = '/content/data/'+acc+'/'
names = subprocess.check_output(['bcftools','query','-l','/content/data/GSE83225/test.vcf'])
names = names.decode("utf-8").split('\n')[:-1]
print(names)
breeds = []
for i in [x for x in os.listdir('/content/data/'+acc) if x.startswith('GSM') and os.path.isdir(f'/content/data/{acc}/{x}')]:
    for f in os.listdir(path+i):
        if f.endswith('.txt'):
            with open(path+i+'/'+f) as text:
                breed = re.findall('breed.+', text.read())[0]
                print(breed)
