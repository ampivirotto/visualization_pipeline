This bash script has several different functions:

Line 2: export path to bcftools plug in

Line 4: change to user specified directory

Line 5: unzip all files in that directory NOTE: should be warning if file is large?

Line 7-21: move all the .idat files to their corresponding sample folder

Line 23-36: move all the metadata .txt files to their corresponding sample folder

Line 38-41: for each sample folder run autoconvert on the idat files within 

Line 45: make a list of all gtc files within the subdirectories NOTE: allow users to select a subset?

Line 48: use bcftools with the gtc2vcf plugin to convert gtc files into a single vcf file, sort, normalize, and index the vcf file