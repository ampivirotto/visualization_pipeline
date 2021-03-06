
ref='/mnt/d/visualization_pipeline/reference/Canis_familiaris.BROADD2.67.dna.toplevel.fa'
export BCFTOOLS_PLUGINS="/home/tuk32868/bin"  ## make sure plugins are set 

cd $direct  ## change to main directory with files 
gunzip *.gz  ## unzip all the files 

for f in *.idat           ## for each idat file in the directory 
do
    filename=${f##*/}          			## Use the last part of a path
    extension=${f##*.}         			## Remove up to the last dot to get extension
    dir=${filename%_*}         			## remove end with GRN or Red	 
    dir=${dir%_*}              			## Remove end with third number
    dir=${dir%_*}              			## Remove end with second number 
    #echo "$filename $dir"
    if [[ -d $dir ]]; then     			## If the directory exists
        mv "$filename" "$dir"/ 			## Move file there.
    elif  [[ ! -e $dir ]]; then 		## if doesn't exist
    	mkdir $dir						## make the directory
    	mv "$filename" "$dir"/ 			## Move file there.
    fi
done

for f in *txt  
do                        			 	## loop through all txt files 
    filename=${f##*/}                   ## Use the last part of a path.
    if [[ $filename == GSM* ]]; then    ## check to see if its a series file 
        extension=${f##*.}              ## Remove up to the last dot.
        dir=${filename%.*}              ## remove extension 
        #echo "$filename $dir"
        if [[ -d $dir ]]; then          ## If the directory exists
            mv "$filename" "$dir"/      ## Move file there.
        elif  [[ ! -e $dir ]]; then     ## if directory doesn't exist
            mkdir $dir                  ## make directory 
            mv "$filename" "$dir"/      ## Move file there.
        fi
    fi
done

for f in */                             ## for each subdirectory 
do
	mono /home/tuk32868/bin/autoconvert/AutoConvert.exe $f $f /mnt/d/visualization_pipeline/illumina_files/$bpm /mnt/d/visualization_pipeline/illumina_files/$egt  ## run autoconvert software
done


### make a list of all files in the directories 
find . -type f | while read name; do [[ $name == *.gtc ]] && echo $name; done > filelist.txt  

### call gtc2vcf using file just used 
bcftools +gtc2vcf --no-version -Ob -b /mnt/d/visualization_pipeline/illumina_files/$bpm -c /mnt/d/visualization_pipeline/illumina_files/$csv -e /mnt/d/visualization_pipeline/illumina_files/$egt -g filelist.txt -f $ref -x $output.sex -v --do-not-check-bpm | bcftools sort -Ou -T ./bcftools-sort.XXXX | bcftools norm --no-version -Ob -o $output.bcf -c x -f $ref && bcftools index -f $output.bcf

