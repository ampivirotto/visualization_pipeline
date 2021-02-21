
ref='/content/reference/canFam3.fa'
#export BCFTOOLS_PLUGINS="/home/tuk32868/bin"  ## make sure plugins are set 

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
	/content/iaap-cli/./iaap-cli gencall $direct/$f $direct/$f /content/illumina_files/$bpm /content/illumina_files/$egt /content/data  ## run iaap-cli software
done


### make a list of all files in the directories 
find . -type f | while read name; do [[ $name == *.gtc ]] && echo $name; done > filelist.txt  

### call gtc2vcf using file just used 
bcftools +gtc2vcf --no-version -Ob -b /content/illumina_files/$bpm -c /content/illumina_files/$csv -e /content/illumina_files/$egt -g $direct/filelist.txt -f $ref -x $direct/$output.sex -v --do-not-check-bpm | bcftools sort -Ou -T ./bcftools-sort.XXXX | bcftools norm --no-version -Ob -o $direct/$output.bcf -c x -f $ref && bcftools index -f $direct/$output.bcf

