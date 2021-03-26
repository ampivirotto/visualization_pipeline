
#ref='/content/reference/canFam2.fa'
ref='/content/visualization_pipeline/command_line/igv.js-flask/igvjs/static/data/public/canFam2.fa'
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
    
    python /content/visualization_pipeline/command_line/change_idat_names2.py $(pwd | grep -P 'data\/\w+' -o | cut -d'/' -f 2)
	/content/iaap-cli/./iaap-cli gencall /content/illumina_files/$bpm /content/illumina_files/$egt ./ -f ./ -g  ## run iaap-cli software
done


### make a list of all files in the directories 
find . -type f | while read name; do [[ $name == *.gtc ]] && echo $name; done > filelist.txt  

### call gtc2vcf using file just used 
bcftools +gtc2vcf --no-version -Ob -b /content/illumina_files/$bpm -c /content/illumina_files/$csv -e /content/illumina_files/$egt -g $direct/filelist.txt -f $ref -x $direct/$output.sex -v --do-not-check-bpm | bcftools sort -Ou -T ./bcftools-sort.XXXX | bcftools norm --no-version -Ob -o $direct/$output.bcf -c x -f $ref && bcftools index -f $direct/$output.bcf
### convert bcf to vcf
bcftools view *.bcf -Ov --output $direct/$output.vcf && printf "Converted bcf to vcf\n"
### bgzip the vcf
bgzip -c *.vcf > $direct/$output.vcf.gz && printf "Bgzipped the vcf\n"
### index bgzipped vcf with tabix
tabix -p vcf *vcf.gz && printf "Indexed the bgzipped vcf with tabix\n"

printf "Done ( ͡° ͜ʖ ͡°)"
#last line doesn't get printed with the subprocess version to start the pipeline, so this line makes sure the real last line gets printed
printf ""
