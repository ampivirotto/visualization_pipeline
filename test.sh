dir='/mnt/d/visuzalization_pipeline/data/GSE93106'
ogdir=''

#cd dir
#gunzip *.gz

#for f in *.idat           # no need to use ls.
#do
    #filename=${f##*/}          # Use the last part of a path.
    #extension=${f##*.}         # Remove up to the last dot.
    #dir=${filename%_*}         # Remove "tv" in front of filename.
    #dir=${dir%_*}              # Remove episode
    #dir=${dir%_*}              # Remove season
    #echo "$filename $dir"
    #if [[ -d $dir ]]; then     # If the directory exists
        #mv "$filename" "$dir"/ # Move file there.
    #elif  [[ ! -e $dir ]]; then
    	#mkdir $dir
    	#mv "$filename" "$dir"/ # Move file there.
    #fi
#done

for f in */
do
	#echo $f
	mono /home/tuk32868/bin/autoconvert/AutoConvert.exe /mnt/d/visualization_pipeline/data/GSE93106/$f /mnt/d/visualization_pipeline/data/GSE93106/$f /mnt/d/visualization_pipeline/data/GSE93106/HumanOmniExpress-24-v1-0-B.bpm /mnt/d/visualization_pipeline/data/GSE93106/HumanOmniExpress_24v1-0_A.egt
done
