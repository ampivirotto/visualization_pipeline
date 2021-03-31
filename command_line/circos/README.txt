Installation: 
	Install page on Circos - http://circos.ca/documentation/tutorials/configuration/installation/
	Commands:
		mkdir software
		mkdir software/circos
		cd software/circos
		# download Circos and place the archive in the directory
		wget http://circos.ca/distribution/circos-current.tgz
		# unpack
		tar xvfz circos-current.tgz
		# make a symlink to current
		ln -s circos-0.69-9/ current
		# delete the tarball, if you want
		# to install GD and perl modules
		sudo apt-get -y install libgd2-xpm-dev
		
	Perl packages:
		cd software/circos/circos-0.69-9/
		bin/circos -modules

		To install from the command line "cpan [library/module]"
		
		Carp
		Config::General
		Data::Dumper
		Digest::MD5
		File::Basename
		File::Spec::Functions
		FindBin
		GD
		GD::Polyline
		Getopt::Long
		Graphics::ColorObject
		IO::File
		List::MoreUtils
		List::Util
		Math::Bezier
		Math::BigFloat
		Math::Round
		Math::VecStat
		Memoize
		POSIX
		Params::Validate
		Pod::Usage
		Readonly
		Regexp::Common
		Set::IntSpan
		Storable
		Time::HiRes
	
	If there's issues installing GD with cpan, 
		apt-get install libgd-dev first and then use perl -MCPAN -e shell


Pipeline code:  
	1. Identify karyotype file and karyotype id from karyotype file saved in pickle file 
	2. Check the VCF format to see if the chromosomes are in # format or chr# format.  If in chr# format, remove the chr and save to NAME_edit.vcf
		a. awk '{gsub(/^chr/,""); print}' your.vcf > no_chr.vcf
		b. https://www.biostars.org/p/98582/
	3. Pull out the SNP density by writing a command line to .sh file and then run using subprocess.  
		a. awk '/^#/ {next} {printf("%s\t%d\n",$1,$2-$2%10000);}' input.vcf | sort | uniq -c | awk '{printf("hs%s\t%s\t%d\t%s\n",$2,$3,$3+10000,$1);}' > vcf.dat
		b. https://www.biostars.org/p/254227/
	4. Read in SNP density file to use bins in other .dat files 
	5. Read in the VCF using allel from scikit.  Set positions, chromosome numbers, and genotypes. 
	6. Using these values, find fst values using fst function.  If there's only one population, this will output an empty fst file. Otherwise, it will output a .dat file using makeDATFile
	7. Using allel get counts of heterozygotes and output into DAT file 
	8. Output circos conf file using outfn, directory, and the karyotype file  

Tutorials Used:
	https://dbsloan.github.io/TS2019/exercises/circos.html
	https://pbgworks.org/sites/pbgworks.org/files/Introduction%20to%20Circos.pdf
	http://circos.ca/documentation/tutorials/recipes/microbial_genomes/configuration
	(Potential R implementation) https://www.royfrancis.com/beautiful-circos-plots-in-r/
