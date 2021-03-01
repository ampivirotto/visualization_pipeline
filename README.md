# visualization_pipeline

Documentation:

	->Requirements:
		flask
		scikit-allel
		GEOparse
		wget
		pickle
		numpy 
		matplotlib
		seaborn
		bcolz
		pandas
		
	
	->GUI Website:
		In the src/ folder we have the instrastructure for a flask webpage.  The subdirectory venv contains the environment variables, the subdirectory app contains the material for the front facing webseries, and the subdirectory pipeline contains the python programs to run with the webpage.  

		To run the website: 
			Activate the virtual environment: 
				venv\Scripts\activate
			Tell Flask where to find the webpage file:
				export FLASK_APP=webpage_test.py
				or for windows set FLASK_APP=webpage_test.py
				OR you can permanently set it by adding FLASK_APP=webpage_test.py to .flaskenv
			Run: 
				On the command line run: flask run 
			Open up the local page in a web browser window.  
			
		Input Fields
			GEO Accession Number:  the series accession number of the study of interest which can be found at Gene Expression Omnibus (GEO)
			Output Label:  output label for resulting files (vcf, bcf, tab) and image files (for any figures)
			ChipType: Illumina Beadchip used in study.  This could be found on the GEO page for the study under platform. If you encounter errors in the pipeline, it may be that the version is incorrect.  
			Visualization: Types of visualization graphs to output. These will be saves as jpgs in the series folder.  
			Subsample: If you would like to just select certain samples from the series.  Once all the files are read in then the user can select the samples of interest.  
			
		TO DO:
			- allow user to select column in metadata
			- allow user to select subset of data 
		
	->Command Line Program:
		Input Fields:
		
		TO DO: 
			- modify web pipeline program to run at the command line 
		
	
	->Pipeline:
		1. Check if series (GSE) has been downloaded yet and create series directory.  
		2. Identify the files needed based on Beadchip type.  
		3. Download all the sample IDAT files and download metadata.
		4. Create bash file and run bash file which will:
			a. Move IDAT files into SAMPLE directories. 
			b. Move metadata into SAMPLE directories.  
			c. Run Autoconvert on all IDAT files > converts to GTC format  
			d. Create filelist of all GTC files (there will be one per sample)
			e. Convert GTC files to BCF format using BCFtools GTC2VCF
			f. Sort, Normalize, and Index BCF file. 
			g. Create VCF format from BCF format.  
		
		TO DO for Pipeline: 
			- compare a single sample to the other samples 
		
	->Types of Visualization: 
		Heterozygosity - Graphs the heterozygosity (percent of 0/1 or 1/0 per variant) across the genome
		PCA - Principal Components Analysis 
		Circos - circular graph that can contain SNP density, heterozygosity, and FST
		Site Frequency Spectrum - Graph of the number of variants for each derived allele count
		Heat Topology - 
		
		TO DO: 
			- allow user to select specific tracts for the circos plots
			- more forms of visualization 
	
	-> References: 
		- https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-i-hello-world
		- Scikit-allel
		- https://github.com/freeseek/gtc2vcf
	
	
	-> Possible Tutorials to Explore: 
		-http://corearray.sourceforge.net/tutorials/SNPRelate/
		-https://nbviewer.jupyter.org/github/rasbt/pattern_classification/blob/master/clustering/hierarchical/clust_complete_linkage.ipynb
		-https://sebastianraschka.com/Articles/heatmaps_in_r.html
		-PCA http://hardingnj.github.io/2018/08/15/making-interactive-pca-plots.html
		-Scikit-allel http://alimanfoo.github.io/2016/06/10/scikit-allel-tour.html
		-https://plantgenomics.univie.ac.at/radseq/snpfiltering/ heatmap2 R function 
		-http://people.virginia.edu/~wc9c/KING/manual.html
		-Saving Plots in R https://www.stat.berkeley.edu/~s133/saving.html

				
 
