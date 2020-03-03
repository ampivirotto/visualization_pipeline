bpm_manifest_file= "/mnt/d/visualization_pipeline/data/canine_snps/CanineHD_B.bpm"
csv_manifest_file= "/mnt/d/visualization_pipeline/data/canine_snps/CanineHD_B.csv"
egt_cluster_file="/mnt/d/visualization_pipeline/data/canine_snps/CanineHD_A.egt"
gtc_list_file="/mnt/d/visualization_pipeline/python_programs/gtclist.txt"
ref="$HOME/res/ncbi-genomes-2020-03-03/GCF_000002285.3_CanFam3.1_genomic.fna" # or ref="$HOME/res/human_g1k_v37.fasta"
out_prefix="canine_snps"
bcftools +gtc2vcf \
  --no-version -Ou \
  -b $bpm_manifest_file \
  -c $csv_manifest_file \
  -e $egt_cluster_file \
  -g $gtc_list_file \
  -f $ref \
  -x $out_prefix.sex | \
  bcftools sort -Ou -T ./bcftools-sort.XXXXXX | \
  bcftools norm --no-version -Ob -o $out_prefix.bcf -c x -f $ref && \
  bcftools index -f $out_prefix.bcf