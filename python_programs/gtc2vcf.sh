bpm_manifest_file="..."
csv_manifest_file="..."
egt_cluster_file="..."
gtc_list_file="..."
ref="$HOME/res/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" # or ref="$HOME/res/human_g1k_v37.fasta"
out_prefix="..."
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