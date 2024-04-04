# Merge PROPHECY Variants

Consensus variants for the PROPHECY study were provided in two files:

1. Autosomes & scaffolds
2. chrX & chrY

These need to be merged & for current purposes, no variants on chrY are required

``` bash
cd /data/tki_bodl_rnaseq/references/create_variant_sets/vcf
micromamba activate bcftools
bcftools merge \
  --output prophecy_unrelated_PASS_MAF0.5_sites-only_merged.vcf.gz \
  --output-type 'z' \
  prophecy_unrelated_PASS_MAF0.5_sites-only.vcf.gz \
  prophecy_unrelated_PASS_MAF0.5_sites-only.chrX-chrY.vcf.gz
bcftools index -t prophecy_unrelated_PASS_MAF0.5_sites-only_merged.vcf.gz
```

After this step, the file information was collected using

``` bash
echo -e "Filename\tNbr Variants\tFile Size \tmd5sum" > fileinfo.txt
find *vcf.gz \
	-type f \
	-exec bash \
	-c 'md=$(md5sum "$0"); siz=$(wc -c <"$0"); ln=$(zcat <"$0" | egrep -vc "^#"); echo ${md} ${ln} $((${siz}/1048576))MB' {} \; |\
	awk -v OFS='\t' '{print $2,$3,$4,$1}' >> fileinfo.txt
```