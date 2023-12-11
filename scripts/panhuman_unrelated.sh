#!/bin/bash
#SBATCH -p tki_bodl
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem=16GB
#SBATCH -o /home/tki_stephenp/TKI/create_variant_sets/logs/%x_%j.out
#SBATCH -e /home/tki_stephenp/TKI/create_variant_sets/logs/%x_%j.err

## This script follows this basic outline
##
## 1. Create a header from an existing file
##    - Remove the FORMAT field as genotypes will be exculded
##    - Add column names as a separate step, inly including the 1st 8 columns
## 2. Step through each autosome
##    - Filter using bcftools to an AF > 0.5, excluding headers & genotypes
##    - Use awk to remove InDels > 50bp
##    - Use egrep to remove SVs
## 3. Add information from chrX
##    - Filter using bcftools to an AF > 0.5, excluding headers & genotypes
##    - Use awk to remove InDels > 50bp
##    - Use egrep to remove SVs
## 4. Compress the VCF
## 5. Index the VCF

micromamba activate bcftools

WD="/data/tki_bodl/pederson/create_variant_sets"
OUTFILE="${WD}/vcf/1000GP_SNV_INDEL_panhuman_unrel_0.5.vcf"
VCFPATH="/data/tki_bodl/data/references/1000GP_NYGC_30x_GRCh38_Phased/"
VCFPRE="${VCFPATH}/1kGP_high_coverage_Illumina"
VCFPOST="filtered.SNV_INDEL_SV_phased_panel"
SUF="vcf.gz"


## Create the header
echo -ne "$(date)\t Creating header for ${OUTFILE}..."
bcftools view -h "${VCFPRE}.chr22.${VCFPOST}.${SUF}" |\
  egrep '^##' |\
  egrep -v "^##FORMAT" > ${OUTFILE}
## Add the column names, excluding the genotypes
zcat "${VCFPRE}.chr22.${VCFPOST}.${SUF}" |\
  egrep -m1 "^#CHROM" |
  sed -r 's/(.+INFO).+/\1/g' >> ${OUTFILE}
echo -e "done"


## Now step through each autosome, filtering correctly, removing genotypes &
## headers, then using awk to remove long InDels and egrep to remove structural 
## variants
for ((i=1; i<=22; i++))
do
	echo -ne "$(date)\tFiltering variants from chr${i}..."
  zcat "${VCFPRE}.chr${i}.${VCFPOST}.${SUF}" |\
	  bcftools view \
      -G \ 
      -H \
      -i '(INFO/AC_EUR_unrel + INFO/AC_EAS_unrel + INFO/AC_AMR_unrel + INFO/AC_SAS_unrel + AC_AFR_unrel) > 0.5*(INFO/AN_EUR_unrel + INFO/AN_EAS_unrel + INFO/AN_AMR_unrel + INFO/AN_SAS_unrel + AN_AFR_unrel)' |\
	 awk 'length($4) <= 50 && length($5) <= 50' |\
	 egrep -v 'HGSV' >> ${OUTFILE}
   echo -e "done"
done

## Adding chrX will require removing the Hemi fields first, then operating the same way
EXCL='INFO/AC_Hemi_EAS,INFO/AC_Hemi_AMR,INFO/AC_Hemi_EUR,INFO/AC_Hemi_AFR,INFO/AC_Hemi_SAS,INFO/AC_Hemi_EUR_unrel,INFO/AC_Hemi_EAS_unrel,INFO/AC_Hemi_AMR_unrel,INFO/AC_Hemi_SAS_unrel,INFO/AC_Hemi_AFR_unrel,INFO/AC_Hemi'
## Note that this file has v2 in it's name. 
## It is also worth noting that in this file, male individuals are changed to 1 genotypes in non-PAR regions if the allele was detected.
## Using 0.5 * the sum of allele numbers is prudent as these were changed as well.
## In the PAR regions, the allele counts will be the same (2504)
echo -ne "$(date)\tFiltering variants from chrX..."
bcftools annotate -x ${EXCL} "${VCFPRE}.chrX.${VCFPOST}.v2.${SUF}" |\
  bcftools view \
    -G \
    -H \
    -i '(INFO/AC_EUR_unrel + INFO/AC_EAS_unrel + INFO/AC_AMR_unrel + INFO/AC_SAS_unrel + AC_AFR_unrel) > 0.5*(INFO/AN_EUR_unrel + INFO/AN_EAS_unrel + INFO/AN_AMR_unrel + INFO/AN_SAS_unrel + AN_AFR_unrel)' |\
  awk 'length($4) <= 50 && length($5) <= 50' |\
  egrep -v 'HGSV' >> ${OUTFILE}
  echo -e "done"

  echo -ne "$(date)\tCompressing VCF..."
  bgzip ${OUTFILE}
  echo -e "done"
  echo -ne "$(date)\tIndexing VCF..."
  tabix -p vcf ${OUTFILE}.gz
  echo -e "done"
  echo -e "$(date)\tAll steps complete"