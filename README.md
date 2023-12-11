# Create Sets of Variants from 1000 Genomes Porject Data

This documents the workflow used for creating the variant sets for populations within the 1000 Genomes datasets.
All files have been sourced from `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV` with the README in that folder containing a description of how files were created.
No Y-chromosome variants are available in this resource and as such, no variants on the Y-chromosome were included.

The complete dataset contains 3202 individuals, of which 2504 were considered to be unrelated.


All SNV and InDels <50bp with an allele frequency > 50% in the respective population will be retained.
These variant sets are being created to supply `STAR` with variants for use with the `STARconsensus` method, as well as for using `transmogR` to create a reference transcriptome for working with `salmon`.
In order to be viable for both InDels need to be restricted to <50bp.
This is not a strict requirement for `salmon` but is required for `STAR`

The sets of variants being created currently include

1. *Panhuman* variants across all populations  
	a. Using the 2504 unrelated individuals only  
  	b. Using all 3202 individuals
2. Variants from participants identified as being from the *African* (AFR) population
   a. Using unrelated individuals within the AFR population  
   b. Using all individuals within the AFR population

It should also be noted that according the [Kaminow et al]() paper, the AFR set of variants contains more unique variants, not included in the panhuman reference set, and as such, this will make an important additional set of variation.