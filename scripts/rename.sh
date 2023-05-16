


bcftools annotate --rename-chrs ${impute_one}/misc/chr_annotations.txt --threads 10 -Oz -o ${impute_one}/resource/${name_without_chr}_chrRenamed.vcf.gz ${one_vcf}