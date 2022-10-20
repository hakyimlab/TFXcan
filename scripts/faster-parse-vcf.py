

import os, sys, re
import subprocess

sample_names = sys.argv[1]
chr_names = list(range(1, 23)) + ['Y']
#print(chr_names)

vcf_folder = '/projects/covid-ct/imlab/data/GEUVADIS/vcf_files'
output_folder = '/projects/covid-ct/imlab/users/temi/projects/TFXcan/output'

# header: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT

with open(sample_names, 'r') as sn:
    for line in sn:
        for chrom in chr_names:
            # each line is an individual
            query = line.split()
            
            # the regions for this individual
            regions = query[1:]
            
            # filter for those regions in that chromosome >> a boolean list is return
            filt_bool = [True if str(chrom) == re.split(':', reg)[0] else False for reg in regions]
            
            if not any(filt_bool):
                print(f'\nchr{chrom} not present. Moving on === \n')
                continue
            
            print(f'\nchr{chrom} present. Extracting variant information === \n')
            
            chrom_regions = [x for x, y in zip(regions, filt_bool) if y == True]
            chrom_regions = ','.join(chrom_regions)
        
        # create the sh call
            call_ = str(f'~/miniconda3/bin/bcftools view -v snps -H -r {chrom_regions} -s {query[0]} {vcf_folder}/ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz > {output_folder}/{query[0]}_chr{chrom}_subsetted_vcf_file.vcf')
            
            #call_ = str(f'~/miniconda3/bin/bcftools view -H -r {chr_name}:{",".join(query[1:])} -s {query[0]} {vcf_folder}/ALL.chr{chr_name}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz > {output_folder}/{query[0]}_subsetted_vcf_file.vcf')
        
            #print(call_)
            
            # call the shell script
            subprocess.run(call_, shell=True)
            subprocess.run(str(f'gzip -f {output_folder}/{query[0]}_chr{chrom}_subsetted_vcf_file.vcf'), shell=True) # removed the txt file
        
        