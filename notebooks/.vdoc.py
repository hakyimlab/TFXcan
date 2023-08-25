# type: ignore
# flake8: noqa
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
vcf_file = '/lus/grand/projects/TFXcan/imlab/data/GEUVADIS/vcf_snps_only/ALL.chr4.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz'
=======
```{r}
glue('conda activate compbio-tools\nbcftools view -H -r {paste0(snps_dt$bcf_r, collapse=\',\')} -s {paste0(ind_names, collapse=",")} /lus/grand/projects/covid-ct/imlab/data/GEUVADIS/vcf_snps_only/ALL.chr4.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz')
#
#
#
vcf_file = '/lus/grand/projects/covid-ct/imlab/data/GEUVADIS/vcf_snps_only/ALL.chr4.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz'
>>>>>>> 8d1e50b0eb3dd8bdfc1db9105ea1984253c9762d
#
#
#
import cyvcf2
import pandas as pd
import os
import numpy as np
#
#
#
os.getcwd()
#
#
#
<<<<<<< HEAD
test_variants = pd.read_table('/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/TFPred_pipeline/analysis/test_variants.txt')
=======
test_variants = pd.read_table('../analysis/test_variants.txt')
>>>>>>> 8d1e50b0eb3dd8bdfc1db9105ea1984253c9762d
test_variants
#
#
#
<<<<<<< HEAD
individuals = pd.read_table('/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/TFPred_pipeline/metadata/valid_subset1KGgenomes_individuals.txt', header=None)[0].to_list()
=======
individuals = pd.read_table('/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/metadata/valid_subset1KGgenomes_individuals.txt', header=None)[0].to_list()
>>>>>>> 8d1e50b0eb3dd8bdfc1db9105ea1984253c9762d
individuals[0:5], len(individuals)
#
#
#
vcf_chr = cyvcf2.cyvcf2.VCF(vcf_file, samples=individuals)
#
#
#
#
def find_variants_in_vcf_file(cyvcf2_object, samples, queries):

    """
    Given a cyvcf2 object and a kipoiseq.Interval object, as well as a list of samples, extract the variants for the samples for the intervals.
    Parameters:
        cyvcf2_object: A cyvcf2 object
        interval_object: a kiposeq.Interval object
            should have `chrom`, `start`, and `end` attributes
        samples: list
            a list of samples: [a, b, c]
            At least one of these samples should be in the vcf file. 
    Returns: dict
        chr: the chromosome
        position: the list of positions
        sample: the samples and variants information
    """

    #n_samples = len(samples)

    # check that samples are in the vcf file
    # if not set(samples).issubset(cyvcf2_object.samples):
    #     raise Exception(f'[ERROR] Fatal. Some samples are not in the VCF file.')

    variants_dictionary = {}

    for query in queries:
        query_dictionary = {}
        query_dictionary['chr'] = query.split(':')[0]
        query_dictionary['positions'] = tuple(variant.POS for variant in cyvcf2_object(query))

        if not query_dictionary['positions']:
            print(f'{query} does not exist in vcf file')
            continue
        else:
            for i, sample in enumerate(samples):
                try:
                    if sample in cyvcf2_object.samples:
                        sample_variants = tuple([variant.genotypes[i][0:2], variant.gt_bases[i].split('|')] for variant in cyvcf2_object(query))
                        # return(sample_variants)
                        sample_alleles = ['_'.join(sample_variants[i][1]) for i in range(len(sample_variants))]
                        query_dictionary[sample] = sample_alleles
                except UserWarning:
                    print(f'[WARNING] {sample} is not in the VCF file.')
                    continue
        variants_dictionary[query] = pd.DataFrame(query_dictionary)

    output = pd.DataFrame(np.vstack(list(variants_dictionary.values())))
    colnames = ['chr', 'position']
    colnames.extend(samples)
    output.columns = colnames
    return(output)
#
#
#
#
queries = [f"{test_variants['chr'].tolist()[i]}:{test_variants['start'].tolist()[i]}-{test_variants['end'].tolist()[i]}" for i in range(test_variants.shape[0])]
queries
```
#
individual_snps = find_variants_in_vcf_file(vcf_chr, samples=individuals, queries=queries)
#
#
#
output_snps = individual_snps.merge(test_variants.loc[:,['chr', 'start', 'snp_id']], left_on=['chr', 'position'], right_on=['chr', 'start'], how='inner')
output_snps.drop('start', axis=1, inplace=True)
output_snps = output_snps[['chr', 'position', 'snp_id'] + [c for c in output_snps if c not in ['chr', 'position', 'snp_id']]]
output_snps.head()
#
#
#
output_snps.to_csv(f'/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/analysis/geuvadis_99_variants.txt', sep='\t', index=False)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
