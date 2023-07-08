


from parsl.app.app import python_app


def create_locus_info(chr_, loci_list):

    vcf_path = str('/lus/grand/projects/TFXcan/imlab/data/GEUVADIS/vcf_snps_only/ALL.{}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz')

    queries = []
    for l in loci_list:
        locus_split = l.split('_')
        qchr = locus_split[0]
        qstart = locus_split[1]
        qend = locus_split[2]
        query = f'{qchr}:{qstart}-{qend}'
        queries.append(query)

    return({'vcf_file': vcf_path.replace('{}', chr_), 'queries': queries})

@python_app
def create_prediction_matrix(chr_, chr_info, individuals, tfpred_matrix, save_dir):
    import re, os
    import pandas as pd
    import cyvcf2

    if not os.path.isfile(chr_info['vcf_file']):
        raise Exception(f"ERROR - {chr_info['vcf_file']} does not exist.")

    def find_variants_in_vcf_file(cyvcf2_object, samples, queries, offset=200000):

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

        import re
        variants_dictionary = {}

        for l in queries:
            query_dictionary = {}
            locus_split = re.split(':|-', l)
            qchr = locus_split[0]
            qstart = int(locus_split[1]) - offset
            qend = int(locus_split[2]) + offset
            query = f'{qchr}:{qstart}-{qend}'

            query_dictionary['chr'] = qchr
            query_dictionary['positions'] = tuple(variant.POS for variant in cyvcf2_object(query))
            query_dictionary['ref'] = tuple(variant.REF for variant in cyvcf2_object(query)) 
            query_dictionary['alt'] = tuple(variant.ALT for variant in cyvcf2_object(query)) 

            if not query_dictionary['positions']:
                print(f'{query} does not exist in vcf file')
                continue
            else:
                for i, sample in enumerate(samples):
                    try:
                        if sample in cyvcf2_object.samples:
                            sample_variants = tuple([variant.genotypes[i][0:2], variant.gt_bases[i].split('|')] for variant in cyvcf2_object(query))
                            # return(sample_variants)
                            sample_alleles = [''.join(sample_variants[i][1]) for i in range(len(sample_variants))]
                            query_dictionary[sample] = sample_alleles
                    except UserWarning:
                        print(f'[WARNING] {sample} is not in the VCF file.')
                        continue
            try:
                vv = pd.DataFrame(query_dictionary)
                colnames = ['chr', 'position', 'ref', 'alt']
                colnames.extend(samples)
                vv.columns = colnames
                vv["locus"] = vv["chr"] + ':' + vv["position"].astype(str) 
                vv = vv.drop(['chr', 'position'], axis=1)
                column_to_move = vv.pop("locus")
                vv.insert(0, "locus", column_to_move)
                variants_dictionary[l] = vv
            except ValueError:
                variants_dictionary[l] = None

        return(variants_dictionary)
    
    def alleles_to_dosages(aa):
        import pandas as pd

        ref = aa[1]
        alt = aa[2]
        bb = pd.Series(aa[3:])
        update_dict = {f'{ref}{a}':1 for a in alt}
        update_dict.update({f'{a}{ref}':1 for a in alt})
        update_dict.update({f'{a*2}':2 for a in alt})
        update_dict.update({f'{ref*2}':0})
        bb = bb.map(update_dict)
        return(bb)
    
    def collect_locus_tfpred_score(locus, tfpred_dt):
        import re
        nlocus = re.sub(':|-', '_', locus)
        return(tfpred_dt.loc[tfpred_dt['locus'] == nlocus].set_index('locus').T.reset_index())

    def create_per_locus_training_matrix(dosages_dt, locus, tfpred_dt):
        locus_tfpred = collect_locus_tfpred_score(locus=locus, tfpred_dt = tfpred_dt)
        result = locus_tfpred.merge(dosages_dt, how='inner', left_on='index', right_on='index')
        return(result)
    
    # read the vcf file + find variants
    chr_vcf_cy = cyvcf2.cyvcf2.VCF(chr_info['vcf_file'], samples=individuals)
    locus_genotypes = find_variants_in_vcf_file(chr_vcf_cy, individuals, chr_info['queries'])
    chr_vcf_cy.close()

    info1 = {k: v.apply(alleles_to_dosages, axis=1).set_index(v['locus']).T.reset_index() for k, v in locus_genotypes.items() if v is not None}
    info2 = {k: create_per_locus_training_matrix(v, k, tfpred_matrix) for k, v in info1.items() if v is not None}

    # saving the data
    chrF = os.path.join(save_dir, chr_)
    if not os.path.isdir(chrF): os.makedirs(chrF, exist_ok = True)
    save_alleles = [v.to_csv(os.path.join(chrF, f'{k}_1KG_alleles.csv'), sep = '\t', index=False) for k, v in locus_genotypes.items()]
    save_dosages = [v.to_csv(os.path.join(chrF, f'{k}_1KG_dosages.csv'), sep = '\t', index=False) for k, v in info1.items()]
    save_prediction_matrices = [v.to_csv(os.path.join(chrF, f'{k}_1KG_prediction_matrix.csv'), sep = '\t', index=False) for k, v in info2.items()]

    print(f"INFO - Allele information for {chr_} have been collected and saved.")

    return(0)