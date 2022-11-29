

def create_region_file(open_vcf, region, subset_vcf_dir, individual, software_paths=[]):
    '''
    Creates a subsetted vcf file per region
    Arguments:
        open_vcf: A vcf that is already opened
        region: region to query, a list [chr, start, end]
        subset_vcf_dir: the directory to save the subsetted file
        individual: the individual to subset the file for
        software_paths: a list of paths to software to use [bcftools, tabix]
    '''

    import kipoiseq
    import subprocess

    path_to_bcftools = software_paths[0]
    path_to_tabix = software_paths[1]

    SEQUENCE_LENGTH = 393216
    
    # Center the interval at the region
    interval = kipoiseq.Interval(region[0], region[1], region[2]).resize(SEQUENCE_LENGTH) # resizing will change the regions
    path = f'{subset_vcf_dir}/{individual}_{interval.chr}_{interval.start}_{interval.end}_subset_genotypes.vcf.gz'
    region_interval = f'{interval.chr}:{interval.start}-{interval.end}'
    view_cmd = f"{path_to_bcftools} view {open_vcf} -r {region_interval} -s {individual} --output-type z --output-file {path} && {path_to_tabix} -p vcf {path}"
    out = subprocess.run(view_cmd, shell=True)

    return {'subset_path':path, 'interval':interval, 'individual':individual, 'region':region}


def extract_individual_sequence(subset_dict, fasta_file_path, fasta_extractor, delete_region=False):

    '''
    Extracts a sequence from a reference for a region with the variants of an individual applied
    Arguments:
        subset_dict: dict - the result of `create_region_file()`
        fasta_file_path: string/path - the path to the fasta file
        fasta_extractor: extractor object - in case there are not variants to apply to that region, this helps extract the reference sequence instead
        delete_region: bool - after the sequence has been applied, should the temporary vcf file be deleted?
    '''
    import kipoiseq
    import warnings
    import os

    kseq_extractor = kipoiseq.extractors.SingleSeqVCFSeqExtractor(fasta_file=fasta_file_path, vcf_file=subset_dict['subset_path'])
    center = subset_dict['interval'].center() - subset_dict['interval'].start

    individual = subset_dict['individual']
    individual_sequences = {}
    for ind in [individual]:

        warnings.filterwarnings('error')
        
        try:
            individual_sequences[ind] = kseq_extractor.extract(interval=subset_dict['interval'], anchor=center, sample_id=ind)
            seq_source = 'var'
        except Warning:
            warnings.simplefilter("always", category=UserWarning)
            print('No variants for this region. Using reference genome.\n')
            individual_sequences[ind] = fasta_extractor.extract(interval=subset_dict['interval'], anchor=[])
            seq_source = 'ref'

    if delete_region == True:
        os.remove(subset_dict['subset_path'])
        os.remove(f"{subset_dict['subset_path']}.tbi")

    return {'sequence':individual_sequences, 'sequence_source':seq_source, 'region':subset_dict['region'][3]}


def create_sequences(sample, region, fasta_file_path, fasta_extractor, open_vcf_file, temporary_vcf_dir, software_paths, script_path):

    # create the region file
    a = create_region_file(open_vcf=open_vcf_file, region=region, subset_vcf_dir=temporary_vcf_dir, individual=sample, software_paths=software_paths)
    print(f'Region file is created')
    try:
        print('Extracting individual sequences')
        b = extract_individual_sequence(subset_dict=a, fasta_file_path=fasta_file_path, fasta_extractor=fasta_extractor, delete_region=True)
        #b['sequence'][sample] = one_hot_encode(b['sequence'][sample])[np.newaxis]
        return(b)
    except ValueError:
        #status = [0, 'NA'] #evaluate_check(is_query_done, b={'sequence_source':'NA'})
        return({'sequence':None, 'sequence_source':'NA', 'region':'NA'})

def create_input(sample, region, fasta_file_path, fasta_extractor, open_vcf_file, temporary_vcf_dir, software_paths):

    # create the region file
    a = create_region_file(open_vcf=open_vcf_file, region=region, subset_vcf_dir=temporary_vcf_dir, individual=sample, software_paths=software_paths)
    print(f'Region file is created')
    try:
        print('Extracting individual sequences')
        b = extract_individual_sequence(subset_dict=a, fasta_file_path=fasta_file_path, fasta_extractor=fasta_extractor, delete_region=True)
        b['sequence'][sample] = one_hot_encode(b['sequence'][sample])[np.newaxis]
        return(b)
    except ValueError:
        #status = [0, 'NA'] #evaluate_check(is_query_done, b={'sequence_source':'NA'})
        return({'sequence':None, 'sequence_source':'NA', 'region':'NA'})
