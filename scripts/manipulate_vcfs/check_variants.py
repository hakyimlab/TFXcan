import os
import numpy as np
import cyvcf2
import kipoiseq
import pandas as pd

working_dir = '/grand/covid-ct/imlab/users/temi/projects/TFXcan/scripts/'
os.chdir(working_dir)

reference_genome = '/grand/covid-ct/imlab/data/hg_sequences/hg38/Homo_sapiens_assembly38.fasta'
vcf_files_dir = '/grand/covid-ct/imlab/data/GEUVADIS/vcf_snps_only'

class FastaStringExtractor:
    def __init__(self, fasta_file):
        import pyfaidx

        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.

        import kipoiseq
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = kipoiseq.Interval(interval.chrom,
                                    max(interval.start, 0),
                                    min(interval.end, chromosome_length),
                                    )
        # pyfaidx wants a 1-based interval
        sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                            trimmed_interval.start + 1,
                                            trimmed_interval.stop).seq).upper()
        # Fill truncated values with N's.
        pad_upstream = 'N' * max(-interval.start, 0)
        pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()

chromosomes = [f'chr{i}' for i in list(range(1, 23))]
chromosomes.extend(['chrX', 'chrY', 'chrM'])

fasta_extractor = FastaStringExtractor(reference_genome)


variants_status = {}

for chr in chromosomes:
    vcf_f = f'{vcf_files_dir}/ALL.{chr}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz'
    if os.path.isfile(vcf_f):
        chr_variants = []
        for i, variant in enumerate(cyvcf2.cyvcf2.VCF(vcf_f)):
            out = []
            out.extend([variant.CHROM, variant.start, variant.end, variant.REF, variant.ALT])
            chr_variants.append(out)
    
    chr_reference_variants_positions = [[v[0], v[1], v[2]] for v in chr_variants]
    print(chr_reference_variants_positions[1:3])

    chr_reference_variants = []
    for vpos in chr_reference_variants_positions:
        slength = len(range(vpos[1], vpos[2]))
        reg_interval = kipoiseq.Interval(vpos[0], vpos[1], vpos[2]).resize(slength)
        reference_region = fasta_extractor.extract(interval=reg_interval, anchor=[])
        chr_reference_variants.append(reference_region)

    chr_reference_variants = ''.join(chr_reference_variants)
    chr_vcf_ref_variants = ''.join([v[3] for v in chr_variants])

    status = chr_reference_variants == chr_vcf_ref_variants
    #print(f'    {status}')

    variants_status[chr] = status



