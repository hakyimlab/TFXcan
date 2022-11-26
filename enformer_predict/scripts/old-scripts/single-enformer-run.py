
import argparse, os, re, sys, h5py, csv, warnings, json
import tensorflow as tf
import numpy as np
import kipoiseq
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--sequence_region', help="Region to predict on", type=str, default='')
parser.add_argument('--sequence_sample', help="Name or id of the individual", type=str, default='')
args = parser.parse_args()

def main():

    # get the path of the script as well as parameters         
    whereis_script = os.path.dirname(sys.argv[0])  
    script_path = os.path.abspath(whereis_script) 

    usage_codes = f'{script_path}/enformer-usage-codes.py'
    personal_enformer = f'{script_path}/personal-enformer.py'

    # import the enformer-usage_codes.py file
    exec(open(usage_codes).read(), globals(), globals())
    exec(open(personal_enformer).read(), globals(), globals())

    # read the parameters file
    with open(f'{script_path}/../metadata/enformer_parameters.json') as f:

        parameters = json.load(f)

        model_path = parameters['model_path']
        fasta_file = parameters['hg38_fasta_file']
        output_dir = parameters['output_dir']
        vcf_file = parameters['vcf_file']
        path_to_bcftools = parameters['path_to_bcftools']
        path_to_tabix = parameters['path_to_tabix']
        temporary_vcf_dir = parameters['temporary_vcf_dir']
        TF = parameters['TF']
        logfile_path = parameters['logfile_path']
        sequence_folder = parameters['sequence_folder']

    #enf_model = tf.saved_model.load(model_path).model
    sequence_extractor = FastaStringExtractor(fasta_file)

    # load the model
    enf_model = tf.saved_model.load(model_path).model

    # I need a log file
    # read in the log file for this individual ; doing this so that the log file is not opened everytime
    logfile_csv = f'{logfile_path}/{args.sequence_sample}_predictions_log.csv'
    if os.path.isfile(logfile_csv):
        open_mode = 'a'
    else:
        open_mode = 'w'

    with open(logfile_csv, open_mode, encoding='UTF8') as running_log_file:
        logwriter = csv.writer(running_log_file)
        if open_mode == 'w':
            logwriter.writerow(['motif', 'individual', 'status', 'sequence_type'])

        # extract sequence
        qsplit = args.sequence_region.split('_')
        region = [qsplit[0], int(qsplit[1]), int(qsplit[2]), args.sequence_region]
        input_dict = create_input(sample=args.sequence_sample, region=region, fasta_file_path=fasta_file, fasta_extractor=sequence_extractor, open_vcf_file=vcf_file, temporary_vcf_dir=temporary_vcf_dir, software_paths=[path_to_bcftools, path_to_tabix])

        # with open(args.sequence_file, 'r') as s:
        #     sequence = s.readline()

        motif_pred = run_predictions(sequence_encoded=input_dict['sequence'][args.sequence_sample], region=args.sequence_region, sample=args.sequence_sample, seq_type=input_dict['sequence_source'], model=enf_model, output_dir=output_dir)

        logwriter.writerow(motif_pred)

        # always flush
        running_log_file.flush()
        os.fsync(running_log_file)


if __name__ == '__main__':
    main()
