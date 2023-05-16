

import vcf
import kipoiseq
from kipoiseq import Interval
from pyliftover import LiftOver
import os, sys
import gzip

lo_1 = LiftOver('/project2/haky/temi/temp_globus/hackathon/hg19ToHg38.over.chain.gz')
lo = LiftOver('/project2/haky/temi/temp_globus/hackathon/hg38ToHg19.over.chain.gz')

TF = 'Gata3'
samples = ['HG00096']
chroms =  [str(l) for l in range(1, 23)]
chroms.append('X')
sample_count = 2
# This loop constructs "individual" vcf files split by gene window of enformer input size

# chroms = [str(x) for x in range(1,23)] + ["X"]
# genes_of_interest =  ["ERAP1"]
for chrom in chroms:
    print('Starting with chr{x} ===\n'.format(x=chrom))
    with open("/project2/haky/temi/projects/TFXcan/processed-data/" + TF + "_motif_regions.txt", "r") as motif_regions:

        vcf_reader = vcf.Reader(filename='/project2/haky/temi/temp_globus/hackathon/ALL.chr' +chrom+ '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz')

        for line in motif_regions:
            p_line = line.strip().split()
            chr_name = p_line[0]
            if chr_name != ('chr'+chrom):
                continue
            
            set_name = str(p_line[6] + p_line[7]) # either TP1 or TN25355

            # check if the bed file exists
            if os.path.isfile("/project2/haky/temi/projects/TFXcan/processed-data/motif-bed/chr" +chrom+ "/chr_" + chrom + "_motif_" + set_name +".bed"):
                print('{sname} bed file exists.\n'.format(sname=set_name))
                continue

            print('{sname} bed file does not exist. Creating one...\n'.format(sname=set_name))
            vcf_writer = vcf.Writer(open("/project2/haky/temi/projects/TFXcan/processed-data/motif-vcf/chr" + chrom + "_motif_" + set_name +".vcf","w"), vcf_reader)       # open output vcf file

            # target_interval = kipoiseq.Interval("chr"+chrom,                # Construct interval coordinates
            #                             int(p_line[1]), # Start 
            #                             int(p_line[2])) # End
            # print(target_interval)
            #window_coords = target_interval.resize(SEQUENCE_LENGTH)         # Modify the length to match intended window size; not needed 
            window_coords_start = p_line[1]
            window_coords_end = p_line[2]
            print('{a}, {b}, {c}, {d}\n'.format(a=set_name, b=chrom, c=window_coords_start, d=window_coords_end))
            if len(lo_1.convert_coordinate("chr"+chrom, int(window_coords_start))) == 0:        # Liftover coordinates to GRCh38
                continue
            if len(lo_1.convert_coordinate("chr"+chrom, int(window_coords_end))) == 0:
                continue
            hg38_start = int(lo_1.convert_coordinate("chr"+chrom, int(window_coords_start))[0][1])
            hg38_end = int(lo_1.convert_coordinate("chr"+chrom, int(window_coords_end))[0][1])

            if hg38_end < hg38_start:
                continue
            # if os.path.exists("/project2/haky/temi/projects/TFXcan/processed-data/motif-vcf/chr" + chrom + "_motif_" + set_name +".vcf"):
            #     print("{sname} done already".format(set_name))
            #     continue
            for record in vcf_reader.fetch(chrom, hg38_start, hg38_end):  
                vcf_writer.write_record(record)
            vcf_writer.close()

            os.system("mkdir -p "+"/project2/haky/temi/projects/TFXcan/processed-data/motif-bed/chr"+chrom)

            with open("/project2/haky/temi/projects/TFXcan/processed-data/motif-vcf/chr" + chrom + "_motif_" + set_name +".vcf","r") as i:                         # Convert vcf files into simplified bed files
                with open("/project2/haky/temi/projects/TFXcan/processed-data/motif-bed/chr" +chrom+ "/chr_" + chrom + "_motif_" + set_name +".bed", "w") as o:
                    
                    for line in i:
                        try:
                            if line.startswith("##"):
                                # print(line)
                                continue

                            p_line = line.strip().split("\t")

                            if line.startswith("#"):
                                o_line = [
                                    p_line[0],
                                    p_line[1],
                                    p_line[3],
                                    p_line[4]
                                    ]

                                sample_bool = []
                                for ind,samp in enumerate(p_line[9:]):
                                    if samp in samples:
                                        sample_bool.append(True)
                                        o_line.append(samp)
                                    else:
                                        sample_bool.append(False)

                                o.write("\t".join(o_line)+"\n")

                            # print(p_line[3])
                            if len(p_line[3]) != 1:
                                continue
                            if len(p_line[4]) != 1:
                                continue

                            all_ref = True
                            for geno in p_line[len(p_line)-sample_count: len(p_line)]:
                                if geno == "0|0":
                                    continue
                                else:
                                    all_ref = False
                                    break
                            if all_ref:
                                continue

                            # print(p_line[0])
                            # print(p_line[1])
                            # print(lo.convert_coordinate("chr"+p_line[0], int(p_line[1])))
                            # print(lo.convert_coordinate("chr"+p_line[0], int(p_line[1]))[0][1])
                            if len(lo.convert_coordinate("chr"+p_line[0], int(p_line[1]))[0]) == 0:
                                continue
                            o_line = [
                                p_line[0],
                                str(lo.convert_coordinate("chr"+p_line[0], int(p_line[1]))[0][1]),
                                p_line[3],
                                p_line[4]
                                ]

                            for ind,samp in enumerate(p_line[9:]):
                                if sample_bool[ind]:
                                    o_line.append(samp)


                            # print(o_line)
                            # print("\t".join(o_line)+"\n")
                            o.write("\t".join(o_line)+"\n")
                        except:
                            continue
            
            # remove the vcf files - too large
            os.remove('/project2/haky/temi/projects/TFXcan/processed-data/motif-vcf/chr{x}_motif_{y}.vcf'.format(x=chrom, y=set_name))
            print('/project2/haky/temi/projects/TFXcan/processed-data/motif-vcf/chr{x}_motif_{y}.vcf has been removed.\n'.format(x=chrom, y=set_name))
