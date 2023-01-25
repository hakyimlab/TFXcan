

sink(glue("{project_dir}/test_run/my_first_sink.txt"))

#prints numbers from 1 to 20
for (i in 1:20)
print(i)
sink()
#unlink(glue("{project_dir}/test_run/my_first_sink.txt")) # deletes the file

closeAllConnections()

print(2 + 2)


text_to_write <- "#!/bin/bash\n
#PBS -l select=1:system=polaris\n
#PBS -l walltime=00:59:00,filesystems=home:grand\n
#PBS -A covid-ct\n
#PBS -q preemptable\n 
#PBS -N scan_genomewide_motifs\n
#PBS -k doe\n
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/logs/scan_genomewide_motifs.out\n
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/logs/scan_genomewide_motifs.err"

fileConn <- file(glue("{project_dir}/test_run/my_first_sink.txt"))
writeLines(text_to_write, fileConn)
close(fileConn)