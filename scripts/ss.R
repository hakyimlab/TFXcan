suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--command"),
        help='The command to run'),
    make_option("--num1", default=4, 
        help = "First number to sum or minus"),
    make_option("--num2", default=5, 
        help="Second number to sum or minus")
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

summ <- function(a, b){
    print(c(a, b))
    return(a + b)
}

minuss <- function(a, b){
    print(c(a, b))
    return(a - b)
}


main <- function(which_command){
    run <- switch(as.character(which_command),
        'summ' = summ(opt$num1, opt$num2),
        stop('The command does not exist.')
    )

    return(run)
}

main(opt$command)



# summ()
# Rscript /beagle3/haky/users/temi/projects/TFXcan/scripts/ss.R $a $b
# Rscript -e "summ()" /beagle3/haky/users/temi/projects/TFXcan/scripts/ss.R --num1 6 --num2 12 
# Rscript /beagle3/haky/users/temi/projects/TFXcan/scripts/ss.R --num1 6 --num2 12 --command "summ"
# Rscript /beagle3/haky/users/temi/projects/TFXcan/scripts/ss.R -e "summ()" $a $b 

# this works
# Rscript -e "source('/beagle3/haky/users/temi/projects/TFXcan/scripts/ss.R'); summ()"

# Rscript -e "source('/beagle3/haky/users/temi/projects/TFXcan/scripts/ss.R'); summ()" $a $b

# Rscript -e "system(paste())"
# but i need a variation of this to work
# a=4
# b=5
# Rscript -e "args <- commandArgs(trailingOnly = TRUE); a <- args[1]; b <- args[2]" -e "source('/beagle3/haky/users/temi/projects/TFXcan/scripts/ss.R'); summ(a, b)" ${a} ${b}