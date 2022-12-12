args = R.utils::commandArgs(trailingOnly=TRUE)

print(args[1])
if(is.na(args[2])){
    print('None')
}