


Yes. I want to train a tree-based ML model that predicts gene expression of a gene based on binding of a TF at the promoter of all protein-coding genes. I also want to interpret this model to construct a gene regulatory network that connects TF binding at promoters to expression. 


The predictors are predicted binding of a TF e.g. GATA1, at the promoter of all 19,000 protein coding genes. The terms or observations are in each individual. So, a matrix of ~400 individual X 19,000 TSS of genes. The outcome variable is a vector of gene expression for a gene in ~400 individuals. 