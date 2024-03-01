# rnaseq_gseapy
This repository utilizes RNAseq differential expression analysis input to perform GSEA in python. The python package used for this pipeline development is [gseapy](https://github.com/zqfang/GSEApy/tree/master)

## Basic Installation ##
Install gseapy package from bioconda or pip. Follow installation instructions from [README](https://github.com/zqfang/GSEApy)

## How to run the pipeline ##

1. Copy rnaseq_gseapy.params from [parameterFile](rnaseq_gseapy.params) to your current folder
2. Make necessary changes in the parameter file. Give path for de_file_edgeR NOTE : you can use other DE tool such as DEseq2 output too. Just assign correct column name in the parameter file. 
3. Execute
```
    python check_libraries_gseapy.py  #This command will save all the libraries from gseapy in a file name library_names_gseapy.csv
    python rnaseq_gseapy.py rnaseq_gseapy.params #This will mainly run enrichr overrepresentatin and ranked gsea analysis
    
```
