library(MatrixEQTL)
time=Sys.time()
time=paste(substring(time,9,10),'_',substring(time,12,13),'-',substring(time,15,16),sep = '')


#methylation data, corsiv locatino data
Args=c('corsiv_methylation_02_15-38.txt',
       'corsiv_location_02_15-38.txt',
       'snp_genotype_25_10-42.txt',
       'snp_location_25_10-42.txt',
  sprintf('mqtl_result_%s.txt',time),
  1,
  1e5)
avmethy_file_name=Args[1]
corsiv_location_file_name=Args[2]
SNP_file_name=Args[3]
snps_location_file_name=Args[4]
# Output file name
output_file_name_cis = Args[5]

# Only associations significant at this level will be saved
pvOutputThreshold_cis = as.numeric(Args[6])
print(pvOutputThreshold_cis)


# Distance for local gene-SNP pairs
cisDist = as.numeric(Args[7])

#matrixEQTL Analysis
useModel = modelLINEAR

# Covariates file name
# Set to character() for no covariates
covariates_file_name = character()



## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(avmethy_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(corsiv_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  useModel = useModel, 
  errorCovariance = numeric(), 
  verbose = TRUE, 
  output_file_name='tranout',
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold = pvOutputThreshold_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = TRUE);

if (nrow(me$cis$eqtls) > 0) {
  assoc <- read.csv(output_file_name_cis, sep="\t", head=T)
  assoc$r2 <- (assoc$t.stat)^2 /(me$param$dfFull + assoc$t.stat^2)
  write.table(assoc, file=output_file_name_cis, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
} else{
  ## remove if no associations ...
  unlink(output_file_name)
}
