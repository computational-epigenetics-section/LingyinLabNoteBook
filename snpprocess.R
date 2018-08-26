library(vcfR)

time=Sys.time()
time=paste(substring(time,9,10),'_',substring(time,12,13),'-',substring(time,15,16),sep = '')

#Arguments:genotype_files:txt containing all sample genotype files, tab deliminated
Args='genotype.txt'

#get genotype data
#read txt containing list of sample genotype file
geno_files<-read.table(Args,sep='\t',header=FALSE)
#read the first sample file
vcf<-read.vcfR(as.character(geno_files[,1]))
print(sprintf('read first vcf file %s',geno_files[,1]))
vcfdf<-cbind(as.data.frame(getFIX(vcf)),as.data.frame(vcf@gt))
print('convert vcf to dataframe')

#initialize id and genotype vector
id<-vector()
genotype<-vector()
chr<-vcfdf$CHROM
pos<-vcfdf$POS
ref<-vcfdf$REF
alt<-vcfdf$ALT
sample_name<-as.character(colnames(vcfdf)[9])

#get id and genotype from first vcf file
id<-paste(chr,'_',pos,sep='')
print(vcfdf)
genotype<-as.vector(as.numeric(substr(vcfdf[9],1,1))+as.numeric(substr(vcfdf[9],3,3)))

print('start extracting genotype and location')
#construct data frame of snp genotype
snpgt<-data.frame(id,genotype)
colnames(snpgt)[2]<-sample_name
#construct data frame of snp location
snploc<-data.frame(id,chr,pos)

#extract the 9th(genotype) column of other sample gt file
for(sample_file in geno_files[2:length(geno_files)]){
  sample_vcf<-read.vcfR(as.character(sample_file))
  sample_df<-cbind(as.data.frame(getFIX(sample_vcf)),as.data.frame(sample_vcf@gt))
  sample_name<-as.character(colnames(sample_df)[9])
  genotype<-as.vector(as.numeric(substr(sample_df[,9],1,1))+as.numeric(substr(sample_df[,9],3,3)))
  snpgt<-cbind(snpgt,genotype)
  colnames(snpgt)[length(snpgt)]<-sample_name
}



#eliminate rows with NA genotype
snp<-cbind(ref,alt,snploc,snpgt)
org<-length(snp$id)
snp<-snp[complete.cases(snp), ]
omNA<-length(snp$id)
print(sprintf('After removing missing values: %s out of %s snps remained',omNA,org))

#eliminate tag snps
snp<-snp[snp$ref %in% c('A','C','T','G')&snp$alt %in% c('A','C','T','G'),]
print(sprintf('After removing non-SNPS:%s out of %s snps remained',length(snp$id),omNA))

#save for plotting
write.table(snp,'snp_ref_alt.txt',sep='\t',row.names=FALSE)
#get refined data
snploc<-as.data.frame(snp[,3:5])
snpgt<-snp[,-c(1,2,3,4,5)]



#write snp data into TXT format
#id samplegenotype
snp_genotype_file_name=sprintf('snp_genotype_%s.txt',time)
write.table(snpgt,snp_genotype_file_name,sep ='\t',row.names = FALSE)

#write snplocation data
snps_location_file_name=sprintf('snp_location_%s.txt',time)
write.table(snploc,'snploc.txt',sep='\t',row.names=FALSE)
