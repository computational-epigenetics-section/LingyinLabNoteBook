time=Sys.time()
time=paste(substring(time,9,10),'_',substring(time,12,13),'-',substring(time,15,16),sep = '')


#Arguments:corsiv_file_name:csv
Args='all_CorRSIVs_06_27_2018.csv'
#Read Corsiv Data in CSV format
corsiv_file_name<-Args
corsiv<-read.csv(corsiv_file_name)

#get corsiv position chr*_*****
#position<-as.numeric(substring(corsiv$Bin.Name,regexpr('_',corsiv$Bin.Name)+1))
#get corsiv chromosome chr*_*****
#chr<-as.numeric(substring(corsiv$Bin.Name,4,regexpr('_',corsiv$Bin.Name)-1))

#bind chromosome, corsiv id, sample avmethylation level
corsiv1<-corsiv[,c(2,1,5:14)]
#average corsiv methylation level by corsiv id
csaverage<-aggregate(.~Uniq_ID, data=corsiv1, FUN=mean)

#bind corsiv id, coriv location
#clusterid<-corsiv[,1]
#corsiv2<-cbind(clusterid,position)

#get corsiv,chr,left start right end
#csleft<-aggregate(position ~ clusterid, data = corsiv2, FUN=function(x) mean(x)-99)
#csright<-aggregate(position ~ clusterid, data = corsiv2, max)
#left<-csleft[,2]
#right<-csright[,2]
left<-substring(corsiv[,37],regexpr(':',corsiv[,37])+1,regexpr('-',corsiv[,37])-1)
right<-substring(corsiv[,37],regexpr('-',corsiv[,37])+1)
corsiv3<-data.frame(corsiv$Uniq_ID,left,right,corsiv$X...2.5k.bp.Genes,corsiv$Overlapping.Genes)
corsiv3<-corsiv3[!duplicated(corsiv3),]
colnames(corsiv3)<-c('id','left','right','gene2.5kb','gene.overlap')
csavmethy<-merge(csaverage,corsiv3, by.x='Uniq_ID',by.y='id')
colnames(csavmethy)[1]<-'id'

#save for plotting use: id,sample corsiv methy, left right 25 gene, overlap gene
write.table(csavmethy,sprintf('corsiv_general_info_%s.txt',time),row.names=FALSE)
print('finish reading corsiv file')

#Write Methylation Level File into TXT format
#id samplenamess
avmethy<-csavmethy[,c(1,3:12)]
avmethy_file=sprintf('corsiv_methylation_%s.txt',time)
write.table(avmethy,avmethy_file,sep ='\t',row.names = FALSE)
avmethy_file_name=paste('./',avmethy_file,sep='')
print('finish writing methylation level file')

#get corsiv location data
#geneid chr left right
corsivloc<-csavmethy[,c(1,2,13,14)]
corsivloc_file=sprintf('corsiv_location_%s.txt',time)
write.table(corsivloc,corsivloc_file,sep='\t',row.names = FALSE)
corsiv_location_file_name=paste('./',corsivloc_file,sep='')
print('finish writing corsiv location file')
