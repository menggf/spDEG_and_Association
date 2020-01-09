## dowload the files from AMP-AD project:
https://www.synapse.org/#!Synapse:syn2580853
##################################################

##Read in genotyping data by format transformation from 'bed' to 'dat'
system("plink --noweb --bfile ../AMP-AD_HBTRC_MSSM_IlluminaHumanHap650Y --recode --tab --out gt")
convert.snp.ped(pedfile="gt.ped",mapfile="gt.map", outfile="gt.dat")
input=load.gwaa.data(phenofile = "gt.ann", genofile = "gt.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")

## QC analysis
qc1 <- check.marker(input, p.level=0)
input1 <- input[qc1$idok , qc1$snpok]
qc2 <- check.marker(input1,  fdr=0.2)
input2 <- input1[qc2$idok, qc2$snpok]

data1.gkin <- ibs(input2[, sample(autosomal(input2), 100000)], weight="freq")
data1.dist <- as.dist(0.5-data1.gkin)
data1.mds <- cmdscale(data1.dist)
strc=as.data.frame(data1.mds)
names(strc)<-c("PC1","PC2")
strc=cbind(strc, clin[row.names(strc),])
library(ggplot2)
plot(data1.mds) 
km <- kmeans(data1.mds, centers=2, nstart=1000)

pa=row.names(input@phdata);
remove.pa=pa[!pa %in%row.names(input2@phdata)]
ge=input@gtdata@snpnames;
remove.ge=ge[!ge %in% input2@gtdata@snpnames]
write.table(remove.ge,"remove_snp.txt",row.names=F,col.names=F,quote=F)
write.table(remove.pa,"remove_patients.txt",row.names=F,col.names=F,quote=F)
