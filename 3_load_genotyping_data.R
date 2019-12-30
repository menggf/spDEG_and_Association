## This script is to load the imputed genotype data
## The files are stored in ./data/gwas directory

input=list()
input[["chr1"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr1.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr2"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr2.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr3"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr3.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr4"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr4.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr5"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr5.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr6"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr6.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr7"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr7.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr8"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr8.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr9"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr9.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr10"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr10.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr11"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr11.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr12"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr12.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr13"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr13.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr14"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr14.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr15"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr15.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr16"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr16.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr17"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr17.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr18"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr18.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr19"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr19.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr20"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr20.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr21"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr21.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")
input[["chr22"]]=load.gwaa.data(phenofile = "data/gwas/gt.ann", genofile = "data/gwas/chr22.dat", force = TRUE, makemap = FALSE, sort = TRUE, id = "TID")


for(i in 1:22){
   print(i)
   zz=paste("chr",i,sep="");
   input2=input[[zz]];
   snpnames=input2@gtdata@snpnames
   input2=input2[,snpnames[!duplicated(snpnames)]]
   qc1 <- check.marker(input2, p.level=0.00001,maf=0.1)
   input2 <- input2[qc1$idok, qc1$snpok]
   input[[zz]]=input2;
}
