library(BiocParallel)
library(edgeR)
library(Rcpp)
library(fpc)
source("util.R")


## readin the HBRRC micorray expression data
da<-read.table("data/AMP-AD_HBTRC_MSSM_Agilent44Karray_PFC_AgeCorrected_all.tsv",header=T)
mx=as.matrix(da[,c(-1:-10)])
tag=!is.na(as.vector(da$gene_symbol))
mx=mx[tag,]
row.names(mx)<-as.vector(da$gene_symbol)[tag]
sd=apply(mx,1,sd,na.rm = T)  
ee=mx[sd>0.05,] # remove the genes with low expression variance

## readin the HBRRC clinical information data
clin=read.table("data/AMP-AD_HBTRC_MSSM_Agilent44Karray_Covariates.txt",header=T,sep="\t")
row.names(clin)<-paste("X",as.vector(clin$TID),sep="")
ad=row.names(subset(clin, DiseaseStatus=="Alzheimer"))
ct=row.names(subset(clin, DiseaseStatus=="Control"))

## AD and normal expression data
ee.ad=ee[,colnames(ee) %in% ad]
ee.ct=ee[,colnames(ee) %in% ct]

####################### single patient DEGs ###############################

library(BiocParallel)
## estimat the mean and sd indicate the expression profiles 
## in normal samples
mu = apply(ee.ct, 1, mean)
sd = apply(ee.ct, 1, sd,na.rm = T)

## estimate the single patient DEGs
pre.deg.hbtrc = sapply(colnames(ee.ad), function(x) {
    w = ee.ad[, x]
    z = (w - mu)/sd
    z[is.na(z)]=0
    return(z)
})
names(pre.deg.hbtrc) <- colnames(ee.ad)   
deg.hbtrc=apply(pre.deg.hbtrc, 2,function(x){
    cutoff=0.05
    p = 2*pnorm(-abs(x))
    nn=length(p[p < cutoff])
    if(nn > 2000){
      cutoff=sort(p)[2000]
    }
    bi=rep(0,length(x))
    bi[p < cutoff & x > 0 ] = 1
    bi[p < cutoff & x < 0 ] = -1
    return(bi)
})
colnames(deg.hbtrc) <-colnames(pre.deg.hbtrc)
row.names(deg.hbtrc) <- row.names(ee.ad)

## estimate the patient modules
library(DEComplexDisease)
res.deg=deg.specific(deg.hbtrc, min.genes=50, min.patients=5, cores=40,overlap=0.75)

## generated patient moules of different size 
pas.list=list()
allpas=colnames(ee.ad)

# module size == 60
for(x in ad[ad %in% names(res.deg)]){
    sig.ges=res.deg[[x]][["genes"]]
    dd=apply(ee.ad[sig.ges,], 2, function(y) sum((ee.ad[sig.ges,x]-y)**2, na.rm=T))
    pas=names(sort(dd)[1:60])
    pas.list[["pdeg60"]][[x]][["pas1"]]= pas;
    pas.list[["pdeg60"]][[x]][["pas2"]]= allpas[!allpas %in% pas];
    pas.list[["pdeg60"]][[x]][["pas3"]]= ct;    
}

# module size == 40
for(x in ad[ad %in% names(res.deg)]){
    sig.ges=res.deg[[x]][["genes"]]
    dd=apply(ee.ad[sig.ges,], 2, function(y) sum((ee.ad[sig.ges,x]-y)**2, na.rm=T))
    pas=names(sort(dd)[1:40])
    pas.list[["pdeg40"]][[x]][["pas1"]]= pas;
    pas.list[["pdeg40"]][[x]][["pas2"]]= allpas[!allpas %in% pas];
    pas.list[["pdeg40"]][[x]][["pas3"]]= ct;
}

# module size == 80
for(x in ad[ad %in% names(res.deg)]){
    sig.ges=res.deg[[x]][["genes"]]
    dd=apply(ee.ad[sig.ges,], 2, function(y) sum((ee.ad[sig.ges,x]-y)**2, na.rm=T))
    pas=names(sort(dd)[1:80])
    pas.list[["pdeg80"]][[x]][["pas1"]]= pas;
    pas.list[["pdeg80"]][[x]][["pas2"]]= allpas[!allpas %in% pas];
    pas.list[["pdeg80"]][[x]][["pas3"]]= ct;
}

# module size == 100
for(x in ad[ad %in% names(res.deg)]){
    sig.ges=res.deg[[x]][["genes"]]
    dd=apply(ee.ad[sig.ges,], 2, function(y) sum((ee.ad[sig.ges,x]-y)**2, na.rm=T))
    pas=names(sort(dd)[1:100])
    pas.list[["pdeg100"]][[x]][["pas1"]]= pas;
    pas.list[["pdeg100"]][[x]][["pas2"]]= allpas[!allpas %in% pas];
    pas.list[["pdeg100"]][[x]][["pas3"]]= ct;
}

# module size == 120
for(x in ad[ad %in% names(res.deg)]){
    sig.ges=res.deg[[x]][["genes"]]
    dd=apply(ee.ad[sig.ges,], 2, function(y) sum((ee.ad[sig.ges,x]-y)**2, na.rm=T))
    pas=names(sort(dd)[1:120])
    pas.list[["pdeg120"]][[x]][["pas1"]]= pas;
    pas.list[["pdeg120"]][[x]][["pas2"]]= allpas[!allpas %in% pas];
    pas.list[["pdeg120"]][[x]][["pas3"]]= ct;
}

ee.list=list()
allpas=colnames(ee.ad)
for(x in ad[ad %in% names(res.deg)]){
    sig.ges=res.deg[[x]][["genes"]]
    ee.list[[x]][["sig"]]=ee.ad[sig.ges,x]    
}


load("clin.rda") # load patient clinical information

#####################################################################################
#					Association studies												#
#																					#
#####################################################################################

## perform association studies using all the AD patients

temp= mclapply(1:22, function(i){
   zz=paste("chr",i,sep="");
   input2=input[[zz]];
   snpnames=input2@gtdata@snpnames
   smps=row.names(phdata(input2))
   grp=rep(-1, length(smps))
   grp[paste("X",smps,sep="") %in% ad ]=1
   grp[paste("X",smps,sep="") %in% ct ]=0
   phdata(input2)$groups=grp
   used.smps=smps[paste("X",smps,sep="") %in% ad | paste("X",smps,sep="") %in% ct] 
   input2.ccf1 <- qtscore( groups ~ pcs[,1]+pcs[,2]+pcs[,3]+pcs[,4]+pcs[,5], 
						  data=input2[used.smps, snpnames] ,quiet=T) 
   # association studies using imputation
   input2.ccf11 <- qtscore(groups ~ pcs[,1]+pcs[,2]+pcs[,3]+pcs[,4]+pcs[,5], 
						   data=input2[used.smps, snpnames], times =1000,quiet=T)   
   lb1=lambda(input2.ccf1)
   preres1=descriptives.scan(input2.ccf1, top=50000)
   res1=subset(preres1, P1df< 1e-3)
   snp=row.names(res1)
   res1=cbind(id=snp, res1, lambda=lb1$estimate, se=lb1$se);
   return(res1)
}, mc.cores=22, mc.preschedule=F)
library(dplyr)
all.gwas=bind_rows(temp)
have.snp=unique(as.vector(all.gwas$id))

library(biomaRt)
snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
ann.snp=getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "ensembl_gene_stable_id","ensembl_type","consequence_type_tv"),filters = c("snp_filter"),values = have.snp,mart = snpmart)
results0=cbind(all.gwas, ann.snp[have.snp, c("gene","chr_name","chrom_start","consequence_type_tv")])


## perform association studies using AD patient modules
a=vector()
b=vector()
for(tp in names(pas.list)){
  for(pa in names(pas.list[[tp]])){
    a=append(a, tp)
    b=append(b, pa)
  }
}
lst=data.frame(id=b,lab=a);
nn=nrow(lst)
gwas.hbtrc=mclapply(1:nn, function(y){
  pat.id =as.vector(as.matrix(lst[y,"id"]));
  tp=as.vector(as.matrix(lst[y,"lab"]));
     
  grp=rep(-1, length(smps))
  pas1=pas.list[[tp]][[pat.id]][["pas1"]];
  pas2=pas.list[[tp]][[pat.id]][["pas2"]];
  pas3=pas.list[[tp]][[pat.id]][["pas3"]];
  pas1=gsub("X","", pas1)
  pas2=gsub("X","", pas2)
  pas3=gsub("X","", pas3)
  rr1=data.frame(id=vector(),Chromosome=vector(), Position=vector(), 
					Strand=vector(), A1=vector(),A2=vector(), N=vector(), 
					effB=vector(),se_effB=vector(), chi2.1df=vector(),
					P1df=vector(),Pc1df =vector(), effAB=vector(), effBB=vector(),
					chi2.2df=vector(),P2df=vector(),lambda=vector(), se=vector())
  rr2=data.frame(id=vector(),Chromosome=vector(), Position=vector(), 
					Strand=vector(), A1=vector(),A2=vector(), N=vector(), 
					effB=vector(),se_effB=vector(), chi2.1df=vector(),
					P1df=vector(),Pc1df =vector(), effAB=vector(), 
					effBB=vector(),chi2.2df=vector(),P2df=vector(),
					lambda=vector(), se=vector())
  rr3=data.frame(id=vector(),Chromosome=vector(), Position=vector(),
					Strand=vector(), A1=vector(),A2=vector(), N=vector(), effB=vector(),
					se_effB=vector(), chi2.1df=vector(),P1df=vector(),Pc1df =vector(), 
					effAB=vector(), effBB=vector(),chi2.2df=vector(),P2df=vector(),
					lambda=vector(), se=vector())
  for(i in 1:22){
    zz=paste("chr",i,sep="");
    input2=input[[zz]];
    snpnames=input2@gtdata@snpnames
    smps=row.names(phdata(input2))
    snpnames=snpnames[!duplicated(snpnames)];
    grp=rep(-1, length(smps)) #  module patient vs control
    grp[smps %in% pas1 ]=1
    grp[smps %in% pas2 ]=2
    grp[smps %in% pas3 ]=0
    phdata(input2)$groups=grp
    used.smps=smps[smps %in% pas1 | smps %in% pas3]
	pcs=pc[[zz]][used.smps,]
    input2.ccf1 <- qtscore(groups ~ pcs[,1]+pcs[,2]+pcs[,3]+pcs[,4]+pcs[,5] , 
							data=input2[used.smps, snpnames] ,quiet=T) 
    input2.ccf11 <- qtscore(groups ~ pcs[,1]+pcs[,2]+pcs[,3]+pcs[,4]+pcs[,5] , 
							data=input2[used.smps, snpnames], times =1000,quiet=T) 
    lb1=lambda(input2.ccf1)
    preres1=descriptives.scan(input2.ccf1, top=10000)
    res1=subset(preres1, P1df< cutoff)
    id1=row.names(res1)
    p.adj1=as.vector(input2.ccf11@results[id1,"P1df"])
    if(nrow(res1) !=0){
      res1=cbind(id=row.names(res1), res1,lambda=lb1$estimate, se=lb1$se, p.adj=p.adj1);
      rr1=rbind(rr1, res1)
    }
    grp=rep(-1, length(smps)) #  module patient vs non-module patients
    grp[smps %in% pas1 ]=1
    grp[smps %in% pas2 ]=0
    grp[smps %in% pas3 ]=2
    phdata(input2)$groups=grp
    used.smps=smps[smps %in% pas1 | smps %in% pas2]
	  pcs=pc[[zz]][used.smps,]
    input2.ccf2 <- qtscore(groups ~ pcs[,1]+pcs[,2]+pcs[,3]+pcs[,4]+pcs[,5] , 
							data=input2[used.smps, snpnames] ,quiet=T) 
    lb2=lambda(input2.ccf2)
    preres2=descriptives.scan(input2.ccf2, top=length(snpnames)-100)
    res2=preres2[id1,]
    if(nrow(res2) !=0){
      res2=cbind(id=row.names(res2), res2, lambda=lb2$estimate, se=lb2$se);
      rr2=rbind(rr2, res2)
    }
    grp=rep(-1, length(smps)) #  module non-module patient vs control
    grp[smps %in% pas1 ]=2
    grp[smps %in% pas2 ]=1
    grp[smps %in% pas3 ]=0
    phdata(input2)$groups=grp  
    used.smps=smps[smps %in% pas2 | smps %in% pas3]
	  pcs=pc[[zz]][used.smps,]
    input2.ccf3 <- qtscore(groups ~ pcs[,1]+pcs[,2]+pcs[,3]+pcs[,4]+pcs[,5] , 
							data=input2[used.smps, snpnames] ,quiet=T) 
    lb3=lambda(input2.ccf3)
    lambda3=lb3$estimate
    se3=lb3$se
    preres3=descriptives.scan(input2.ccf3, top=length(snpnames)-100)
    res3=preres3[id1,]
    if(nrow(res3) !=0){
      res3=cbind(id=row.names(res3), res3, lambda=lb3$estimate, se=lb3$se);
      rr3=rbind(rr3, res3)
    }
  }
  rr1=rr1[order(rr1$P1df, decreasing=F),]
  rr2=rr2[order(rr2$P1df, decreasing=F),]
  rr3=rr3[order(rr3$P1df, decreasing=F),]
  list(res1=rr1,res2=rr2, res3=rr3, pa=pat.id, tp=tp, pas1=paste(pas1, collapse=","), pas2=paste(pas2, collapse=","), pas3=paste(pas3, collapse=","))
}, mc.preschedule=F, mc.cores=23)


outcome=data.frame() # tranform list as data.frame
cutoff=0.05 # set the 
for(i in 1:length(gwas.hbtrc) ){
  if(class(gwas.hbtrc[[i]])=="try-error")
    next()
  sub.res1=subset(gwas.hbtrc[[i]][["res1"]],p.adj < cutoff)
  snp=as.vector(sub.res1$id);
  if(length(snp) == 0)
     next()
  sub.res2= gwas.hbtrc[[i]][["res2"]][snp,]
  sub.res3= gwas.hbtrc[[i]][["res3"]][snp,]
  op=as.data.frame(cbind(id=snp, pid=snp, effB=as.vector(sub.res1$effB), lambda=as.vector(sub.res1$lambda),se=as.vector(sub.res1$se), P1df1=as.vector(sub.res1$P1df),P1df2=as.vector(sub.res2$P1df),P1df3=as.vector(sub.res3$P1df), p.adj=as.vector(sub.res1$p.adj), pa=gwas.hbtrc[[i]][["pa"]], tp=gwas.hbtrc[[i]][["tp"]], n1=length(strsplit(gwas.hbtrc[[i]][["pas1"]],",")[[1]]),n2=length(strsplit(gwas.hbtrc[[i]][["pas2"]],",")[[1]]),n3=length(strsplit(gwas.hbtrc[[i]][["pas3"]],",")[[1]])))
  outcome=rbind(op, outcome)
}
outcome$lambda=as.numeric(as.vector(outcome$lambda))
new.outcome=subset(outcome,lambda < 1.1  )
have.snp=unique(as.vector(new.outcome$id))

library(biomaRt)
snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
ann.snp=getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "ensembl_gene_stable_id","ensembl_type","consequence_type_tv"),filters = c("snp_filter"),values = have.snp,mart = snpmart)
results=cbind(new.outcome, ann.snp[have.snp, c("gene","chr_name","chrom_start","consequence_type_tv")])
results=results[order(as.numeric(as.vector(results$P1df1))),]

