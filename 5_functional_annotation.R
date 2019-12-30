
gwas.gene=unique(as.vector(results$gene)) # AD risk gene
gwas.gene=gwas.gene[!is.na(gwas.gene) & gwas.gene!="-" & gwas.gene!=""]

library(clusterProfiler)
library(org.Hs.eg.db)
go <- enrichGO(gene         = gwas.gene,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.1)
                
pdf("go_gwas.pdf",height=5,width=7)
plot(dotplot(go))
dev.off()

sel.snp=unique(as.vector(results$pid))
rel.snp=sapply(sel.snp, function(z){
    if(length(grep(":",z)) == 0 )
      return(z)
    zz=strsplit(z,":")[[1]]
    if(length(grep("rs\\d+",z))!=0)
      return(zz[1])
    cr=as.numeric(zz[1])
    loc=as.numeric(zz[2])
    sub.snp.ann = snp.list[[cr]]
    snp.id=as.vector(sub.snp.ann$V3)[which.min(abs(as.vector(sub.snp.ann$V2) - loc))]
    return(snp.id)
  })
names(rel.snp) <- sel.snp
smps=row.names(phdata(input[["chr1"]]))
smps=smps[smps!="15995" & smps!="16105" & smps!="21926" & smps!="22001"]
mx=matrix(nrow=0, ncol=length(smps))
colnames(mx) <- smps
mx=as.data.frame(mx)
for(i in 1:22){
  print(i)
 zz=paste("chr",i,sep="");
 input2=input[[zz]]; 
 sel=sel.snp[sel.snp %in% colnames(gtdata(input2))];
 if(length(sel)==0)
   next()
 sub.input2=input2[smps, sel]
 t(as.genotype(sub.input2@gtdata))->gt
 mx=rbind(mx, gt)
}

mx=as.matrix(mx)
temp=rel.snp[row.names(mx)]
row.names(mx)<- temp
colnames(mx)<-paste("X",colnames(mx),sep="" )

exp=cbind(ee.ad,ee.ct)
used.smps=intersect(colnames(exp), colnames(mx));
new.mx=mx[, used.smps]
new.mx2=apply(new.mx,c(1,2), function(x){
  if(is.na(x))
    return(NA)
  if(x == "1/1")
    return(0)
  if(x == "1/2" | x == "2/1")
    return(1)
  if(x == "2/2")
    return(2)  
})

new.exp=exp[,used.smps]
n.ges=nrow(new.exp);

used.snp=row.names(new.mx2)
preused.ges=row.names(new.exp)
da<-read.table("../pasgwas/ensembl_gene.txt",header=T,sep="\t")
row.names(da)<-as.vector(da$ensembl)
used.ges= as.vector(da[preused.ges, ]$gene)
#row.names(new.exp)<-used.ges
tag=grep("CHR", as.vector(da$chr))
genloc=unique(da[-1*tag,c(2,3,4,5)])


names(genloc)<-c("gene","chr_name","left","right")
genloc=genloc[!duplicated(genloc$gene),]
row.names(genloc)<-as.vector(genloc$gene)

snploc=unique(results[,c("id","chr_name","chrom_start")])
names(snploc)<-c("snpid","chr","pos")
row.names(snploc)<-as.vector(snploc$snpid)
snploc=snploc[!is.na(snploc$pos),]
library(MatrixEQTL)
et=Matrix_eQTL_main(
    snps=SlicedData$new()$CreateFromMatrix(new.mx2),
    gene=SlicedData$new()$CreateFromMatrix(new.exp),
    output_file_name = "myeqtl.txt",
    pvOutputThreshold = 1e-2,
    useModel = modelLINEAR,
    output_file_name.cis = "myeqtl_cis.txt",
    pvOutputThreshold.cis = 0.01,
    snpspos = snploc,
    genepos = genloc,
    cisDist = 1e6,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

cis=et$cis$eqtls
trans=subset(et$trans$eqtls, pvalue< 1e-5)

has=as.vector(results$id)  
ann.trans=sapply(has, function(x) paste(as.vector(subset(trans, snps==x,)$gene), collapse="," ))
ann.cis=sapply(has, function(x) paste(as.vector(subset(cis, snps==x,)$gene), collapse="," ))

final=cbind(results, trans=ann.trans, cis=ann.cis)

write.table(final, "final_hbtrc.txt", row.names=F, quote=F, sep="\t")
save(list=ls(), file="mid_hbtrc.rda")

###############################################################
cols=rainbow(12)
smps=colnames(mx)
for(i in 1:nrow(final)){
  snp=as.vector(as.matrix(final[i,"id"]))
  pa =as.vector(as.matrix(final[i,"pa"]))
  tp =as.vector(as.matrix(final[i,"tp"]))
  p =as.vector(as.matrix(final[i,"P1df1"]))
  pas1=pas.list[[tp]][[pa]][["pas1"]]
  pas2=pas.list[[tp]][[pa]][["pas2"]]
  pas3=pas.list[[tp]][[pa]][["pas3"]]
  pas1=pas1[pas1 %in% smps]
  pas2=pas2[pas2 %in% smps]
  pas3=pas3[pas3 %in% smps]
  
  labs=c(rep("Module",length(pas1)), rep("non-Module",length(pas2)), rep("Control",length(pas3)))
  gt=mx[snp, c(pas1,pas2,pas3)]
  gt[is.na(gt)]<-"1/1"
  tb=table(gt,labs)
  tb1=as.matrix(table(gt,labs))
  tb0=table(labs)
  tb2=tb1[,c(1,2)]
  tb3=tb2
  tb3[,2]=tb1[,2]+tb1[,3]
  colnames(tb3)<-c("Control","AD")
  p2=chisq.test(tb3*100)$p.value
  tb[,1]=tb[,1]/tb0[1]
  tb[,2]=tb[,2]/tb0[2]
  tb[,3]=tb[,3]/tb0[3]
  tb=as.matrix(tb)
  lbs=paste(c("Module","non-Module","Control"),"\n(n=",tb0[c("Module","non-Module","Control")],")",sep="")
  pdf(paste("pic_hbtrc/",snp,"_",pa,"_",tp,".pdf",sep=""),width=5,height=5)
  barplot(tb[,c("Module","non-Module","Control")], col=cols[c(4,8,2)],names.arg=lbs,main=paste("Module p=",signif(as.numeric(p),digits=3), sep=""))
  dev.off()
}

#######  generated a unique one ##############################


for(pa in names(ee.list)){
  temp=subset(final, pa==pa)
  ee.list[[pa]][["result"]]=temp
}
ee.list.hbtrc=ee.list;
save(ee.list.hbtrc, file="ee.list.hbtrc.rda")
save(ee.ad, input, ct, file="temp_hbtrc_eva.rda") 

#######################################################################
load("ee.list.rnaseq.rda")
nn=length(ee.list.rnaseq)

load("imp_gwas2.hbtrc.rda")
#load("snp.list.rda")
library(biomaRt)
outcome2=data.frame()

for(i in 1:length(gwas2.hbtrc) ){
  if(class(gwas2.hbtrc[[i]])=="try-error")
    next()
  sub.res1=subset(gwas2.hbtrc[[i]][["res1"]],P1df < 1e-4)
  snp=as.vector(sub.res1$id);
  if(length(snp) == 0)
     next()
  
  new.snp=sapply(snp, function(z){
    if(length(grep(":",z)) == 0 )
      return(z)
    zz=strsplit(z,":")[[1]]
    if(length(grep("rs\\d+",z))!=0)
      return(zz[1])
    cr=as.numeric(zz[1])
    loc=as.numeric(zz[2])
    sub.snp.ann = snp.list[[cr]]
    snp.id=as.vector(sub.snp.ann$V3)[which.min(abs(as.vector(sub.snp.ann$V2) - loc))]
    return(snp.id)
  })
   op=as.data.frame(cbind(id=new.snp, pid=snp, effB=as.vector(sub.res1$effB), lambda=as.vector(sub.res1$lambda),se=as.vector(sub.res1$se), P1df1=as.vector(sub.res1$P1df), pa=gwas2.hbtrc[[i]][["pa"]], tp=as.vector(sub.res1$lambda)))
  outcome2=rbind(op, outcome2)
}
  
have.snp=unique(as.vector(outcome2$id))
snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
test.snp=getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "ensembl_gene_stable_id","ensembl_type"),filters = c("snp_filter"),values = have.snp,mart = snpmart)
id <- read.table("../pasgwas/ensembl_gene.txt",header=T,sep="\t")
row.names(id)<-as.vector(id$ensembl)
test.snp=cbind(test.snp, gene=as.vector(id[as.vector(test.snp$ensembl_gene_stable_id),]$gene))
new.test.snp=data.frame()
for(snp in unique(as.vector(test.snp$refsnp_id))){
  tt=subset(test.snp, refsnp_id==snp)
  if(nrow(tt)==1){
    new.test.snp=rbind(tt,new.test.snp)
  }else{
    gg=unique(as.vector(tt$gene))
    gg=gg[!is.na(gg)]
    tt$gene=paste(gg,collapse=",")
    new.test.snp=rbind(tt[1,],new.test.snp)
  }
}
row.names(new.test.snp) <-as.vector(new.test.snp$refsnp_id);
has.snp=as.vector(outcome2$id)
results2=cbind(outcome2, new.test.snp[has.snp, c("gene","chr_name","chrom_start")])

###########################################################################################################
#check overlap

results.rnaseq=ee.list.rnaseq[[3]][["result"]]
results.hbtrc=ee.list.hbtrc[[3]][["result"]]

new.results.rnaseq0=data.frame()
for(chr in c(1:22)){
    pos1= as.vector(subset(results.rnaseq, chr_name==chr,)$chrom_start)
    temp=subset(results2 , chr_name==chr,);
    pos2= as.vector(temp$chrom_start)
    #print(length(pos2))
    dis=sapply(pos2, function(x) min(abs(x-pos1),  na.rm = TRUE)) < 100000
    labs=rep("No", length(dis))
    labs[dis]<-"Yes"
    new.results.rnaseq0=rbind( cbind(temp, findit=labs),new.results.rnaseq0)
}


new.results.rnaseq=data.frame()
for(chr in c(1:22)){
    pos1= as.vector(subset(results.hbtrc, chr_name==chr,)$chrom_start)
    temp=subset(results.rnaseq , chr_name==chr,);
    pos2= as.vector(temp$chrom_start)
    #print(length(pos2))
    dis=sapply(pos2, function(x) min(abs(x-pos1),  na.rm = TRUE)) < 100000
    labs=rep("No", length(dis))
    labs[dis]<-"Yes"
    new.results.rnaseq=rbind( cbind(temp, findit=labs),new.results.rnaseq)
}

new.results.hbtrc=data.frame()
for(chr in c(1:22)){
    pos1= as.vector(subset(results.rnaseq, chr_name==chr,)$chrom_start)
    temp=subset(results.hbtrc , chr_name==chr,);
    pos2= as.vector(temp$chrom_start)
    #print(length(pos2))
    dis=sapply(pos2, function(x) min(abs(x-pos1),  na.rm = TRUE)) < 100000
    labs=rep("No", length(dis))
    labs[dis]<-"Yes"
    new.results.hbtrc=rbind( cbind(temp, findit=labs),new.results.hbtrc)
}

write.table(new.results.hbtrc, "test.txt",row.names=F,quote=F,sep="\t")
write.table(new.results.rnaseq, "test2.txt",row.names=F,quote=F,sep="\t")


########### simulation evlaution

load("imp_gwas_hbtrc3_3.rda")
outcome.simu1=data.frame()
cutoff=0.05
for(i in 1:length(gwas.hbtrc) ){
  if(class(gwas.hbtrc[[i]])=="try-error")
    next()
  #sub.res1=subset(gwas.hbtrc[[i]][["res1"]],Pc1df < 1e-5)
  sub.res1=subset(gwas.hbtrc[[i]][["res1"]],as.numeric(p.adj) < 0.05)
  snp=as.vector(sub.res1$id);
  if(length(snp) == 0)
     next()
  new.snp=sapply(snp, function(z){
    if(length(grep(":",z)) == 0 )
      return(z)
    zz=strsplit(z,":")[[1]]
    if(length(grep("rs\\d+",z))!=0)
      return(zz[1])
    cr=as.numeric(zz[1])
    loc=as.numeric(zz[2])
    sub.snp.ann = snp.list[[cr]]
    snp.id=as.vector(sub.snp.ann$V3)[which.min(abs(as.vector(sub.snp.ann$V2) - loc))]
    return(snp.id)
  })
   op=as.data.frame(cbind(id=new.snp, pid=snp, effB=as.vector(sub.res1$effB), lambda=as.vector(sub.res1$lambda),se=as.vector(sub.res1$se), P1df1=as.vector(sub.res1$P1df), p.adj=as.vector(sub.res1$p.adj), pa=gwas.hbtrc[[i]][["pa"]], tp=gwas.hbtrc[[i]][["tp"]]))
  outcome.simu1=rbind(op, outcome.simu1)
}
outcome.simu1$lambda=as.numeric(as.vector(outcome.simu1$lambda))
new.outcome.simu1=subset(outcome.simu1,lambda < 1.1  )

load("imp_gwas_hbtrc4_0_1.rda")
outcome.simu2=data.frame()
cutoff=0.05
for(i in 1:length(gwas.hbtrc) ){
  if(class(gwas.hbtrc[[i]])=="try-error")
    next()
  #sub.res1=subset(gwas.hbtrc[[i]][["res1"]],Pc1df < 1e-5)
  sub.res1=subset(gwas.hbtrc[[i]][["res1"]],as.numeric(p.adj) < 0.05)
  sub.res2=gwas.hbtrc[[i]][["res2"]][row.names(sub.res1),]
  sub.res3=gwas.hbtrc[[i]][["res3"]][row.names(sub.res1),]
  snp=as.vector(sub.res1$id);
  if(length(snp) == 0)
     next()
  new.snp=sapply(snp, function(z){
    if(length(grep(":",z)) == 0 )
      return(z)
    zz=strsplit(z,":")[[1]]
    if(length(grep("rs\\d+",z))!=0)
      return(zz[1])
    cr=as.numeric(zz[1])
    loc=as.numeric(zz[2])
    sub.snp.ann = snp.list[[cr]]
    snp.id=as.vector(sub.snp.ann$V3)[which.min(abs(as.vector(sub.snp.ann$V2) - loc))]
    return(snp.id)
  })
   op=as.data.frame(cbind(id=new.snp, pid=snp, effB=as.vector(sub.res1$effB), lambda=as.vector(sub.res1$lambda),se=as.vector(sub.res1$se), P1df1=as.vector(sub.res1$P1df), P1df2=as.vector(sub.res2$P1df),P1df3=as.vector(sub.res3$P1df),p.adj=as.vector(sub.res1$p.adj), pa=gwas.hbtrc[[i]][["pa"]], tp=gwas.hbtrc[[i]][["tp"]]))
  outcome.simu2=rbind(op, outcome.simu2)
}
outcome.simu2$lambda=as.numeric(as.vector(outcome.simu2$lambda))
new.outcome.simu2=subset(outcome.simu2,lambda < 1.1  )

have.snp=unique(as.vector(new.outcome.simu2$id))

snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
test.snp=getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "ensembl_gene_stable_id","ensembl_type","consequence_type_tv"),filters = c("snp_filter"),values = have.snp,mart = snpmart)
id <- read.table("../pasgwas/ensembl_gene.txt",header=T,sep="\t")
row.names(id)<-as.vector(id$ensembl)
test.snp=cbind(test.snp, gene=as.vector(id[as.vector(test.snp$ensembl_gene_stable_id),]$gene))
new.test.snp=data.frame()
for(snp in unique(as.vector(test.snp$refsnp_id))){
  tt=subset(test.snp, refsnp_id==snp)
  if(nrow(tt)==1){
    new.test.snp=rbind(tt,new.test.snp)
  }else{
    gg=unique(as.vector(tt$gene))
    gg=gg[!is.na(gg)]
    tt$gene=paste(gg,collapse=",")
    new.test.snp=rbind(tt[1,],new.test.snp)
  }
}
row.names(new.test.snp) <-as.vector(new.test.snp$refsnp_id);
has.snp=as.vector(new.outcome.simu2$id)
results.simu2=cbind(new.outcome.simu2, new.test.snp[has.snp, c("gene","chr_name","chrom_start","consequence_type_tv")])
results.simu2$P1df1=as.numeric(as.vector(results.simu2$P1df1))
results.simu2=results.simu2[order(results.simu2$P1df1),]

outcome2=data.frame()
cutoff=0.05
for(i in 1:length(gwas.hbtrc) ){
  if(class(gwas.hbtrc[[i]])=="try-error")
    next()
  sub.res1=subset(gwas.hbtrc[[i]][["res1"]],Pc1df < 1e-3)
  #sub.res1=subset(gwas.hbtrc[[i]][["res1"]],p.adj < cutoff)
  snp=as.vector(sub.res1$id);
  if(length(snp) == 0)
     next()
  sub.res2= gwas.hbtrc[[i]][["res2"]][snp,]
  sub.res3= gwas.hbtrc[[i]][["res3"]][snp,]
  new.snp=sapply(snp, function(z){
    if(length(grep(":",z)) == 0 )
      return(z)
    zz=strsplit(z,":")[[1]]
    if(length(grep("rs\\d+",z))!=0)
      return(zz[1])
    cr=as.numeric(zz[1])
    loc=as.numeric(zz[2])
    sub.snp.ann = snp.list[[cr]]
    snp.id=as.vector(sub.snp.ann$V3)[which.min(abs(as.vector(sub.snp.ann$V2) - loc))]
    return(snp.id)
  })
  if(is.null(gwas.hbtrc[[i]][["fd.ceradsc"]]))
    gwas.hbtrc[[i]][["fd.ceradsc"]]=1
   op=as.data.frame(cbind(id=new.snp, pid=snp, effB=as.vector(sub.res1$effB), lambda=as.vector(sub.res1$lambda),se=as.vector(sub.res1$se), P1df1=as.vector(sub.res1$P1df),P1df2=as.vector(sub.res2$P1df),P1df3=as.vector(sub.res3$P1df), p.adj=as.vector(sub.res1$p.adj), pa=gwas.hbtrc[[i]][["pa"]], tp=gwas.hbtrc[[i]][["tp"]], n1=length(strsplit(gwas.hbtrc[[i]][["pas1"]],",")[[1]]),n2=length(strsplit(gwas.hbtrc[[i]][["pas2"]],",")[[1]]),n3=length(strsplit(gwas.hbtrc[[i]][["pas3"]],",")[[1]])))
  outcome2=rbind(op, outcome2)
}
outcome2$lambda=as.numeric(as.vector(outcome2$lambda))
new.outcome2=subset(outcome2,lambda < 1.1  )
#new.outcome2=subset(outcome2,tp != "dif" & tp != "sig" )
have.snp=unique(as.vector(new.outcome2$id))

library(biomaRt)
snpmart = useMart(biomart = "ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
test.snp=getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "ensembl_gene_stable_id","ensembl_type","consequence_type_tv"),filters = c("snp_filter"),values = have.snp,mart = snpmart)
id <- read.table("../pasgwas/ensembl_gene.txt",header=T,sep="\t")
row.names(id)<-as.vector(id$ensembl)
test.snp=cbind(test.snp, gene=as.vector(id[as.vector(test.snp$ensembl_gene_stable_id),]$gene))
new.test.snp=data.frame()
for(snp in unique(as.vector(test.snp$refsnp_id))){
  tt=subset(test.snp, refsnp_id==snp)
  if(nrow(tt)==1){
    new.test.snp=rbind(tt,new.test.snp)
  }else{
    gg=unique(as.vector(tt$gene))
    gg=gg[!is.na(gg)]
    tt$gene=paste(gg,collapse=",")
    new.test.snp=rbind(tt[1,],new.test.snp)
  }
}
row.names(new.test.snp) <-as.vector(new.test.snp$refsnp_id);
has.snp=as.vector(new.outcome2$id)
results2=cbind(new.outcome2, new.test.snp[has.snp, c("gene","chr_name","chrom_start","consequence_type_tv")])
results2=results2[order(as.numeric(as.vector(results2$P1df1))),]




###################################################################################################

has.snp=unique(as.vector(new.results$pid))
smps=intersect(row.names(phdata(input[["chr12"]])) ,row.names(phdata(input[["chr16"]])))
mx=matrix(nrow=0, ncol=length(smps))
colnames(mx) <- smps
mx=as.data.frame(mx)
for(i in 1:22){
  print(i)
 zz=paste("chr",i,sep="");
 input2=input[[zz]]; 
 sel=has.snp[has.snp %in% colnames(gtdata(input2))];
 if(length(sel)==0)
   next()
 sub.input2=input2[smps, sel]
 t(as.genotype(sub.input2@gtdata))->gt
 mx=rbind(mx, gt)
}
mx=as.matrix(mx)
ids=row.names(mx)
snps=as.vector(new.results[ids, ]$id)
pp1=as.vector(results[ids, ]$P1df1)
pp2=as.vector(results[ids, ]$P1df2)
pp3=as.vector(results[ids, ]$P1df3)
pas=as.vector(results[ids, ]$pa)
tp=as.vector(results[ids, ]$tp)
gene=as.vector(results[ids, ]$gene)

cols=rainbow(12)
smps=colnames(mx)
for(i in 1:nrow(mx)){
  snp=snps[i];
  id=ids[i]
  pas1=pas.list[[tp[i]]][[pas[i]]][["pas1"]]
  pas2=pas.list[[tp[i]]][[pas[i]]][["pas2"]]
  pas3=pas.list[[tp[i]]][[pas[i]]][["pas3"]]
  pas1=sub("X","",pas1)
  pas2=sub("X","",pas2)
  pas3=sub("X","",pas3)
  
  x=mx[i,smps %in% c(pas1,pas2,pas3)]
  labs=rep("unknown",length(x))
  new.smps=names(x)
  labs[new.smps %in% pas1]="Module"
  labs[new.smps %in% pas2]="non-Module"
  labs[new.smps %in% pas3]="Control"
  tb=table(x, labs)
  cs=colSums(tb);
  new.tb=sapply(1:ncol(tb), function(x) tb[,x]/cs[x])
  colnames(new.tb) <- colnames(tb)
  
  nns=paste(c("Module","non-Module","Control"),"(n=",c(length(pas1), length(pas2),length(pas3)),")",sep="")
  pdf(paste("pic_hbtrc/",snp,"_",tp[i],".pdf",sep=""), width=6,height=6)
  barplot(new.tb[,c("Module","non-Module","Control")],legend.text=row.names(tb),col=cols[c(4,8,2)], main=paste(snp,", ",gene[i],", p=", signif(as.numeric(pp1[i]),digits=2), sep=""), names.arg=nns)
  dev.off()
}

ad.gene=c("CR1","BIN1","INPP5D","HLA-DRB1","TREM2","CD2AP","NYAP1","EPHA1","PTK2B","CLU","SPI1","MS4A2","PICALM","SORL1","FERMT2","SLC24A4","ABCA7","APOE","CASS4","ECHDC3","ACE","MCF2C","NME8","MAPT","APP")


ad.ge=read.table("ad_gene.txt",header=T,sep="\t")

opp=data.frame()
for(i in 1:23){
  print(i)
  ge=as.vector(ad.ge[i,"symbol"])
  sub.results2=results2[!is.na(as.vector(results2$gene)) & as.vector(results2$gene) ==ge,]
  #sub.results2=sub.results2[ grep(ge,as.vector(sub.results2$gene)),]
  
  if(nrow(sub.results2)==0){
    cr=as.vector(ad.ge[i,"chr"])
    from=as.vector(ad.ge[i,"start"])
    to=as.vector(ad.ge[i,"end"])
    sub.results2=subset(results2, chr_name==cr & chrom_start > from -  500000 & chrom_start < to + 50000)
  }
  if(nrow(sub.results2)==0)
    next()
  opp=rbind(sub.results2[which.min(as.vector(sub.results2$P1df1)),],opp)
}

smps=row.names(phdata(input[["chr1"]]))
smps=smps[smps!="15995" & smps!="16105" & smps!="21926" & smps!="22001"]
mx=matrix(nrow=0, ncol=length(smps))
colnames(mx) <- smps
mx=as.data.frame(mx)
has.snp=unique(as.vector(opp$pid))
for(i in 1:22){
  print(i)
 zz=paste("chr",i,sep="");
 input2=input[[zz]]; 
 snpnames=input2@gtdata@snpnames
 sel=has.snp[has.snp %in% snpnames]
 if(length(sel)==0)
   next()
 sub.input2=input2[smps, sel]
 t(as.genotype(sub.input2@gtdata))->gt
 mx=rbind(mx, gt)
}

mx=as.matrix(mx)
colnames(mx)<-paste("X",colnames(mx),sep="" )


cols=rainbow(12)
smps=colnames(mx)
for(i in 1:nrow(opp)){
  snp=as.vector(as.matrix(opp[i,"pid"]))
  new.snp=as.vector(as.matrix(opp[i,"id"]))
  pa =as.vector(as.matrix(opp[i,"pa"]))
  tp =as.vector(as.matrix(opp[i,"tp"]))
  p =as.vector(as.matrix(opp[i,"P1df1"]))
  pas1=pas.list[[tp]][[pa]][["pas1"]]
  pas2=pas.list[[tp]][[pa]][["pas2"]]
  pas3=pas.list[[tp]][[pa]][["pas3"]]
  pas1=pas1[pas1 %in% smps]
  pas2=pas2[pas2 %in% smps]
  pas3=pas3[pas3 %in% smps]
  ge=as.vector(as.matrix(opp[i,"gene"]))
    
  labs=c(rep("Module",length(pas1)), rep("non-Module",length(pas2)), rep("Control",length(pas3)))
  gt=mx[snp, c(pas1,pas2,pas3)]
  gt[is.na(gt)]<-"1/1"
  tb=table(gt,labs)
  tb1=as.matrix(table(gt,labs))
  tb0=table(labs)
  tb2=tb1[,c(1,2)]
  tb3=tb2
  tb3[,2]=tb1[,2]+tb1[,3]
  colnames(tb3)<-c("Control","AD")
  p2=chisq.test(tb3*100)$p.value
  tb[,1]=tb[,1]/tb0[1]
  tb[,2]=tb[,2]/tb0[2]
  tb[,3]=tb[,3]/tb0[3]
  tb=as.matrix(tb)
  lbs=paste(c("Module","non-Module","Control"),"\n(n=",tb0[c("Module","non-Module","Control")],")",sep="")
  pdf(paste("pic_hbtrc2/",new.snp,"_",pa,"_",tp,".pdf",sep=""),width=5,height=5)
  barplot(tb[,c("Module","non-Module","Control")], col=cols[c(4,8,2)],names.arg=lbs,main=paste(new.snp,", ",ge,", p=",signif(as.numeric(p),digits=3), sep=""))
  dev.off()
}


opp2=data.frame()
for(i in 1:23){
  print(i)
  ge=as.vector(ad.ge[i,"symbol"])
  sub.results0=results0[!is.na(as.vector(results0$gene)) & as.vector(results0$gene) ==ge,]
  #sub.results2=sub.results2[ grep(ge,as.vector(sub.results2$gene)),]
  
  if(nrow(sub.results0)==0){
    cr=as.vector(ad.ge[i,"chr"])
    from=as.vector(ad.ge[i,"start"])
    to=as.vector(ad.ge[i,"end"])
    sub.results0=subset(results0, chr_name==cr & chrom_start > from -  500000 & chrom_start < to + 50000)
  }
  if(nrow(sub.results0)==0)
    next()
  print(dim(sub.results0))
  opp2=rbind(sub.results0[which.min(as.vector(sub.results0$P1df)),],opp2)
}




####################################
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

gg.manhattan <- function(df, threshold, hlight, col, ylims, title){
  # format df
  df.tmp <- df %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
    mutate( is_annotate=ifelse(SNP %in% hlight, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  ggplot(df.tmp, aes(x=BPcum, y=-log10(P))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), size=1) +
    scale_color_manual(values = rep(col, 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
    
    # add plot and axis titles
    labs(x = "") +
    
    # add genome-wide sig and sugg lines
    geom_hline(yintercept = -log10(sig)) +
    geom_hline(yintercept = -log10(sugg), linetype="dashed") +
    # Add highlighted points
    geom_point(data=subset(df.tmp, is_highlight=="yes"), color="green", size=1) +
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(gene)), size=4, force=5,box.padding = unit(0.1, "lines"), point.padding = unit(0.1, "lines"),arrow = arrow(length = unit(0.01, 'npc')) ) +
    # Custom the theme:
    theme_bw(base_size = 15) +
    theme( 
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}


dd1=results2[,c("id","chr_name","chrom_start","P1df1","gene")]
names(dd1)<-c("SNP", "CHR", "BP", "P","gene")
dd1$BP=as.numeric(as.vector(dd1$BP))
dd1$CHR=as.numeric(as.vector(dd1$CHR))
dd1$P=as.numeric(as.vector(dd1$P))
dd1=dd1[order(dd1$P),]
new.dd1=dd1[!duplicated(dd1$SNP),]

library(wesanderson)
mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B") # chr color palette
mypalette <- wes_palette(n=5, name="FantasticFox1")
mysnps <- as.vector(subset(new.dd1, P< 5e-8)$SNP)[1:10]
mysnps=c("rs4788579","rs146624252","rs3867593","rs113337484","rs9912864","rs769450","rs80167208","rs34233526","rs11253483","rs72129870")
mygene=c("IST1","PDE1A","ZBTB4","AKIRIN2","NTN1","APOE","BOC","MTHFD1L","LARP4B","VAV3")

sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line

p1=gg.manhattan(new.dd1, threshold=5e-8, hlight=mysnps, col=mypalette, ylims=c(4,9), title="")




dd2=results0[,c("id","chr_name","chrom_start","P1df","gene")]
names(dd2)<-c("SNP", "CHR", "BP", "P","gene")
dd2$BP=as.numeric(as.vector(dd2$BP))
dd2$CHR=as.numeric(as.vector(dd2$CHR))
dd2$P=as.numeric(as.vector(dd2$P))
dd2=dd2[order(dd2$P),]
new.dd2=dd2[!duplicated(dd2$SNP),]


mysnps2 <- c("rs2405283","rs769450")
sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line
p2=gg.manhattan(new.dd2, threshold=1e-6, hlight=mysnps2, col=mypalette, ylims=c(3,9), title="")

library(gridExtra )
pdf("mh.pdf",width=10,height=7)
grid.arrange(p2, p1, nrow = 2)
dev.off()

########## go for spDEG
library(clusterProfiler)
library(org.Hs.eg.db)
seed.gene=ad[ad %in% names(res.deg)]
go.deg=lapply(ad[ad %in% names(res.deg)], function(x){
    sig.ges=res.deg[[x]][["genes"]]
    go <- enrichGO(gene         = sig.ges,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.1)
  return(go@result)
})

names(go.deg) <- seed.gene
ip=read.table("input.txt",header=T,sep="\t")


tem=vector()
for(x in names(go.deg)){
 tem=append(tem, as.vector(subset(go.deg[[x]], p.adjust < 0.01)$Description) )
}
tem=unique(tem)

output=data.frame()

for(x in names(go.deg)){
   subx =subset(go.deg[[x]], p.adjust < 0.05 )
   subx=subx[subx$Description %in% tem, c("Description","p.adjust")] 
   if(nrow(subx)==0)
     next()
   gwg=as.vector(subset(ip, seed.pat.==x)$gene)
   if(length(gwg)==0)
     next()
   for(gg in gwg){
     output=rbind(data.frame(gene=gg, subx), output)
   }
}
output=read.table("output.txt",sep="\t",header=T)
library(igraph)
graph_from_data_frame(output, directed = TRUE)->gg


##############

load("tag.rda")
tag=paste(as.vector(results$pa), as.vector(results$tp),sep="/")
results.braak=results[tag %in% tag.braak,]
results.Atrophy=results[tag %in% tag.Atrophy,]
new.results.braak=merge.results(results.braak)
new.results.Atrophy=merge.results(results.Atrophy)



###
new.results=merge.results2(results)
tag=paste(as.vector(results$pa), as.vector(results$tp),sep="#")
cc.results=sapply(unique(tag), function(x){
  sub.results=new.results[tag==x,]
  return(length(unique(as.vector(sub.results$group))))
})


#### evaluation based on simulation