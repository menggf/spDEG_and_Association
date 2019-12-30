## merge the SNPs within a genomic region
merge.results<-function(dfs){
  chrs=unique(as.vector(dfs$chr_name))
  chrs=chrs[!is.na(chrs)]
  tmp=data.frame()
  for(chr in chrs){
    sub.dfs=subset(dfs,chr_name==chr)
    sub.dfs=sub.dfs[order(as.numeric(as.vector(sub.dfs$chrom_start))),]
    cc= as.numeric(as.vector(sub.dfs$n1))
    pos=as.numeric(as.vector(sub.dfs$chrom_start))
    pv= as.numeric(as.vector(sub.dfs$P1df1))
    
    mm=abs(sapply(pos,function(x) pos-x))
    lbs=rep(1, length(pos))
    if(length(pos)> 1){
      p=1
      q=1
      for(i in 2:length(pos)){
        if(mm[p, i] > 5e4){
          p=i
          q=q+1
        }
        lbs[i]=q
      }
      wh=vector()
      for(i in unique(lbs)){
        sub.pv=pv
        sub.pv[lbs!=i]=1
        wh=append(wh, which.min(sub.pv))
      }
    }
    grp=paste(chr,lbs,sep="#")
    cc.lbs=table(grp)
    new.dfs=cbind(sub.dfs, group=grp,cc=as.vector(cc.lbs[grp]))[wh,]
    tmp=rbind(new.dfs,tmp)
  }
  return(tmp)
}

## merge the SNPs within a genomic region
merge.results2<-function(dfs){
  chrs=unique(as.vector(dfs$chr_name))
  chrs=chrs[!is.na(chrs)]
  tmp=data.frame()
  for(chr in chrs){
    sub.dfs=subset(dfs,chr_name==chr)
    if(nrow(sub.dfs)==0)
      next()
    sub.dfs=sub.dfs[order(as.numeric(as.vector(sub.dfs$chrom_start))),]
    pos=as.numeric(as.vector(sub.dfs$chrom_start))
    pv=as.numeric(as.vector(sub.dfs$P1df1))
    mm=abs(sapply(pos,function(x) pos-x))
    lbs=rep(1, length(pos))
    if(length(pos)> 1){
      p=1
      q=1
      for(i in 2:length(pos)){
        if(mm[p, i] > 1e4){
          p=i
          q=q+1
        }
        lbs[i]=q
      }
      wh=vector()
      for(i in unique(lbs)){
        sub.pv=pv
        sub.pv[lbs!=i]=1
        wh=append(wh, which.min(sub.pv))
      }
    }
    grp=paste(chr,lbs,sep="#")
    cc.lbs=table(grp)
    new.dfs=cbind(sub.dfs, group=grp,cc=as.vector(cc.lbs[grp]))
    tmp=rbind(new.dfs,tmp)
  }
  return(tmp)
}
