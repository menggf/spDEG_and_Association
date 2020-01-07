
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

