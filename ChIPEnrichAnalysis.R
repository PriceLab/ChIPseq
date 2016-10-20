# find enrichments of genes associated with Smad3 binding 
# for GO terms, TReNA modules, HD genes

library(chipenrich)
library(trena)
library( edgeR )
library( doBy )

# 1. Use ChIP-enrich to count peaks within 10kb of each TSS

# note: chipenrich() calculates GO enrichment scores.
# However, these scores seem to not work well with our
# data, since they do not distinguish genes with high regulatory
# potnential (several Smad3 binding sites) 
# from those with low regulatory potential (just 1-2 binding sites)
 
setwd("/proj/price1/sament/chipseq/Smad3")
peakfile="Smad3_peaks.reproducible_fdr01+10reads.bed"
peaks = read.table( peakfile )
peakres = read.csv("Smad3.MACS.edgeR.results.2016-6-29.csv")
counts = peakres[,2:5]
libsize = read.table( "libsize.txt" )
cpm = cpm( counts , lib.size = libsize[1:4,1] )
chrom = gsub( "\\:(.*)" , "" , peakres[,1] )
pos = gsub( "(.*)\\:" , "" , peakres[,1] )
start = gsub( "\\-(.*)" , "" , pos )
end = gsub( "(.*)\\-" , "" , pos )
peakloc = data.frame( chrom , start , end )


enrich = 
chipenrich( 
   peaks = peakfile , 
   genome = "mm9" ,
   locusdef = "10kb" ,
   genesets = "reactome" )

# we will focus on the genes with the most Smad3 binding sites
# specifically, those with peak count >2 S.D. above the mean count within 10kb

peakgenes = enrich$peaks
x = cbind( peakloc , cpm )
peakgenes = merge( peakgenes , x , by.x = 2:4 , by.y = 1:3 )
mean_cpm = rowMeans( peakgenes[,14:17] )
genesum = summaryBy( mean_cpm~peakgenes.gene_symbol , 
   data = data.frame( peakgenes$gene_symbol , mean_cpm ) ,
   fun.names = sum )
cpm.gene = genesum[,2]
names(cpm.gene) = genesum[,1]
colnames(genesum) = c("gene_symbol" , "cpm")
genesum = merge( genes.mmu , genesum , by.x = 2 , by.y = 1 , all.x = T )
genesum$cpm[ is.na( genesum$cpm) ] = 0
genesum$logcpm = log(genesum$cpm+1)
genesum$cpm.zscore = ( genesum$logcpm - mean(genesum$logcpm) ) / sd(genesum$logcpm)

peaks_per_gene = enrich$peaks_per_gene
peaks_per_gene = merge( genes.mmu , peaks_per_gene , by = 1 , all.x = T )
peaks_per_gene$num_peaks[ is.na(peaks_per_gene$num_peaks) ] = 0
peaks_per_gene$log_num_peaks = log( peaks_per_gene$num_peaks+1 )
peaks_per_gene$num_peaks.zscore = 
   (peaks_per_gene$log_num_peaks-mean(peaks_per_gene$log_num_peaks)) / sd(peaks_per_gene$log_num_peaks)

peaks_per_gene = merge( peaks_per_gene , genesum , by.x = 2 , by.y = 1 )
 
peaks_per_gene$cpmXcount = peaks_per_gene$num_peaks.zscore * peaks_per_gene$cpm.zscore

save( enrich , genes.mmu , peaks_per_gene , file="chipenrich_results_smad3.RData" )

# 2. Find enrichments with the Smad3 module predicted by TReNA

setwd("/proj/price1/CHDI/users/sament/allelic_series_striatum")
print(load("trn_2015_nov.RData"))
trn = trimTRN( trn0 )

smad3.trena = trn[ , "Smad3" , drop = F ]

genes = genes.mmu$SYMBOL

a = peaks_per_gene[ peaks_per_gene$num_peaks.zscore > 2 , 1 ]
b = rownames(smad3.trena)[ smad3.trena != 0 ]
t = table( genes %in% a , genes %in% b )
fisher.test( t )
#        FALSE  TRUE
#  FALSE 21524   875
#  TRUE    625    55
#p-value = 1.32e-06
#odds ratio 
#  2.165 

c = peaks_per_gene[ peaks_per_gene$num_peaks > 0 , 1 ]
t = table( genes %in% c , genes %in% b )
fisher.test( t )
# p = 2.8e-84
# odds ratio = 4.33

# 3. GO enrichment analysis with Fisher's exact test

load("/proj/price1/sament/resources/goterms_mmu_2015-11-4.RData")
n = length(genesets)
size = sapply( 1:n , function(x) length( genesets[[x]] ) )
sets = genesets[ size <= 500 & size >= 10 ]
a = peaks_per_gene[ peaks_per_gene$num_peaks.zscore > 2 , "SYMBOL" ]
n = length(sets)
p = rep( NA , n )
est = rep( NA , n )
nGenes = rep(NA , n )
setGenes = rep( NA , n )
all.genes = genes.mmu$SYMBOL
for( i in 1:n ) {
  set = sets[[i]]
  tab = table( all.genes %in% a , all.genes %in% set )
  if( any( nrow(tab) != 2 , ncol(tab) != 2 ) ) next
  test = fisher.test( tab , alternative = "greater" )
  p[i] = test$p.value
  est[i] = test$estimate
  nGenes[i] = tab[2,2]
  setGenes[i] = paste( intersect( a , set ) , collapse = ";" , sep = ";" )
  if( i/100 == round(i/100) ) cat(i,"\n")
}
size = sapply( 1:n , function(x) length( sets[[x]] ) )
fisher.res = data.frame( set = names(sets) , size , nGenes , est , p )
fisher.res$q =  p.adjust( p , method = "BH" )
fisher.res$setGenes = setGenes

fisher.res[ fisher.res$q < 0.01 , 1:3 ]

write.csv( fisher.res , file="GO_enrichments.top_Smad3_targets.2016-8-1.csv" , row.names=F 


# 4. enrichment for HD DEGs

setwd("/proj/price1/CHDI/users/sament/allelic_series_striatum")
print(load("edgeR.lrt.RData"))
edgeR = res

edgeR.genes = edgeR$gene
fdr = grep("FDR",colnames(edgeR))
fdr = edgeR[,fdr]
pval = grep("PValue",colnames(edgeR))
pval = edgeR[,pval]
logfc = grep("logFC",colnames(edgeR))
logfc = edgeR[,logfc]

# define DEGs bidirectionally ( both up- and down-regulated )
 
degs.fdr05 = list()
degs.pval01 = list()
for( i in 1:ncol(fdr) ) {
  fdr05 = fdr[,i] < 0.05
  fdr05.genes = edgeR.genes[ fdr05 ]
  degs.fdr05[[i]] = fdr05.genes
  p01 = pval[,i] < 0.01
  p01.genes = edgeR.genes[ p01 ]
  degs.pval01[[i]] = p01.genes
}
names(degs.fdr05) = names(degs.pval01) = gsub("FDR.","",colnames(fdr))

all.genes = intersect( edgeR.genes , genes.mmu$SYMBOL )
a = intersect( all.genes , peaks_per_gene[ peaks_per_gene$num_peaks.zscore > 2 , "SYMBOL" ] )
n = ncol(fdr)
p.fdr05 = rep( 1 , n )
est.fdr05 = rep( 0 , n )
nGenes.fdr05 = rep( 0 , n )
p.p01 = rep(1,n)
est.p01 = rep(0,n)
nGenes.p01 = rep(0,n)
for( i in 1:n ) {
  set.fdr05 = degs.fdr05[[i]]
  set.p01 = degs.pval01[[i]]
  t1 = table( all.genes %in% a , all.genes %in% set.fdr05 )
  if( all( nrow(t1) == 2 , ncol(t1) == 2 ) ) {
    test = fisher.test( t1 , alternative="greater" )
    p.fdr05[i] = test$p.value
    est.fdr05[i] = test$estimate
    nGenes.fdr05[i] = t1[2,2]
  }
  t2 = table( all.genes %in% a , all.genes %in% set.p01 )
  if( all( nrow(t2) == 2 , ncol(t2) == 2 )) {
    test = fisher.test( t2 , alternative="greater" )
    p.p01[i] = test$p.value
    est.p01[i] = test$estimate
    nGenes.p01[i] = t2[2,2]
  }
}

edgeR.fisher.res.bidirectional = data.frame(
  condition = names(degs.fdr05) ,
  nDEGs.fdr05 = sapply( 1:n , function(x) length(degs.fdr05[[x]]) ) ,
  nModDEGs.fdr05 = nGenes.fdr05 ,
  est.fdr05 ,
  p.fdr05 ,
  q.fdr05 = p.adjust( p.fdr05 ) ,
  nDEGs.p01 = sapply( 1:n , function(x) length(degs.pval01[[x]]) ) ,
  nModDEGs.p01 = nGenes.p01 ,
  est.p01 ,
  p.p01 ,
  q.p01 = p.adjust( p.p01 )
)

setwd("/proj/price1/sament/chipseq/Smad3")
write.csv( edgeR.fisher.res.bidirectional , row.names=F ,
   file = "fishertest.Smad3.bidirectional_degs.csv")

# down-regulated genes

degs.fdr05 = list()
degs.pval01 = list()
for( i in 1:ncol(fdr) ) {
  fdr05 = fdr[,i] < 0.05 & logfc[,i] < 0
  fdr05.genes = edgeR.genes[ fdr05 ]
  degs.fdr05[[i]] = fdr05.genes
  p01 = pval[,i] < 0.01 & logfc[,i] < 0
  p01.genes = edgeR.genes[ p01 ]
  degs.pval01[[i]] = p01.genes
}
names(degs.fdr05) = names(degs.pval01) = gsub("FDR.","",colnames(fdr))

all.genes = intersect( edgeR.genes , genes.mmu$SYMBOL )
a = intersect( all.genes , peaks_per_gene[ peaks_per_gene$num_peaks.zscore > 2 , "SYMBOL" ] )
n = ncol(fdr)
p.fdr05 = rep( 1 , n )
est.fdr05 = rep( 0 , n )
nGenes.fdr05 = rep( 0 , n )
p.p01 = rep(1,n)
est.p01 = rep(0,n)
nGenes.p01 = rep(0,n)
for( i in 1:n ) {
  set.fdr05 = degs.fdr05[[i]]
  set.p01 = degs.pval01[[i]]
  t1 = table( all.genes %in% a , all.genes %in% set.fdr05 )
  if( all( nrow(t1) == 2 , ncol(t1) == 2 ) ) {
    test = fisher.test( t1 , alternative="greater" )
    p.fdr05[i] = test$p.value
    est.fdr05[i] = test$estimate
    nGenes.fdr05[i] = t1[2,2]
  }
  t2 = table( all.genes %in% a , all.genes %in% set.p01 )
  if( all( nrow(t2) == 2 , ncol(t2) == 2 )) {
    test = fisher.test( t2 , alternative="greater" )
    p.p01[i] = test$p.value
    est.p01[i] = test$estimate
    nGenes.p01[i] = t2[2,2]
  }
}

edgeR.fisher.res.down = data.frame(
  condition = names(degs.fdr05) ,
  nDEGs.fdr05 = sapply( 1:n , function(x) length(degs.fdr05[[x]]) ) ,
  nModDEGs.fdr05 = nGenes.fdr05 ,
  est.fdr05 ,
  p.fdr05 ,
  q.fdr05 = p.adjust( p.fdr05 ) ,
  nDEGs.p01 = sapply( 1:n , function(x) length(degs.pval01[[x]]) ) ,
  nModDEGs.p01 = nGenes.p01 ,
  est.p01 ,
  p.p01 ,
  q.p01 = p.adjust( p.p01 )
)

write.csv( edgeR.fisher.res.down , row.names=F ,
   file = "fishertest.Smad3.down_degs.csv")

all.fdr05.down = unique( unlist( degs.fdr05 ))
t = table( all.genes %in% all.fdr05.down , all.genes %in% a )
fisher.test( t )
t 

# up-regulated genes

degs.fdr05 = list()
degs.pval01 = list()
for( i in 1:ncol(fdr) ) {
  fdr05 = fdr[,i] < 0.05 & logfc[,i] > 0
  fdr05.genes = edgeR.genes[ fdr05 ]
  degs.fdr05[[i]] = fdr05.genes
  p01 = pval[,i] < 0.01 & logfc[,i] > 0
  p01.genes = edgeR.genes[ p01 ]
  degs.pval01[[i]] = p01.genes
}
names(degs.fdr05) = names(degs.pval01) = gsub("FDR.","",colnames(fdr))

all.genes = intersect( edgeR.genes , genes.mmu$SYMBOL )
a = intersect( all.genes , peaks_per_gene[ peaks_per_gene$num_peaks.zscore > 2 , "SYMBOL" ] )
n = ncol(fdr)
p.fdr05 = rep( 1 , n )
est.fdr05 = rep( 0 , n )
nGenes.fdr05 = rep( 0 , n )
p.p01 = rep(1,n)
est.p01 = rep(0,n)
nGenes.p01 = rep(0,n)
for( i in 1:n ) {
  set.fdr05 = degs.fdr05[[i]]
  set.p01 = degs.pval01[[i]]
  t1 = table( all.genes %in% a , all.genes %in% set.fdr05 )
  if( all( nrow(t1) == 2 , ncol(t1) == 2 ) ) {
    test = fisher.test( t1 , alternative="greater" )
    p.fdr05[i] = test$p.value
    est.fdr05[i] = test$estimate
    nGenes.fdr05[i] = t1[2,2]
  }
  t2 = table( all.genes %in% a , all.genes %in% set.p01 )
  if( all( nrow(t2) == 2 , ncol(t2) == 2 )) {
    test = fisher.test( t2 , alternative="greater" )
    p.p01[i] = test$p.value
    est.p01[i] = test$estimate
    nGenes.p01[i] = t2[2,2]
  }
}

edgeR.fisher.res.up = data.frame(
  condition = names(degs.fdr05) ,
  nDEGs.fdr05 = sapply( 1:n , function(x) length(degs.fdr05[[x]]) ) ,
  nModDEGs.fdr05 = nGenes.fdr05 ,
  est.fdr05 ,
  p.fdr05 ,
  q.fdr05 = p.adjust( p.fdr05 ) ,
  nDEGs.p01 = sapply( 1:n , function(x) length(degs.pval01[[x]]) ) ,
  nModDEGs.p01 = nGenes.p01 ,
  est.p01 ,
  p.p01 ,
  q.p01 = p.adjust( p.p01 )
)

write.csv( edgeR.fisher.res.up , row.names=F ,
   file = "fishertest.Smad3.up_degs.csv")

all.fdr05.up = unique( unlist( degs.fdr05 ))
t = table( all.genes %in% all.fdr05.up , all.genes %in% a )
fisher.test( t )
t


