##  Create QQ plots, Volcano Plots, HeatMap of sign of logFC,
##  tSNPE or PCA Cluster based on sign of log FC
##  SHARED SIGNIFICANT GSs
library(qvalue)
library(ggpubr)
library(annotables)

source("/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Code/Summary_funcs.r")

GSEADir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Results/GSEA/"
GSEADirGTEx <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/GSEA/"
ResultDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Results/"
SummaryDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Summary/"
PlotDir <- paste0(ResultDir,"Plots/")

GSEAOutDir <- paste0(ResultDir,"GSEA/")
GSEAGSFiles <- c("/home/askol/bin/GSEA_genesets/Custom/Hormone_Immune_Custom.gmt",
                "/home/askol/bin/GSEA_genesets/msigdb_v6.2_GMTs/msigdb.v6.2.symbols.gmt")

setwd(SummaryDir)


## projects <- c("TCGA-SKCM", "TCGA-THCA", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LAML")
projects = read.table(file = "/gpfs/data/stranger-lab/askol/TCGA/TCGA_Target_projects.txt",
    header=TRUE, as.is=TRUE, sep="\t")
projects = projects[grep("TCGA",projects[,1]),1]

skip.cancers <- c("TCGA-BRCA","TCGA-OV", "TCGA-CESC", "TCGA-PRAD", "TCGA-TGCT",
                  "TCGA-UCS", "TCGA-UCEC")
projects <- projects[!projects %in% skip.cancers] 

proj.tis.match <- rbind(c("TCGA-SKCM", "Skin_Not_Sun_Exposed_Suprapubic"),
                        c("TCGA-SKCM", "Skin_Sun_Exposed_Lower_leg"),
                        c("TCGA-THCA", "Thyroid"),
                        c("TCGA-LIHC", "Liver"),
                        c("TCGA-LUAD", "Lung"),
                        c("TCGA-LAML", "Cells_EBV-transformed_lymphocytes"),
                        c("TCGA-LIHC", "Skin_Not_Sun_Exposed_Suprapubic"),
                        c("TCGA-LUAD", "Thyroid"))
                        
tmp <- collect.results(projects, ResultDir)

logFC.tcga <- tmp$logFC
logPs.tcga <- tmp$logPs
Qs.tcga <- tmp$Qs
geneInfo.tcga <- tmp$GeneInfo

geneInfo <- grch38
geneInfo <- select(geneInfo, -end, -description)
## REMOVE GENES THAT MAP TO MULTIPLE CHROMOSOMES ##
geneInfo <- removeMultiMapGenes(geneInfo)

## ####################### ##
##
## BRING IN GTEX TISSUES
##
## ####################### ##

GSEADirGtex <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/GSEA/"

tissues <- read.table(file = paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                          "data/support_files/all_v8_tissues_both_sexes.txt"), as.is=T)
tissues <- unlist(tissues)
ResultDirGtex <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/"

tmp <- collect.results.gtex(tissues, ResultDirGtex)

logFC.gtex <- tmp$logFC
logPs.gtex <- tmp$logPs
Qs.gtex <- tmp$Qs
geneInfo.gtex <- tmp$GeneInfo

saveRDS(list(logFC.tcga = logFC.tcga, logPs.tcga = logPs.tcga,
             Qs.tcga = Qs.tcga, logFC.gtex = logFC.gtex, logPs.gtex = logPs.gtex,
             Qs.gtex = Qs.gtex), file = "Summary.RDS")
if (0){
    tmp <- loadRDS("Summary.RDS")
    logFC.tcga = tmp$logFC.tcga;  logPs.tcga = tmp$logPs.tcga;
    Qs.tcga = tmp$Qs.tcga logFC.gtex = tmp$logFC.gtex; logPs.gtex = tmp$logPs.gtex;
    Qs.gtex = tmp$Qs.gtex
}

#####  PLOT 1: NUMBER OF DE GENES ######
## GET DISTRIBUTION OF THE NUMBER OF DE GENES ##
get.noDE.genes(Qs.tcga, Qs.gtex,  project, tissues, PlotDir)

## WHICH GENES ARE COMMONLY DE BETWEEN GTEX AND TCGA ##


## ## WHICH GENES ARE SHARED ACROSS CANCER
## META ANALYSIS AND FDR > .10
## PLOT NUMBER OF CANCERS/TISSUES SHARED PER GENE
## AND
## PLOT NUMBER GENES FROM EACH CHROMOSOME WITH MULTIPLE CANCERS/TISSUES/GENE

file.pre = "tcga"

tmp <- get.noDE.cancers.for.gene(Qs.tcga, qthresh= .1, geneInfo,
                                       Ntrunc = 5, min.studies = 5,
                                       ResultDir, PlotDir, file.pre)
noDE.tcga <- tmp$NoDE
t <- tmp$t; t.nox <- tmp$t.nox
t
## mean number of genes different between obs and expected ##
x.ind = c(21,22)
## mean absolute differences in observered relative to expected ##
mean(abs(t$obs[,1]-t$exp[,1])/abs(t$exp[,1])) * 100
## mean absolute difference in observed relative to expected in X & Y
mean(abs(t$obs[x.ind,1] - t$exp[x.ind,1])/abs(t$exp[x.ind,1])) * 100
## mean absolute difference in observed relative to expected in autosome
mean(abs(t$obs[-x.ind,1] - t$exp[-x.ind,1])/abs(t$exp[-x.ind,1])) * 100

t.nox

file.pre = "tcga_all"
tmp <- get.noDE.cancers.for.gene(Qs.tcga, qthresh= .1, geneInfo,
                                       Ntrunc = 5, min.studies = 1,
                                 ResultDir, PlotDir, file.pre)
noDE.tcga <- tmp$NoDE
t1 <- tmp$t
## mean number of genes different between obs and expected ##
mean(abs(t1$obs[,1]-t1$exp[,1])/abs(t1$exp))
t1.nox <- tmp$t.nox
t1
t1.nox
x.ind = c(23,24)
## mean absolute differences in observered relative to expected ##
mean(abs(t1$obs[,1]-t1$exp[,1])/abs(t1$exp[,1])) * 100
## mean absolute difference in observed relative to expected in X & Y
mean(abs(t1$obs[x.ind,1] - t1$exp[x.ind,1])/abs(t1$exp[x.ind,1])) * 100
## mean absolute difference in observed relative to expected in autosome
mean(abs(t1$obs[-x.ind,1] - t1$exp[-x.ind,1])/abs(t1$exp[-x.ind,1])) * 100



## total number of DE genes ##
t <- table(noDE.tcga$n)
tot = sum(t[names(t)!="0"])
## NUMBER OF SEX DE GENES IN MORE THAN 1 CANCER ##
gt1 = sum(t[names(t)%in%c("0","1") == FALSE])
gt1/tot

gt12 <- sum(t[as.numeric(names(t)) > 12])
gt12/tot

file.pre = "gtex"
tmp <- get.noDE.cancers.for.gene(Qs.gtex, qthresh= .1, geneInfo,
                                       Ntrunc = 5, min.studies = 10,
                                       ResultDir, PlotDir, file.pre)
noDE.gtex <- tmp$NoDE
t <- tmp$t.gtex; t.nox <- tmp$t.gtex.nox


## PLOT THE DISTRIBUTION OF DE GENES ON THE X CHROMOSOME ##
x.tcga <- plot.de.on.x(logPs.tcga, logFC.tcga, geneInfo.tcga, PlotDir, file.pre="tcga")
x.gtex <- plot.de.on.x(logPs.gtex, logFC.gtex, geneInfo.gtex, PlotDir, file.pre="gtex")
    

## FIND OUT WHICH GENES ARE COMMONLY DE IN BOTH CANCER AND NORMAL TISSUE ##
tmp <- get.common.common.DE.genes(noDE.gtex, noDE.tcga,
                                  Qs.tcga,  Qs.gtex,
                                  thresh.gtex = 10, thresh.tcga=5,
                                  ResultDir, PlotDir)
Ns <- tmp$N; Qs <- tmp$Qs

## REMOVE GENES NOT OBSERVED IN AT LEAST 50% OF EACH POPULATION ##
N.max.gtex <- max(Ns$N.gtex)
N.max.tcga <- max(Ns$N.tcga)
Ns <- filter(Ns, N.gtex >= floor(N.max.gtex/2) &
                  N.tcga >= floor(N.max.tcga/2))
Ns.dist <- get.n.dist(Ns)


## INTERESTING GENES ##
ind <- which(Ns.dist$interest.gt == 1| Ns.dist$interest.gt2 == 1 |
             Ns.dist$interest.lt == 1| Ns.dist$interest.lt2 == 1)
out <- Ns.dist[ind,]
out <- select(out, symbol, chr, start, n.tcga, n.tcga.imp, N.tcga, n.gtex, n.gtex.imp, N.gtex, prob.g.g.t.ge, prob.t.g.g.ge)

## WRITE OUT GENES THAT ARE OUTLIERS WITH RESPECT TO BE FOUND MORE OR LESS IN TCGA
## THAN IN GTEX
file <- paste0(ResultDir, "DiscordantNoDE_TCGA.v.GTEX.txt")
write.table(file = file, out, quote=FALSE, row.names=F, col.names=T)

file <- paste0(PlotDir,"Dist_DE_Conditional.pdf")

plot.cond.de.distribution(Ns.dist, file=file)

## PLOT PCA BASED ON DE IN TCGA ##
file <- paste0(PlotDir,"TCGA_DE_PCA.pdf")
plot.PCA(logPs.tcga, file, geneInfo.tcga)
plot.PCA(logPs.tcga, file, geneInfo.tcga, exclX=T)

## PLOT PCA BASED ON DE IN GTEX ##
file <- paste0(PlotDir,"GTEx_DE_PCA.pdf")
plot.PCA(logPs.gtex, file, geneInfo.gtex)
plot.PCA(logPs.gtex, file, geneInfo.gtex, exclX=T)


## ############################################################ ##
##                 GSEA ANALYSIS SUMMARY                        ##
## ############################################################ ##

## BRING IN RESULTS OF GSEA ANALYSIS ##
gsea <- get.fgsea.results(projects, GSEADir)

make.cytoscape.out(gsea)

## HOW MANY GENE SETS HAVE Q-VALUES LESS THAN 0.05 using ALL CHROMOSOMES ##
t <- gsea %>% filter(padj <= 0.05, nox==0) %>% select(pathway, proj) %>% table()
sum(t)
mean(t)
colSums(t)
median(colSums(t))
dim(t)
sort(table(rowSums(t)), decreasing=T)


## HOW MANY GENE SETS HAVE Q-VALUES LESS THAN 0.05 using AUTOSOMES ONLY ##
t <- gsea %>% filter(padj <= 0.05, nox==1) %>% select(pathway, proj) %>% table()
sum(t)
colSums(t)
median(colSums(t))
dim(t)
sort(table(rowSums(t)), decreasing=T)

## 
## DETERMINE OVERLAP OF CONTRIBUTING GENES
plot.file <- paste0(PlotDir, "GSEA_gene_geneset_OL.pdf")
tmp <- gsea.OL(gsea, plot.file = plot.file, topN=100, LEcutoff=3, NoX=0)
gsea.OL.X <- tmp$gseaOL
genesContrib <- tmp$genesContrib

## MEAN NUMBER OF CONTRIBUTING GENES IN TOP 100 GENES SETS WITH AND WITHOUT DUPES 
sol <- summarize.OL(genesContrib)

## DETERMINE OVERLAP OF CONTRIBUTING GENES WHEN NOT INCLUDING X ##
plot.file <- paste0(PlotDir, "GSEA_gene_geneset_OL_NoX.pdf")
tmp <- gsea.OL(gsea, plot.file = plot.file, topN=100, LEcutoff=3, NoX=1)
gsea.OL.noX.out <- tmp$gsea.OL
genesContrib.noX <- tmp$genesContrib

## MEAN NUMBER OF CONTRIBUTING GENES IN TOP 100 GENES SETS WITH AND WITHOUT DUPES 
sol.nox <- summarize.OL(genesContrib.noX)


##### REPORT THE TOP 10 GENE SETS IN EACH CANCER WITH AND WITHOUT X GENES #####
file <- paste0(ResultDir,"/Tables/Top10GSEA_X.txt")
gseaTbl <- reportTop10(gsea, NoX = 0, TopRank=10, maxlogP = 39.99, file = file)

file <- paste0(ResultDir,"/Tables/Top10GSEA_NoX.txt")
gseaTblX <-reportTop10(gsea, NoX = 1, TopRank=10, maxlogP = 39.99, file = file)

## CREATE OUTPUT FOR CYTOSCAPE ##
outDir <- paste0(ResultDir, "GSEA/CytoscapeInFiles/")
make.cytoscape.out(gsea, outDir)

#### REPORT THE OVERLAP OF GENESETS ACROSS CANCERS ##
## ALSO CREATES FILE TO USE IN CYTOSCAPE TO EXPLORE CONNECTION BETWEEN
## MOST SHARED GENESETS
gs.ol <- genesetOL(gsea, NoLE=5, qThresh=0.05, maxGS=500, ResultDir) 

## BRING GENE SET CLUSTER MEMBERSHIP OF MAJOR CLUSTERS BACK IN ##
## THESE WERE CREATED IN CYTOSCAPE
file <- "../Results/GSEA/Shared_NoX_gt6shared default node.csv"
gsclust <- read.table(file = file, as.is=T, header=T, sep=",")
cl.genes <- clust.genes(gsclust)
cl.genesets <- cl.genes$cluster.genesets
cl.genes <- cl.genes$cluster.genes

out.file <- paste0(ResultDir,"GSEA_Cluster_GeneSets_NoX_gt6shared.txt")
write.table(file = out.file, cl.genesets, quote=F, row.names=F,
            col.names=T)

out.file <- paste0(ResultDir,"GSEA_Cluster_Genes_NoX_gt6shared.txt")
## keep only those genes that are present in 1/4 or more gene sets in cluster ##
cl.genes.25 <- filter(cl.genes, count > ceiling(.25*NGeneSets))
write.table(file = out.file, cl.genes.25, quote=F, row.names=F, col.names=T)


## DETERMINE THE DISTRIBUTION OF GENES IN LEGENES OF ALL GENE SETS, SIGNIFICANT
## OR NOT
gene.gs.dist <- gs.gene.concat(gsea)

## DETERMINE THE FREQUENCY OF GENES IN THE ENRICHED COMMON GENE SET HIGH-LEVEL
## CLUSTERS
cl.genes <- clust.genes(gsclust)
cl.genes <- cl.genes$cluster.genes %>% filter(cluster != 2)

## DETERMINE IF GENES IN ENRICHED COMMON GENE SET HIGH-LEVEL CLUSTERS ARE
## PRESENT DISPROPORTIONALLY IN THESE GENE-SETS THEN IN OTHERS (E.G. ARE THE
## GENES THAT ARE COMMON AMONGST THESE HIGH GENESETS THERE BECAUSE THEY
## DEFINE THE HIGH-LEVEL FUNCTIONS OR BECAUSE THEY ARE UBIQUITOUS IN ALL
## GENE SETS
gene.dists <- gene.dist.test(gene.gs.dist, cl.genes)


## CREATE A BINARY HEAT MAP FOR EACH HIGH-LEVEL CLUSTER WITH
## PROJECT ON THE X, AND GENE SET ON THE Y
plot.file <- paste0(ResultDir, "Plots/Cluster_genesets_HighLevel.pdf")
plot.cluster.proj.geneset(cl.genesets, gsea, file = plot.file)


## CREATE A BINARY HEAT MAP FOR EACH HIGH-LEVEL CLUSTER WITH
## PROJECT ON X AND GENE ON Y RANKED FROM MOST TO LEAST COMMON (TO GENE
## SETS IN CLUSTER

plot.file <- paste0(ResultDir, "Plots/Cluster_gene_DE_HighLevel.pdf")
cluster.proj.genes.in.gs(cl.genes.25, logPs.tcga, file = plot.file)




## BRING IN RESULTS OF GSEA ANALYSIS FOR GTEX ##
gsea.gtex <- get.fgsea.results(tissues, GSEADirGtex)

file <- paste0(PlotDir, "GTEx_qqplots.pdf")
get.median.p.value(gsea.gtex, file)

## HOW MANY GENE SETS HAVE Q-VALUES LESS THAN 0.05 using ALL CHROMOSOMES ##
t <- gsea.gtex %>% filter(padj <= 0.05, nox==0) %>% select(pathway, proj) %>% table()
sum(t)
mean(t)
tiss.gsea.count <- colSums(t)
median(colSums(t))
dim(t)
sort(table(rowSums(t)), decreasing=T)


## HOW MANY GENE SETS HAVE Q-VALUES LESS THAN 0.05 using AUTOSOMES ONLY ##
t <- gsea.gtex %>% filter(padj <= 0.05, nox==1) %>% select(pathway, proj) %>% table()
sum(t)
tiss.gsea.count.nox <- colSums(t)
median(colSums(t))
dim(t)
sort(table(rowSums(t)), decreasing=T)
t <- gsea.gtex %>% filter(padj <= 0.05, nox==1) %>% select(pathway, proj) %>% table()
sum(t)
colSums(t)
dim(t)

























## Preform quantile to quantile comparison ##
## What is the average rank of the genes TCGA in the xth quantile of GTEx ##
## And vice-versa
quants <- c(.001, .01, .05, .10, .15, .20, .25)
QQ.check <- quant.quant.comparo(proj.tis.match, logPs.tcga, logPs.gtex, geneInfo, quants)


QQ.check[QQ.check$quant == 0.001 & QQ.check$x.txt == "inclX",]

file = paste0(PlotDir, "quantile_TCGA_GTEX.pdf")
QQ.plot(QQ.check, file)


## COMPARE TWO BRAIN REGIONS ##
proj.tis.match.gtex <- rbind(c("Brain_Cortex", "Brain_Caudate_basal_ganglia"),
                         c("Brain_Cortex", "Lung"),
                         c("Brain_Caudate_basal_ganglia", "Lung"))


quants <- c(.001, .01, .05, .10, .15, .20, .25)
QQ.check.gtex <- quant.quant.comparo(proj.tis.match.gtex, logPs.gtex, logPs.gtex, geneInfo, quants)

file = paste0(PlotDir, "quantileGTEX_GTEX.pdf")
QQ.plot(QQ.check.gtex, file=file)

QQ.check.brain$label <-  paste(QQ.check.brain$ref, QQ.check.brain$comp, sep=":")


## EXAMINE THE DISTRIBUTION OF SURVIVAL P-VALUES FOR GENES FOUND TO BE DE ##
## 1) IN BOTH CANCER AND GTEX, 2) CANCER ONLY, 3) GTEX ONLY
## ALSO
## 1) PERFORM ABOVE USING ALL DE GENES, 2) GENES UPREGULATED IN MALES, 3) GENE UPREG IN FEMALES
## ALSO
## SURVAL PERFORMED IN 1) FULL SAMPLE, 2) MALES ONLY, 3) FEMALES ONLY
##
RsltDir <- "/gpfs/data/stranger-lab/askol/TCGA2/Survival/Results/"
PlotDir <- paste0(ResultDir, "Plots/")

rslts <- compile.surv.results(projects, RsltDir)

ps  <- rslts[[1]]
qs <- rslts[[2]]
coef <- rslts[[3]]
## remove genes with no symbols
for (i in 1:length(ps)){

    ps[[i]] <- filter(ps[[i]], !is.na(SYMBOL))
    qs[[i]] <- filter(qs[[i]], !is.na(SYMBOL))
    coef[[i]] <- filter(coef[[i]], !is.na(SYMBOL))
}

niter = 500
## DE RANKS BASED ON GTEX DE ANALYSIS ##
## quant = .001 ##
no.x.ind.p <- which(ps[[1]]$chr %in% c(1:22))
no.x.ind.gtex <- which(geneInfo.gtex$chr %in% c(1:22))
no.x.ind.tcga <- which(geneInfo.tcga$chr %in% c(1:22))

##################
##              ##
## quant = 0.001 ##
##              ##
##################
quant = 0.001
file <- paste0(PlotDir,"SurvPvals_gtex.tcga_",quant,".pdf")
surv.summary.001 <- analyze.ranked.survival(tissues, projects,
                                            ps, qs,
                                            logPs.gtex, logFC.gtex,  
                                            logPs.tcga, logFC.tcga,
                                            quant=0.001,
                                            niter=niter, both=FALSE,
                                            file)

## no X
file <- paste0(PlotDir,"SurvPvals_gtex.tcga_noX",quant,".pdf")
surv.summary.001 <- analyze.ranked.survival(tissues, projects,
                                            ps, qs,
                                            logPs.gtex[no.x.ind.gtex,], logFC.gtex[no.x.ind.gtex,],  
                                            logPs.tcga[no.x.ind.tcga,], logFC.tcga[no.x.ind.tcga,],
                                            quant=0.001,
                                            niter=niter, both=FALSE,
                                            file)

## DE RANKS BASED ON TCGA DE ANALYSIS ##
file <- paste0(PlotDir,"SurvPvals_tcga.tcga_",quant,".pdf")
surv.sum.tcga.tcga.001 <- analyze.ranked.survival(projects, projects,
                                                  ps, qs,
                                                  logPs.tcga, logFC.tcga,
                                                  logPs.tcga, logFC.tcga,
                                                  quant=0.001,
                                                  niter=niter, both = FALSE,
                                                  file)

## no X ##
file <- paste0(PlotDir,"SurvPvals_tcga.tcga_noX_",quant,".pdf")
surv.sum.tcga.tcga.001 <- analyze.ranked.survival(projects, projects,
                                                  ps, qs,
                                                  logPs.tcga[no.x.ind.tcga,], logFC.tcga[no.x.ind.tcga,],
                                                  logPs.tcga[no.x.ind.tcga,], logFC.tcga[no.x.ind.tcga,],
                                                  quant=0.001,
                                                  niter=niter, both = FALSE,
                                                  file)


##################
##              ##
## quant = 0.01 ##
##              ##
##################
quant = .01
## USE GENES THAT ARE ONLY DE IN BOTH GTEX AND TCGA ##
file <- paste0(PlotDir,"SurvPvals_gtex.tcga_",quant,".pdf")
surv.summary.01 <- analyze.ranked.survival(tissues, projects,
                                           ps, qs,
                                           logPs.gtex, logFC.gtex, 
                                           logPs.tcga, logFC.tcga,
                                           quant=quant,
                                           niter=niter, both=FALSE,
                                           file)

## NO X ##
file <- paste0(PlotDir,"SurvPvals_gtex.tcga_noX",quant,".pdf")
surv.summary.01 <- analyze.ranked.survival(tissues, projects,
                                           ps, qs,
                                           logPs.gtex[no.x.ind.gtex,], logFC.gtex[no.x.ind.gtex,], 
                                           logPs.tcga[no.x.ind.tcga,], logFC.tcga[no.x.ind.tcga,],
                                           quant=quant,
                                           niter=niter, both=FALSE,
                                           file)

file <- paste0(PlotDir,"SurvPvals_tcga.tcga_",quant,".pdf")
surv.sum.tcga.tcga.01 <- analyze.ranked.survival(projects, projects,
                                                 ps, qs,
                                                 logPs.tcga, logFC.tcga,
                                                 logPs.tcga, logFC.tcga,
                                                 quant=quant,
                                                 niter=niter, both=FALSE,
                                                 file)
## NO X##
file <- paste0(PlotDir,"SurvPvals_tcga.tcga_noX_",quant,".pdf")
surv.sum.tcga.tcga.01 <- analyze.ranked.survival(projects, projects,
                                                 ps, qs,
                                                 logPs.tcga[no.x.ind.tcga,], logFC.tcga[no.x.ind.tcga,],
                                                 logPs.tcga[no.x.ind.tcga,], logFC.tcga[no.x.ind.tcga,],
                                                 quant=quant,
                                                 niter=niter, both=FALSE,
                                                 file)

file <- paste0(PlotDir,"SurvPvals_gtex.tcga_inboth_",quant,".pdf")
surv.summary.DE.both.01 <- analyze.ranked.survival(tissues, projects,
                                                   ps, qs,
                                                   logPs.gtex, logFC.gtex,
                                                   logPs.tcga, logFC.tcga,
                                                   quant=quant,
                                                   niter=niter, both=TRUE,
                                                   file)
## NO X ##
file <- paste0(PlotDir,"SurvPvals_gtex.tcga_inboth_noX_",quant,".pdf")
surv.summary.DE.both.01 <- analyze.ranked.survival(tissues, projects,
                                                   ps, qs,
                                                   logPs.gtex[no.x.ind.gtex,], logFC.gtex[no.x.ind.gtex,],
                                                   logPs.tcga[no.x.ind.tcga,], logFC.tcga[no.x.ind.tcga,],
                                                   quant=quant,
                                                   niter=niter, both=TRUE,
                                                   file)

