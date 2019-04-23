##  Create QQ plots, Volcano Plots, HeatMap of sign of logFC,
##  tSNPE or PCA Cluster based on sign of log FC
##  SHARED SIGNIFICANT GSs
library(qvalue)
library(ggpubr)

source("/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Code/Summary_funcs.r")

GSEADir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Results/GSEA/"
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


#####  PLOT 1: NUMBER OF DE GENES ######
## GET DISTRIBUTION OF THE NUMBER OF DE GENES ##
get.noDE.genes(Qs.tcga, Qs.gtex, project, tissues, PlotDir)

## WHICH GENES ARE COMMONLY DE BETWEEN GTEX AND TCGA ##


## ## WHICH GENES ARE SHARED ACROSS CANCER
## META ANALYSIS AND FDR > .10
## PLOT NUMBER OF CANCERS/TISSUES SHARED PER GENE
## AND
## PLOT NUMBER GENES FROM EACH CHROMOSOME WITH MULTIPLE CANCERS/TISSUES/GENE

file.pre = "tcga"
noDE.tcga <- get.noDE.cancers.for.gene(Qs.tcga, qthresh= .1, geneInfo.tcga,
                                       Ntrunc = 5, min.studies = 5,
                                       ResultDir, PlotDir, file.pre)
file.pre = "gtex"
noDE.gtex <- get.noDE.cancers.for.gene(Qs.gtex, qthresh= .1, geneInfo.gtex,
                                       Ntrunc = 5, min.studies = 10,
                                       ResultDir, PlotDir, file.pre)

## PLOT THE DISTRIBUTION OF DE GENES ON THE X CHROMOSOME ##
x.tcga <- plot.de.on.x(logPs.tcga, logFC.tcga, geneInfo.tcga, PlotDir, file.pre="tcga")
x.gtex <- plot.de.on.x(logPs.gtex, logFC.gtex, geneInfo.gtex, PlotDir, file.pre="gtex")
    

## FIND OUT WHICH GENES ARE COMMONLY DE IN BOTH CANCER AND NORMAL TISSUE ##
tmp <- get.common.common.DE.genes(noDE.gtex, noDE.tcga, Qs.tcga, geneInfo.tcga,
                           Qs.gtex, geneInfo.gtex,
                           thresh.gtex = 10, thresh.tcga=5,
                           ResultDir, PlotDir)
Ns <- tmp$N; Qs <- tmp$Qs

Ns.dist <- get.n.dist(Ns)

file <- paste0(PlotDir,"Dist_DE_Conditional.pdf")
plot.cond.de.distribution(Ns.dist, file=file)

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

