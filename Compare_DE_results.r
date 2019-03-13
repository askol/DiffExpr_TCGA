ResultsDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Results/"
SummaryDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Summary/"
GSEADir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Results/GSEA/"
source("/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Code/Compare_DE_results_funcs.r")

## GET NECESSARY ANALYSIS OUTPUT ##
compare.results <- function(project){
    
    SummaryDir <- paste0(SummaryDir, project,"/")
    if (!dir.exists(SummaryDir)){
        dir.create(SummaryDir)
    }
    setwd(SummaryDir)
    
    ## DETERMINE NUMBER OF MALE / FEMALE SUBJECTS
    ##Ns <- check.sample.size(project)
    
    rslts <- get.results(project, ResultsDir)
    if (length(grep("SYMBOL", names(rslts))) == 0){
        rslts$SYMBOL = rslts$gene
    }

    ## GET RANKINGS FOR GSEA ANALYSIS ##
    gsea.scores <- make.gsea.scores(rslts)
    gsea.scores.nox <- make.gsea.scores(rslts[rslts$chr != "X",])

    RankDir <- paste0(GSEADir, gsub("-","_",project), "/Ranks/")
    if (!dir.exists(RankDir)){
        dir.create(RankDir, recursive=T)
    }
           
    GSEAOutDir <- paste0(GSEADir, gsub("-","_",project), "/")
    
    RankFile <- paste0(RankDir, project, ".rnk")
    write.table(file = RankFile, gsea.scores, quote=F, row.names=F, col.names=F, sep="\t")
    
    RankFile <- paste0(RankDir, project, "_nox.rnk")
    write.table(file = RankFile, gsea.scores.nox, quote=F, row.names=F, col.names=F, sep="\t")
    
    ## RUN GSEA ##
    run.gsea(project, RankDir, GSEAOutDir,
             GSDir = "/home/askol/bin/GSEA_genesets/msigdb_v6.2/",
             geneSet = "msigdb.v6.2.symbols.gmt")
    
    run.gsea(project, RankDir, GSEAOutDir,
             GSDir = "/home/askol/bin/GSEA_genesets/Custom/",
             geneSet = "Hormone_Immune_Custom.gmt")
    
    process.gsea.output(project, GSEAOutDir, "msigdb.v6.2.symbols.gmt", qThresh=0.25)
    process.gsea.output(project, GSEAOutDir, "Hormone_Immune_Custom.gmt", qThresh=0.25)
    
}

## ---- MAIN -----##

source("/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Code/Compare_DE_results_funcs.r")
args <- commandArgs(TRUE)
project = args[1]
compare.results(project)



