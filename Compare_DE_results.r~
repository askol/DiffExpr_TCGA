ResultsDir <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/"
SummaryDir <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/Summary/"
GSEADir <- "/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Results/GSEA/"
source("/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Code/Compare_DE_results_funcs.r")

## GET NECESSARY ANALYSIS OUTPUT ##
compare.results <- function(tissue){

    SummaryDir <- paste0(SummaryDir, tissue,"/")
    if (!dir.exists(SummaryDir)){
        dir.create(SummaryDir)
    }
    setwd(SummaryDir)

    if (0){
        
    ## DETERMINE NUMBER OF MALE / FEMALE SUBJECTS
    Ns <- check.sample.size(tissue)
    
    rslts <- get.results(tissue, ResultsDir)

    ## GET RANKINGS FOR GSEA ANALYSIS ##
    gsea.scores <- make.gsea.scores(rslts)
    gsea.scores.nox <- make.gsea.scores(rslts[rslts$chr != "X",])

    RankDir <- paste0(GSEADir, gsub("-","_",tissue), "/Ranks/")
    if (!dir.exists(RankDir)){
        dir.create(RankDir, recursive=T)
    }

}
    GSEAOutDir <- paste0(GSEADir, gsub("-","_",tissue), "/")
if (0){
    RankFile <- paste0(RankDir, tissue, ".rnk")
    write.table(file = RankFile, gsea.scores, quote=F, row.names=F, col.names=F, sep="\t")
    
    RankFile <- paste0(RankDir, tissue, "_nox.rnk")
    write.table(file = RankFile, gsea.scores.nox, quote=F, row.names=F, col.names=F, sep="\t")
    
    ## RUN GSEA ##
    run.gsea(tissue, RankDir, GSEAOutDir,
             GSDir = "/home/askol/bin/GSEA_genesets/msigdb_v6.2/",
             geneSet = "msigdb.v6.2.symbols.gmt")

    run.gsea(tissue, RankDir, GSEAOutDir,
             GSDir = "/home/askol/bin/GSEA_genesets/Custom/",
             geneSet = "Hormone_Immune_Custom.gmt")

}
    process.gsea.output(tissue, GSEAOutDir, "msigdb.v6.2.symbols.gmt", qThresh=0.25)
    process.gsea.output(tissue, GSEAOutDir, "Hormone_Immune_Custom.gmt", qThresh=0.25)
    
}

## ---- MAIN -----##

source("/gpfs/data/stranger-lab/askol/GTEx/DiffExpression/Code/Compare_DE_results_funcs.r")
args <- commandArgs(TRUE)
tissue = args[1]
compare.results(tissue)



