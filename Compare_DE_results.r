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

    rslts <- select(rslts, -gene)
    
    ## REMOVE GENES W/O NAMES ##
    rslts <- filter(rslts, !is.na(SYMBOL))

    ## REMOVE MULTIMAPPED GENES ##
    dupeGenes <- unique(rslts$SYMBOL[duplicated(rslts$SYMBOL)])

    for (gene in dupeGenes){

        t <- dplyr::filter(rslts, SYMBOL == gene) %>% select(chr) %>% table()
        r <- filter(rslts, SYMBOL == gene) %>% select(start) %>% range(na.rm=T)
        r <- abs(diff(r))
        ## remove gene if on more than one chromosome or different start sites are
        ## > 500k apart ##
        if (length(t) > 1 | r > 500000){
            rslts <- filter(rslts, SYMBOL != gene)
        }
    }
   
    ## GET RANKINGS FOR GSEA ANALYSIS ##    
    gsea.scores <- make.gsea.scores(rslts)
    gsea.scores.nox <- make.gsea.scores(rslts[rslts$chr %in% c("X","Y") == FALSE,])

           
    GSEAOutDir <- paste0(GSEADir, gsub("-","_",project), "/")

    if (0){
        
        RankDir <- paste0(GSEADir, gsub("-","_",project), "/Ranks/")
        if (!dir.exists(RankDir)){
            dir.create(RankDir, recursive=T)
        }
        
        RankFile <- paste0(RankDir, project, ".rnk")
        write.table(file = RankFile, gsea.scores, quote=F, row.names=F, col.names=F, sep="\t")
        
        RankFile <- paste0(RankDir, project, "_nox.rnk")
        write.table(file = RankFile, gsea.scores.nox, quote=F, row.names=F, col.names=F, sep="\t")
    }
    
    ## CREATE MSIG GENE SET OBJECT FOR FGSEA ##
    msig.file <- "/home/askol/bin/GSEA_genesets/msigdb_v6.2/msigdb.v6.2.symbols.gmt"
    out.file <- paste0(GSEAOutDir,project,"_fGSEA_summary.txt")
    rslts <- run.fgsea(project, gsea.scores, msig.file, out.file)

    out.file <- paste0(GSEAOutDir,project,"_fGSEA_noX_summary.txt")
    rslts.nox <- run.fgsea(project, gsea.scores.nox, msig.file, out.file)
        
}

if (0){
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



