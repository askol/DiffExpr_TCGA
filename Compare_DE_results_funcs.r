library(edgeR)
library(limma)
library(dplyr)
library(survival)
library(survminer)
library(WGCNA)
library(pscl)
library(reshape2)


get.results <- function(project, ResultsDir){

    file <- paste0(ResultsDir,project,"_DE_limma.rslts")
    rslts <- read.table(file = file, as.is=T, header=T, quote="", sep="\t")

    return(rslts)
}

make.gsea.scores <- function(rslts){


    scores <- data.frame(gene = rslts$SYMBOL)
    
    ## scores$score <- sign(rslts$logFC - 1)*(-log10(rslts$P.Value))
    scores$score <- sign(rslts$logFC)*(-log10(rslts$P.Value))
    scores <- scores[order(scores$score, decreasing=T),]
    return(scores)
    
}

summarize.gsea <- function(x, col, adj.N=""){

    ## METRICS FOR TFS
    col.ind <- which(names(x) == col)

    ## arac data has TF named Regulator ##
    names(x) <- gsub("Regulator","TF", names(x))    
    
    stats.tf <- calc.stats(x[, col.ind], x$TF)
    stats.targ <- calc.stats(x[, col.ind], x$Target)

    rm.ind <- which(stats.targ$gene %in% stats.tf$gene)
    if (length(rm.ind)>0){
        stats.targ <- stats.targ[-rm.ind,]
    }
    
    stats.tf$type <- "TF"
    stats.targ$type <- "Target"
    
    stats <- rbind(stats.tf, stats.targ)
    
    ## CREATE A SIGNED RANK BASED ON P.ADJ ##
    ## Ranks for p.max ##
    stats$max.rank <- sign.rank(stats$max.p.adj)
    stats$min.rank <- sign.rank(stats$min.p.adj)
    stats$any.rank <- sign.rank(
        (1 - 2*(stats$min.p.adj<stats$max.p.adj)) *
        apply(stats[,c("min.p.adj","max.p.adj")], 1, min))
    stats$mean.rank <- sign.rank(sign(stats$mean) * stats$mean.p.adj)

       return(stats)                        
}


run.gsea <- function(project, RankDir, GSEAOutDir, GSDir, geneSet){

    ## geneSet  =  "msigdb.v6.2.symbols.gmt" or "Hormone_Immune_Custom.gmt"
    RankFiles <- dir(pattern=".rnk", RankDir)                     

    geneSetName <- gsub("\\.v6.+","", geneSet)
    geneSetName <- gsub("\\.gmt", "", geneSet)

    for (RankFile in RankFiles){ ## one includes x genes, the other not

        print(paste("Working on ",geneSet," and ",RankFile))
        
        rpt.lab <- geneSetName
        if (length(grep("nox", RankFile)==1) ){
            rpt.lab <- paste0(rpt.lab,"_nox")
        }
            
        params <- rbind(c("rnk",
                          paste0(RankDir, RankFile)),
                        c("out",
                          paste0(GSEAOutDir, rpt.lab)))
        param.file <- paste0("~/GSEA_parameters_", gsub("-","_",project), ".txt")
        write.table(file = param.file, params, quote=F, col.names=F,
                    row.names=F, sep="\t")
        cmd <- paste0("java -cp /home/askol/bin/gsea-3.0.jar -Xmx5000m ",
                      "xtools.gsea.GseaPreranked ",
                      "-param_file ", param.file,
                      " -gmx ", GSDir, geneSet,
                      " -norm meandiv -nperm 1000 ",
                      "-scoring_scheme weighted -rpt_label ", rpt.lab,
                      " -create_svgs false -make_sets true -plot_top_x 100 ",
                      "-rnd_seed timestamp ",
                      "-set_max 500 -set_min 15 -zip_report false  -gui false")
        print(cmd)
        system(cmd)
    }
    
}
    
process.gsea.output <- function(project, GSEAOutDir,
                                geneSet = "msigdb.v6.2.symbols.gmt", qThresh=0.25){

    ## READ IN RESULTS AND SORT BY NORMALIZED ES ##

 
    for (suff in c("","_nox")){

          gsea <- c()
          
          print(paste0("Working on project: ",project, " with ",suff))

          for (sign in c("pos","neg")){
              
              geneSetName <- gsub("\\.gmt", "", geneSet)
              
              out.dir <- paste0(GSEAOutDir, geneSetName,suff,"/")
              
              dirs <- list.dirs(out.dir)
              dir <- most.recent(dirs[-1])
              files <- dir(dir, pattern="gsea")
              files <- files[grep("xls",files)]
              file <- paste0(dir,"/",files[grep(sign,files)])
              if (length(files) == 0){
                  print(paste0("Didn't find file: ",file))
                  next
              }
              out <- read.table(file = file, as.is=T, header=T, sep="\t")
              if (nrow(out) == 0){
                  print(paste0("No data found in file ",file))
                  next
              }
              nr <- nrow(out)
              qpass.ind <- which(out$FDR.q.val < qThresh)
              if (length(qpass.ind)>0){
                  gsea <- rbind(gsea, cbind(geneSetName, out[qpass.ind,]))
              }
          }        
          
          ## WRITE OUT RESULTS ##
          out.file <- paste0(GSEAOutDir, project, "_", geneSet,suff,"_GSEA_summary.txt")
          
          ## RETURN EMPTY FILE IS NO GS PASS QTHRESH ##
          if (is.null(dim(gsea))){
              system(paste0("touch ",out.file))
              print(paste0("No gene sets passed FDR q-value threshold of ",qThresh))           
          }else{
              
              ## REMOVE UNWANTED COLUMNS AND ROWS ##
              rm.cols <- which(names(gsea) %in% c("GS.DETAILS","X"))
              gsea <- gsea[, -rm.cols]
              
              ## REMOVE SETS WITH FDR.P.VAL == NA ##
              rm.ind <- which(is.na(gsea$FDR.q.val))
              if (length(rm.ind) > 0){
                  gsea <- gsea[-rm.ind,]
              }   
              
              gsea <- gsea[order(gsea$FDR.q.val),]
              
              write.table(file = out.file, gsea, quote=F, row.names=F, col.names=T, sep="\t")
              print(paste0("Wrote GSEA summary in ",out.file))
          }
      }
}

  
most.recent <- function(dirs){

    Ds <- c()
    for (dir in dirs){
        
        i <- file.info(dir)$mtime

        Ds <- rbind(Ds, c(dir, i, as.character(i)))
    }

    ord <- order(Ds[,2], decreasing=T)

    return(Ds[ord[1], 1])
}

    
        
    
