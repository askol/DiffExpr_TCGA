library(edgeR)
library(limma)
library(dplyr)
library(survival)
library(survminer)
library(WGCNA)
library(pscl)
library(reshape2)
library(fgsea, lib.loc="/home/askol/R/x86_64-pc-linux-gnu-library/3.5")


get.results <- function(project, ResultsDir){

    file <- paste0(ResultsDir,project,"_DE_voom.rslts")
    rslts <- read.table(file = file, as.is=T, header=T, quote="", sep="\t")

    return(rslts)
}

make.gsea.scores <- function(rslts){


    scores <- data.frame(gene = rslts$SYMBOL)
    
    ## ADJUST FOR INFINITIES # BREAK TIES WITH ABS(LOGFC)
    inf.ind <- which(is.infinite(rslts$P.Value))
    if (length(inf.ind) > 0){
        rnk <- rank(-1*abs(rslts$logFC[inf.ind]))
        mx <- max(rslts$P.value[-inf.ind])
        rslts$P.value[inf.ind] <- max + rnk
    }
    
    scores$score <- sign(rslts$logFC)*(-log10(rslts$P.Value))
    ## AVERAGE SCORES OF GENES THAT HAVE MORE THAN ONE ENTRY ##
    scores <- scores %>% group_by(gene) %>%
        summarize(score = mean(score, na.rm = TRUE))

    scores <- arrange(scores, desc(score))
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

run.fgsea <- function(project, gsea.scores, msig.file, out.file){
    
    pathways <- gmtPathways(msig.file)   
    scores <- unlist(gsea.scores$score);
    names(scores) <- gsea.scores$gene

    print("Running fgseaRes. . .")
    fgseaRes <- fgseaMultilevel(pathways = pathways, 
                      stats = scores,
                      minSize=15,
                      maxSize=500)

    fgseaRes <- as_tibble(fgseaRes)
    fgseaRes$Ncontrib <- sapply(fgseaRes$leadingEdge, length)

    print("Calculating tags, list and signal. . .")
    fgseaResEx <- calc.tags(fgseaRes, scores)
    
    fgseaResEx$LEgenes <- sapply(fgseaResEx$leadingEdge, paste, collapse=";")
    fgseaResEx <- select(fgseaResEx, -leadingEdge)
    fgseaResEx <- arrange(fgseaResEx, padj, pval)

    write.table(file = out.file, fgseaResEx, quote=F, row.names=F, col.names=T)
    print(paste0("Wrote fgsea results to ",out.file))
    
    return(fgseaRes)    
}


calc.tags <- function(fgseaRes, scores){
    
    scores <- sort(scores)
    gnames <- factor(names(scores), levels=names(scores))
    
    o <- apply(fgseaRes, 1, function(x) tags.list.sig(x, gnames))
    o <- do.call(rbind, o)

    fgseaRes <- as_tibble(cbind(fgseaRes, o))
    return(fgseaRes)
}

tags.list.sig <- function(x, gnames){
  
    n <- x$size
    N <-length(gnames)
    legenes <- factor(unlist(x$leadingEdge), levels = levels(gnames))
    if (sign(x$ES) == -1){
        n.peak <- max(as.numeric(legenes))
        le.n.med <- median(as.numeric(legenes))
    }else{
        n.peak <- N-min(as.numeric(legenes))
         le.n.med <- N - median(as.numeric(legenes))
    }
    tags <- length(legenes)/n
    lst <- n.peak/N
    signal <- tags*(1 - lst)*N/(N-n)

    return(data.frame(tags = tags, lst = lst, signal = signal, n.peak, le.n.med))
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


make.fgsea.pathway <- function(msig.file){

    out.file <- paste0(msig.file,".RDS")
    if (file.exists(out.file)){

        out <- readRDS(out.file)
        return(out)
    }else{

        d <- readLines(msig.file)
        d <- lapply(d, function(x){gs = gsub("\t.+","", x);
                                   genes = unlist(strsplit(x,split="\t"))[-c(1:2)];
                                   eval(parse(text = paste("list(",gs," = genes)")))})
        saveRDS(d, out.file)
        return(d)
    }
        
}        
    
    

get.msigdb.genes <- function(msig.file){

    d <- read.table(file = msig.file, as.is=T, header=F, fill=T)

    d <- d[,-c(1,2)]
    d <- unlist(c(d))
    d <- unique(d)
    ind <- which(d == "")
    d <- d[-ind]
    return(d)
}
