## FUNCTIONS TO HELP SUMMARIZE THE RESULTS FROM GTEX
## DE ANALYSIS

DataDir <- "/gpfs/data/stranger-lab/askol/TCGA2/Coexpression/Data/"

source("~/Code/qq_plot.r")
library(calibrate)
library(WGCNA) ## for the faster cor function
library(reshape2)
library(RColorBrewer)
library(qvalue)
library(dplyr)
library(annotables)


collect.results <- function(projects, ResultDir){
    
    logFC <- c()
    logPs <- c()
    logQs <- c()
    geneInfo <- c()
    for (project in projects){
        
        print(paste0("Working on project : ",project))
        rslt <- load.results(project, ResultDir)
        dupe.ind <- which(duplicated(rslt$SYMBOL))
        if (length(dupe.ind) > 0){
            rslt <- rslt [-dupe.ind,]
        }
        print(paste0("Removing ",length(dupe.ind), " duplicate gene"))
        logFCtmp <- rslt[,c("SYMBOL", "logFC")]
        logPstmp <- rslt[,c("SYMBOL", "P.Value")]
        names(logFCtmp)[2] = names(logPstmp)[2] = project        

        gi <-  rslt[, c("ENSEMBL","SYMBOL", "chr", "start","end","biotype")]
        if (project == projects[1]){
            logFC <- logFCtmp
            logPs <- logPstmp
            geneInfo <- gi
        }else{            
            logFC <- merge(logFC, logFCtmp, by="SYMBOL", all=T)
            logPs <- merge(logPs, logPstmp, by="SYMBOL", all=T)
            geneInfo <- rbind(geneInfo, gi[!gi$SYMBOL %in% geneInfo$SYMBOL,])
        }
    }

    ## CREATE FDR ##
    Qs <- make.Qs(logPs)

    ## take log10 of Ps
    logPs[,-1] = -log10(logPs[,-1])

    geneInfo <- geneInfo[match(logFC$SYMBOL, geneInfo$SYMBOL),]
    return(list(logFC = logFC, logPs = logPs, Qs = Qs, GeneInfo = geneInfo))
}


collect.results.gtex <- function(tissues, ResultDir){

    make.Qs <- function(logPs){
        Qs <- logPs
        for (i in 2:ncol(logPs)){
            p <- logPs[,i]
            miss.ind <- which(is.na(p))
            p <- p[-miss.ind]
            small.ind<- which(p < 10^-30)
            p[small.ind] <- 10^-30
            q <- qvalue(p)$qvalues
            Qs[-miss.ind,i] <- q
        }
        return(Qs)
    }

    load.results <- function(tissue, ResultDir){
        
        file <- paste0(ResultDir, tissue, "_DE_limma.rslts")
        
        if (!file.exists(file)){
            print(paste0("Results file", file, "not found!"))
            return(0)
        }
        rslt <- read.table(file = file, header=T, as.is=T)
        
        return(rslt)
    }
    
    logFC <- c()
    logPs <- c()
    logQs <- c()
    geneInfo <- c()

    for (tissue in tissues){
        
        print(paste0("Working on tissue : ",tissue))
        rslt <- load.results(tissue, ResultDir)
        rslt <- filter(rslt, !duplicated(rslt$ensgene))

        logFCtmp <- rslt[,c("ensgene", "logFC")]
        logPstmp <- rslt[,c("ensgene", "P.Value")]
        names(logFCtmp)[2] = names(logPstmp)[2] = tissue        

        gi <-  select(rslt, ensgene, symbol, chr, start, biotype)
        if (tissue == tissues[1]){
            logFC <- logFCtmp
            logPs <- logPstmp
            geneInfo <- gi
        }else{            
            logFC <- merge(logFC, logFCtmp, by="ensgene", all=T)
            logPs <- merge(logPs, logPstmp, by="ensgene", all=T)
            geneInfo <- rbind(geneInfo, gi[!gi$ensgene %in% geneInfo$ensgene,])
        }
    }

    ## CREATE FDR ##
    Qs <- make.Qs(logPs)

    ## take log10 of Ps
    logPs[,-1] = -log10(logPs[,-1])

    geneInfo <- geneInfo[match(logFC$ensgene, geneInfo$ensgene),]
    return(list(logFC = logFC, logPs = logPs, Qs = Qs, GeneInfo = geneInfo))
}


load.results <- function(project, ResultDir){

    file <- paste0(ResultDir, project, "_DE_voom_nosvs.rslts")

    if (!file.exists(file)){
        print(paste0("Results file", file, "not found!"))
        return(0)
    }
    rslt <- read.table(file = file, header=T, as.is=T, sep="\t", quote="")
    if (length(grep("SYMBOL",names(rslt))) == 0){
        rslt$SYMBOL <- rslt$gene
    }
    return(rslt)
}

get.noDE.genes <- function(Qs.tcga, Qs.gtex, projects, tissues, PlotDir){
    
    ## GET TCGA SAMPLE SIZES ##
    Ns.tcga <- get.sample.sizes(projects)
    Ns.gtex <- get.gtex.sample.sizes(tissues)   

    ## TOTAL NUMBER OF DE GENES
    Nsig.DE.tcga <-  colSums(Qs.tcga[,-1] < .05, na.rm=T)
    Nsig.DE.tcga <- data.frame(project=names(Nsig.DE.tcga), NSig = Nsig.DE.tcga)
    Nsig.DE.tcga <- merge(Nsig.DE.tcga, Ns.tcga, by=1)
    Nsig.DE.tcga <- mutate(Nsig.DE.tcga, N.m = as.numeric(as.character(N.m)),
                           N.f = as.numeric(as.character(N.f)),
                           N.ave = (N.m + N.f)/2) %>%
                               arrange(desc(N.ave))
    file <- paste0(ResultDir, "NDE_tcga.txt")
    write.table(file = file, Nsig.DE.tcga, quote=F, row.names=F, col.names=T)
    
    Nsig.DE.gtex <-  colSums(Qs.gtex[,-1] < .05, na.rm=T)
    Nsig.DE.gtex <- data.frame(project=names(Nsig.DE.gtex), NSig = Nsig.DE.gtex)
    Nsig.DE.gtex <- merge(Nsig.DE.gtex, Ns.gtex, by=1)
    Nsig.DE.gtex <- mutate(Nsig.DE.gtex, N.m = as.numeric(as.character(N.m)),
                           N.f = as.numeric(as.character(N.f)),
                           N.ave = (N.m + N.f)/2) %>%
                               arrange(desc(N.ave))
    file <- paste0(ResultDir, "NDE_gtex.txt")
    write.table(file = file, Nsig.DE.gtex, quote=F, row.names=F, col.names=T)

    file = paste0(PlotDir,"NoDEGenes_tcga.pdf")
    Ns.tcga <- plotNDE(Nsig.DE.tcga, file = file)

    file = paste0(PlotDir, "NoDEGenes_tcga_trunc.pdf")
    trash <- plotNDE(Nsig.DE.tcga, file=file, max.count = 2000)
    
    file = paste0(PlotDir,"NoDEGenes_gtex.pdf")
    Ns.gtex <- plotNDE(Nsig.DE.gtex, file = file)

}

plotNDE <- function(Ns, file, max.count = ""){

    pdf(file, width = 12, height=8)
    
    ## ADJUST MAX COUNT IF REQUESTED ##
    if (max.count != ""){
        Ns <- mutate(Ns, NSig = ifelse(NSig > max.count, max.count, NSig))
    }
    
    rng.de <- range(Ns$NSig)
    rng.N <- range(c(Ns$N.m, Ns$N.f))
    Ns$N.m.sc <- diff(rng.de)/diff(rng.N)*(Ns$N.m - rng.N[1]) + rng.de[1]
    Ns$N.f.sc <- diff(rng.de)/diff(rng.N)*(Ns$N.f - rng.N[1]) + rng.de[1]
    Ns$N.ave.sc <- diff(rng.de)/diff(rng.N)*(Ns$N.ave - rng.N[1]) + rng.de[1]
    mult.fact <- max(c(Ns$N.m, Ns$N.f))/max(Ns$NSig)
    
    ind <- which(names(Ns) == "tissue")
    if (length(ind)>0){ names(Ns)[ind] = "project" }

    ord <- order(Ns$N.ave, decreasing=T)
    Ns$project <- factor(Ns$project, levels=Ns$project[ord])

 
    
    p <- ggplot(Ns, aes(x = project)) +
                geom_bar(aes(y=NSig), stat="identity", fill="steelblue") +
                geom_point(aes(y=N.m/mult.fact), color = "blue", size=2.5)+
                geom_point(aes(y=N.f/mult.fact), color = "red", size=2.5)+
                scale_y_continuous(sec.axis = sec_axis(~.*mult.fact,
                                       name = "Sample Size")) +
                theme_minimal() +
                theme(text=element_text(size=10),
                      axis.text.x=element_text(angle=60, hjust=1)) +
                          ylab("Number Sex DE Genes (q-value < 0.05)") +
                              xlab("")
    
    print(p)
    
    dev.off()

    return(Ns)
}

get.noDE.cancers.for.gene <- function(Qs, qthresh, geneInfo,
                                      Ntrunc, min.studies, 
                                      ResultDir, PlotDir, file.pre){

    NoDE <- Qs
    names(NoDE) <- gsub("ensgene", "ENSEMBL",names(NoDE))
    names(geneInfo) <- gsub("ensgene", "ENSEMBL",names(geneInfo))
    if (any(names(NoDE)=="SYMBOL")){
        NoDE$gene <- NoDE$SYMBOL
        geneInfo$gene <- geneInfo$SYMBOL
    }else {
        NoDE$gene <- NoDE$ENSEMBL
        geneInfo$gene <- geneInfo$ENSEMBL
    }
   
    NoDE$n <- rowSums(NoDE[,-1] <= qthresh, na.rm=T)
    NoDE$N <- rowSums(!is.na(NoDE[,-1]))
    
    NoDE <- select(NoDE, gene, n, N)
  
    NoDE <- right_join(geneInfo, NoDE, by="gene") %>% arrange(desc(n))
    
    ## distribution of number of studies with same de gene ##
    file = paste0(ResultDir,file.pre,"_nodestudiespgene.txt")
    write.table(file = file, NoDE, quote=F, row.names=F, col.names=T)

    file = paste0(PlotDir,file.pre,"_nodestudiespgene.pdf")
    plot.nostudies(NoDE, file, trunc=Ntrunc, min.studies=min.studies)        
    
    return(NoDE)
}

get.common.common.DE.genes <- function(N.gtex, N.tcga, Qs.tcga, geneInfo.tcga,
                                       Qs.gtex, geneInfo.gtex,
                                       thresh.gtex = 10, thresh.tcga=5,
                                       ResultDir, PlotDir){

    gi <- grch38[,c("symbol","ensgene","chr","start")]
    N <- select(N.gtex, ENSEMBL, n, N)
    names(N) <- gsub("n$","n.gtex", names(N))
    names(N) <- gsub("N$","N.gtex", names(N))
    N <- left_join(N, gi, by = c("ENSEMBL"="ensgene"))
    N <- inner_join(N, N.tcga[,c("ENSEMBL","n","N")], by="ENSEMBL")
    names(N) <- gsub("n$","n.tcga", names(N))
    names(N) <- gsub("N$","N.tcga", names(N))

    ## WHICH GENES HAVE MORE THAN thresh.gtex and thresh.tcga DE studies ##
    N <- mutate(N, de.both = ifelse(n.gtex >= thresh.gtex & n.tcga >= thresh.tcga, 1, 0))
    N.out <- N %>% filter(de.both == 1) %>% select(symbol, chr, start, n.tcga, N.tcga,
                              n.gtex, N.gtex) %>% arrange(desc(n.tcga))
    file <- paste0(ResultDir, "DE_in_TCGAge5_GTEX_gte10.txt")
    write.table(file = file, N.out, quote=F, row.names=F, col.names=T)
    print(paste0("Wrote ",file))
    
    ## FIND GENES COMMONLY DE IN TCGA BUT NOT GTEX AND VISE VERSA ##
    N <- N %>% mutate(de.tcga.n.gtex = ifelse(n.tcga >= thresh.tcga & n.gtex <= 2, 1, 0))
    N <- N %>% mutate(de.gtex.n.tcga = ifelse(n.gtex >= thresh.gtex & n.tcga <= 2, 1, 0))
    N.out <-  N %>% filter(de.tcga.n.gtex == 1 | de.gtex.n.tcga==1) %>%
        select(symbol, chr, start, n.tcga, N.tcga,
                              n.gtex, N.gtex) %>% arrange(desc(n.tcga))
    file <- paste0(ResultDir, "DE_different_TCGAge5_GTEX_gte10.txt")
    write.table(file = file, N.out, quote=F, row.names=F, col.names=T)
    print(paste0("Wrote ",file))

    ## Are genes that are commonly DE across cancers also commonly DE
    ## expressed across GTEx tissues
    ## REQUIRE OBSERVATIONS IN AT LEAST 80% OF STUDIES
    n.max.gtex <- max(N$N.gtex); n.max.tcga <- max(N$N.tcga)
    N <- mutate(N, incl = ifelse(N.gtex >= 0.80*n.max.gtex & N.tcga >= 0.80*n.max.tcga &
                       (n.gtex > 0 | n.tcga>0),
                       1, 0))
    N <- mutate(N, n.gtex.imp = n.gtex*N.gtex/n.max.gtex) %>%
        mutate(n.tcga.imp = n.tcga*N.tcga/n.max.tcga)

    reg <- lm(n.gtex.imp ~ n.tcga.imp + I(n.tcga.imp*n.tcga.imp) + I(n.tcga.imp^3), data = subset(N, incl==1))
    N$stdR <- NA
    N$stdR[N$incl==1] <- stdres(reg)

    
    ## CHECK CORRELATION OF SIGNIFICANCE BETWEEN CANCERS AND BETWEEN CANCERS AND GTEX
    Qs <- left_join(Qs.gtex, gi, by = "ensgene")
    Qs <- inner_join(Qs, Qs.tcga, by=c("symbol"="SYMBOL"))
    tcga.names <- names(Qs.tcga)[-1]
    gtex.names <- names(Qs.gtex)[-1]

    file <- paste0(PlotDir, "DE_p_corr_TCGA.pdf")
    corr.tcga <- plot.de.corr(Qs, names1=tcga.names, file=file)
    file <- paste0(PlotDir, "DE_p_corr_GTEX.pdf")
    corr.gtex <- plot.de.corr(Qs, gtex.names, file=file)
    file <- paste0(PlotDir, "DE_p_corr_TCGA_GTEX.pdf")
    corr.tcga.gtex <- plot.de.corr(Qs, tcga.names, gtex.names, file=file)

    file <- paste0(ResultDir, "DE_p_corr_TCGA_best.txt")
    corr.tcga.best.match <- best.match(corr.tcga, n=3)
    write.table(file = file, corr.tcga.best.match, quote=F, row.names=F, col.names=T)

    file <- paste0(ResultDir, "DE_p_corr_GTEX_best.txt")
    corr.gtex.best.match <- best.match(corr.tcga, n=3)
    write.table(file = file, corr.gtex.best.match, quote=F, row.names=F, col.names=T)

    file <- paste0(ResultDir, "DE_p_corr_TCGA_GTEX_best.txt")
    corr.tcga.gtex.best.match <- best.match(corr.tcga.gtex, n=3)
    write.table(file = file, corr.tcga.gtex.best.match, quote=F, row.names=F, col.names=T)

    
    return(list(N=N, Qs = Qs))

}

best.match <- function(C, n){

    best <- c()
    type <- ""
    for (i in 1:2){
        if (i==1){ corr <- C$corr; type="all" }else{ corr = C$corr.x; type="X" }
        corr <- mutate(corr, study=as.character(study)) %>%
            mutate(study2 = as.character(study2))
        for (proj in unique(corr$study)){
            
            tmp <- filter(corr, study==proj) %>% filter(study != study2) %>%
                mutate(rnk = rank(-cor)) %>%
                    filter(rnk <= n) %>% arrange(rnk) %>% select(-rnk)
            best <- rbind(best, c(type, as.character(tmp$study[1]),
                                  as.character(tmp$study2), tmp$cor))
        }
    }
    colnames(best) <- c("genome","study", "best.match", "best.match2",
                        "best.match2", "cor", "cor","cor")
    return(best)
}

quant.quant.comparo <- function(proj.tis.match, logPs1, logPs2, geneInfo, quants){
    
    QQ.check <- c()
    N.genes <-c()
    for (i in 1:nrow(proj.tis.match)){

        c <- merge(logPs1[,c("SYMBOL",proj.tis.match[i,1])],
                   logPs2[,c("SYMBOL",proj.tis.match[i,2])],
                   by="SYMBOL", all.x=F, all.y=F)
        
        N.genes <- rbind(N.genes, c(proj.tis.match[i,1], proj.tis.match[i,2],
                                    sum(!is.na(c[,2])),sum(!is.na(c[,3])),
                                    sum(rowSums(is.na(c))==0)))
        rm <- which(rowSums(is.na(c)) > 0)
        c <- c[-rm,]
        
        x.genes <- geneInfo$SYMBOL[geneInfo$chr %in% c("X","Y")]
        x.ind <- which(c$SYMBOL %in% x.genes )
        c.nox <- c[-x.ind,]
        
        x.txt <- "inclX"
        for (j in 1:2){
            if (j == 2){ c = c.nox; x.txt = "noX"}
            cutoffs.1 <- quantile(c[,2], probs=1-quants)
            cutoffs.2 <- quantile(c[,3], probs=1-quants)
            
            mean.quant.1<- get.quant(c, ref.col = 2, comp.col=3, cutoffs.1, quants)
            mean.quant.2 <- get.quant(c, ref.col = 3, comp.col=2, cutoffs.2, quants)
            
            QQ.check <- rbind(QQ.check,
                              rbind( data.frame(ref = proj.tis.match[i,1],
                                                comp = proj.tis.match[i,2],
                                                x.txt, mean.quant.1),
                                    data.frame(ref=proj.tis.match[i,2],
                                               comp=proj.tis.match[i,1],
                                               x.txt, mean.quant.2)))
        }
    }
    return(QQ.check)
}

QQ.plot <- function(QQ.check, file){
    
    QQ.check$label <- paste(QQ.check$ref, QQ.check$comp, sep=":")
    ind.tcga <- grep("TCGA", QQ.check$ref)
    ind.gtex <- grep("TCGA", QQ.check$ref, invert=TRUE)
    
    pdf(file = file, width=12, height=6)
    
    ggplot(QQ.check[ind.tcga,], aes(x=quant, y=obs.quant, color=label, shape=x.txt)) +
        geom_point() + geom_line() + xlab("quantile searched") +
            ylab("observed mean quantile")

    ggplot(QQ.check[ind.gtex,], aes(x=quant, y=obs.quant, color=label, shape=x.txt)) +
        geom_point() + geom_line() + xlab("quantile searched") +
            ylab("observed mean quantile")
    dev.off()

    print(paste0("Create plot ",file))
}


get.quant <- function(c, ref.col, comp.col, cutoffs, quants){

    obs.quants <- c()
    for (i in 1:length(quants)){

        quant = quants[i]
        cutoff <- cutoffs[i]
        ind <- which(c[,ref.col] > cutoff)
        obs.quant <- mean(rank(c[,comp.col])[ind]) / nrow(c)
        max.quant <- mean(nrow(c) - 0:(length(ind)-1))/nrow(c)
        exp.quant <- mean(nrow(c):(nrow(c)-length(ind))/2)/nrow(c)
        obs.quants <- rbind(obs.quants, c(quant, obs.quant, exp.quant, max.quant, cutoff))

    }
    obs.quants <- as.data.frame(obs.quants)
    names(obs.quants) <- c("quant", "obs.quant","exp.quant","min.quant","cutoff")
    return(obs.quants)
}

compile.surv.results <- function(projects, RsltDir){

    sexes <- c("all","male","female")
    ps <- coef <- list()
    for (project in projects){

        print(paste("Working on project", project))
        tmp <- get.surv.rslts(project, RsltDir)
        
        if (length(ps) == 0){
             for (sex in sexes){
                 
                 coef[[sex]] <- tmp[[1]][[sex]][,1:6]
                 names(coef[[sex]])[ncol(coef[[sex]])] <- project
                 ps[[sex]] <- tmp[[2]][[sex]][,1:6]
                 names(ps[[sex]])[ncol(ps[[sex]])] <- project
             }
         }else{
             
             for (sex in sexes){
                 if (length(tmp[[1]][[sex]]) == 1){ ## happens with na only
                     
                     ps[[sex]] <- cbind(ps[[sex]], NA)
                     coef[[sex]] <- cbind(coef[[sex]], NA)
                     names(coef[[sex]])[ncol(coef[[sex]])] <- project
                     names(ps[[sex]])[ncol(ps[[sex]])] <- project
                     next
                 }
                 
                 coef[[sex]] <- merge(coef[[sex]], tmp[[1]][[sex]][,c("ENSEMBL","gene.coef")],
                                      by = "ENSEMBL", all = T)
                 names(coef[[sex]])[ncol(coef[[sex]])] <- project
                 ps[[sex]] <- merge(ps[[sex]], tmp[[2]][[sex]][,c("ENSEMBL","gene.p")],
                                    by="ENSEMBL", all = T)
                 names(ps[[sex]])[ncol(ps[[sex]])] <- project
             }
         }
    }
    
    qs <- ps
    for (i in 1:length(qs)){
        
        qs[[i]][,projects] <- apply(qs[[i]][,projects], 2, function(x) make.q(x))
    }
    
    return(list(ps, qs, coef))
}

get.surv.rslts <- function(project, RsltDir){
    ps <- list()
    coef <- list()
    for (sex in c("all","male","female")){
        file <- paste0(RsltDir, project,"_lcpm.invnorm.covs_",sex,"_COEF.txt")
        if (!file.exists(file)){
            print(paste("No results for ", project, " : ",sex))
            ps[[sex]] <- NA
            coef[[sex]] <- NA            
        }else{

            coef[[sex]] <- read.table(file, as.is=T, header=T)

            file <- paste0(RsltDir, project,"_lcpm.invnorm.covs_",sex,"_PS.txt")
            ps[[sex]] <- read.table(file, as.is=T, header=T)
        }
    }

    return(list(coef, ps))
}

make.q <- function(ps){

    if (sum(is.na(ps)) == length(ps)){
        qs <- rep(NA, length(ps))
    }else{
        not.na.ind <- which(!is.na(ps))
        qs <- rep(NA, length(ps))
        qs[not.na.ind] <- qvalue(na.omit(ps))$qvalues
    }
    return(qs)
}

make.Qs <- function(logPs){
    Qs <- logPs
    for (i in 2:ncol(logPs)){
        p <- logPs[,i]
        miss.ind <- which(is.na(p))
        if (length(miss.ind) > 0){
            p <- p[-miss.ind]
        }
        small.ind<- which(p < 10^-30)
        p[small.ind] <- 10^-30
        q <- qvalue(p)$qvalues
        if (length(miss.ind) > 0){
            Qs[-miss.ind,i] <- q
        }else{
            Qs[,i] <- q
        }
    }
    return(Qs)
}

analyze.ranked.survival <- function(tissues, projects, ps, qs,
                                    logPs.gtex, logFC.gtex,
                                    logPs.tcga, logFC.tcga,
                                    quant=.001, niter=100, both=FALSE,
                                    PlotFile)
    {

        print(paste("Working on analysis using quant =",quant,
                    "and ",niter, "iterations to deterine 2.5% survival quantiles"))
        if (both == TRUE){
            print("Requiring gene to be significantly DE in both GTEx and TCGA")
        }

        surv.summary <- summarize.surv.sign(tissues, projects,
                                            ps, qs,
                                            logPs.gtex, logFC.gtex,
                                            logPs.tcga, logFC.tcga,
                                            quant=quant,
                                            niter = 500, both=both)
        mid <- surv.summary[[1]]
        lb <- surv.summary[[2]]
        ub <- surv.summary[[3]]      
        
        ## CREATE PLOT OF MEDIAN P VALUES AND EXPECTED MEDIAN P VALUES WITH +/- ##
        file <- gsub("\\.pdf", "_medp\\.pdf", PlotFile)
        plot.tcga.gtex.surv(mid, lb, ub, val="p.med", ylab="Median P-value",
                            file = file )
                    
        ## REPEAT AGAIN FOR MIN P VALUE ##
        file <- gsub("\\.pdf", "_minp\\.pdf", PlotFile)        
        plot.tcga.gtex.surv(mid, lb, ub, val="p.min",
                            ylab = "-log10(Minimum P-value)",
                            file= file )
        

        return(surv.summary)
    }
        
summarize.surv.sign <- function(tissues, projects,
                                ps, qs,
                                logPs.gtex, logFC.gtex,
                                logPs.tcga, logFC.tcga,
                                quant=.001,
                                niter=100, both){

    ## TISSUES are the studies from which the genes with DE in the quant
    ## quantile will be identified

    ## PROJECTS are the studies from which the survival statistics will be extracted
    ## from for the DE identified genes

    ## PS are the survival p-value summaries for the PROJECTS studies
    ## QS are the survival q-value summaries for the PROJECTS studies

    ## LOGPS.GTEX are the -log10(p-values) of DE from the TISSUE studies
    ## LOGFS.GTEX are the log(FC) of male versus female expression
    
    ## QUANT is the DE p-value quantile from logPs.gtex from which genes will be
    ## selected

    ## NITER is the number of iterations used to determine the 2.5% tail of the
    ## survival p-value distribution

    ## BOTH is true if the genes have to be DE in both gtex and tcga
    ##

    sexes <- c("all", "male", "female")
    out <- c()
    lb <- lb.alt <- c()
    ub <- ub.alt <- c()
           
    for (tissue in tissues){
        
        for (project in projects){
            
            P.t <- select(logPs.tcga, SYMBOL, eval(project)) %>% filter(!is.na(eval(project)))
            FC.t <- select(logFC.tcga, SYMBOL, eval(project)) %>% filter(!is.na(eval(project)))
            
            P.g <- select(logPs.gtex, SYMBOL, eval(tissue)) %>% filter(!is.na(eval(tissue)))
            FC.g <- select(logFC.gtex, SYMBOL, eval(tissue)) %>% filter(!is.na(eval(tissue)))            

            if (sum(tissue == project)){
                logPs <- P.t
                logFC <- FC.t
            }else{
                logPs <- inner_join(P.t, P.g, by="SYMBOL")
                logFC <- inner_join(FC.t, FC.g, by="SYMBOL")
            }

            
            thresh.val <- ceiling(nrow(logFC) * quant)

            ind <- order(logPs[,tissue], decreasing=T)

            male.ind <- ind[which(ind %in% which(logFC[,tissue] > 0))][1:thresh.val]
            female.ind <- ind[which(ind %in% which(logFC[,tissue] < 0))][1:thresh.val]
            
            genes.m <- logPs$SYMBOL[male.ind]
            genes.f <- logPs$SYMBOL[female.ind]
            genes <- unique(c(genes.m, genes.f))

            if (both == TRUE){
                
                ## in tcga and gtex 
                
                ind.t <- order(logPs[,project], decreasing=T)
                male.ind <- ind[which(ind %in% which(logFC[,project] > 0))][1:thresh.val]
                female.ind <- ind[which(ind %in% which(logFC[,project] < 0))][1:thresh.val]
                
                genes.m <- intersect(logPs$SYMBOL[male.ind], genes.m)
                genes.f <- intersect(logPs$SYMBOL[female.ind], genes.f)
                genes <- unique(c(genes.m, genes.f))              
            }         
            
            
            for (sex in sexes){  ## TCGA survival subset (all, male, female)

                print(paste(project, tissue, sex))
                ## gtex genes in quantile all, male and female upregulated seperately ##
                gene.ind <- which(ps[[sex]]$SYMBOL %in% genes)
                gene.ind.m <- which(ps[[sex]]$SYMBOL %in% genes.m)
                gene.ind.f <- which(ps[[sex]]$SYMBOL %in% genes.f)
                
                P <- filter(ps[[sex]], SYMBOL %in% logPs$SYMBOL) 
                Q <- filter(qs[[sex]], SYMBOL %in% logFC$SYMBOL)

                ## CALCULATE NULL VALUES (RANDOMLY CHOOSE THE SAME NUMBER OF GENES ##
                sim.all <- sim.surv.summary(P, Q, project, size = length(gene.ind), niter=niter)
                sim.m <-  sim.surv.summary(P, Q, project, size = length(gene.ind.m), niter=niter)
                sim.f <- sim.surv.summary(P, Q, project, size = length(gene.ind.f), niter=niter)
                
                out <- rbind(out,
                             c(project, tissue, sex, "either",
                               surv.summary(P, Q, project, gene.ind)),
                             c(project, tissue, sex, "null.either",
                               sim.all[[1]]),
                             c(project, tissue, sex, "m",
                               surv.summary(P, Q, project, gene.ind.m)),
                             c(project, tissue, sex, "null.m",
                               sim.m[[1]]),
                             c(project, tissue, sex, "f",
                               surv.summary(P, Q, project, gene.ind.f)),
                             c(project, tissue, sex, "null.f",
                               sim.f[[1]])
                             )

                lb <- rbind(lb, c(project, tissue, sex, "null.either", sim.all[[2]]),
                            c(project, tissue, sex, "null.m", sim.m[[2]]),
                            c(project, tissue, sex, "null.f", sim.f[[2]]))

                ub <- rbind(ub, c(project, tissue, sex, "null.either", sim.all[[3]]),
                            c(project, tissue, sex, "null.m", sim.m[[3]]),
                            c(project, tissue, sex, "null.f", sim.f[[3]]))

                lb.alt <- rbind(lb.alt, c(project, tissue, sex, "null.either", sim.all[[4]]),
                            c(project, tissue, sex, "null.m", sim.m[[4]]),
                            c(project, tissue, sex, "null.f", sim.f[[4]]))

                ub.alt <- rbind(ub.alt, c(project, tissue, sex, "null.either", sim.all[[5]]),
                            c(project, tissue, sex, "null.m", sim.m[[5]]),
                              c(project, tissue, sex, "null.f", sim.f[[5]]))
            }                        
        }
    }

    out <- as.data.frame(out)
    lb <- as.data.frame(lb)
    ub <- as.data.frame(ub)
    lb.alt <- as.data.frame(lb.alt)
    ub.alt <- as.data.frame(lb.alt)
    names(out) <- names(lb) <- names(ub) <- names(lb.alt) <- names(ub.alt) <-
        c("project", "gtex.tis", "tcga.set",
          "gtex.upreg.in","p.mean",
          "p.med","p.min", "q.mean", "q.med", "q.min")
    
    return(list(out, lb, ub, lb.alt, ub.alt))
}

surv.summary <- function(P, Q, project, gene.ind){
    
    if (length(gene.ind) == 0){
        out <- rep(NA, 6)
    }else{
        
        P <- P[gene.ind , project]
        if (sum(!is.na(P)) > 0){
            Q <- Q[gene.ind , project]
            mean.p <- mean(P, na.rm=T)
            med.p <- median(P, na.rm=T)
            min.p <- min(P, na.rm=T)
            
            mean.q <- mean(Q, na.rm=T)
            med.q <- median(Q, na.rm=T)
            min.q <- min(Q, na.rm=T)
            
            out <- c(mean.p, med.p, min.p, mean.q, med.q, min.q)      
        }else{

            out <- rep(NA, 6)
        }
    }            
    return(out)
}

sim.surv.summary <- function(P, Q, project, size, niter=100, ci=0.025){

    inds <- which(!is.na(P[,project]))
    mid <- lb <- ub <- c()
    if (length(inds) == 0){
        mid <- rep(NA, 6)
        lb <- lb.alt <- rep(NA, 6)
        ub <- ub.alt <- rep(NA, 6)

    }else{
    
        store <- c()
        f <- function(){ ind <- sample(inds, size=size, replace=F);
                         surv.summary(P, Q, project, ind)}
        store <- replicate(niter, f())           

        lb <- apply(store, 1, quantile, probs = ci) 
        ub <- apply(store, 1, quantile, probs = 1 - ci)
        mid <- rowMeans(store)

        sds <- sqrt(apply(store, 1, var))
        adj <- qnorm(1-ci)
        lb.alt <- mid - sds*adj
        ub.alt <- mid + sds*adj
    }
    return(list(mid, lb, ub, lb.alt, ub.alt))
}
                                           


plot.tcga.gtex.surv <- function(mid, lb, ub, val="p.med",
                                ylab = "Median P-value", file){

    d <- as.tibble(mid)
    d <- select(d, project, gtex.tis, tcga.set, gtex.upreg.in, eval(val))
    d <- rename(d, val = eval(val))
    d <- mutate(d, val = as.numeric(as.character(val)))
    d$null <- grepl("null", d$gtex.upreg.in)
    d$gtex.upreg.in <- gsub("null\\.","",d$gtex.upreg.in)

    lb <- as.tibble(lb)
    lb <- select(lb, project, gtex.tis, tcga.set, gtex.upreg.in, eval(val))
    lb <- rename(lb, lb = eval(val))
    lb <- mutate(lb, lb = as.numeric(as.character(lb)))
    lb$gtex.upreg.in <- gsub("null\\.","",lb$gtex.upreg.in)
    
    ub <- as.tibble(ub)
    ub <- select(ub, project, gtex.tis, tcga.set, gtex.upreg.in, eval(val))
    ub <- rename(ub, ub = eval(val))
    ub <- mutate(ub, ub = as.numeric(as.character(ub)))
    ub$gtex.upreg.in <- gsub("null\\.","",ub$gtex.upreg.in)

    if (length(grep("min", val))==1){
        d <- mutate(d, val = -log10(val))
        lb <- mutate(lb, lb = -log10(lb))
        ub <- mutate(ub, ub = -log10(ub))
    }
    
    d <- full_join(d, lb, by=c("project","gtex.tis","tcga.set","gtex.upreg.in"))
    d <- full_join(d, ub, by=c("project","gtex.tis","tcga.set","gtex.upreg.in"))
    d$gtex.tis <- gsub("\\-", "_", d$gtex.tis)
    ind <- which(d$null==FALSE)
    d$lb[ind] <- d$ub[ind] <- d$val[ind]

    tcga.labs <- c(all = "Both", male = "M", female= "F")
    gtex.labs <- c(Skin_Not_Sun_Exposed_Suprapubic = "SkNE",
                   Skin_Sun_Exposed_Lower_leg = "SkE",
                   Thyroid = "Thr",
                   Liver = "Liv",
                   Lung = "Lng",
                   Brain_Cortex = "BrCx",
                   Brain_Caudate_basal_ganglia = "BrCau",
                   Cells_EBV_transformed_lymphocytes = "LCL",
                   TCGA_SKCM = "SKCM",  TCGA_THCA = "THCA",
                   TCGA_LIHC = "LIHC", TCGA_LUAD = "LUAD",
                   TCGA_LAML = "LAML" )
    
                  
    
    p <- list()
    pdf(file = file, width=14, height=8)
    for (proj in projects){

        tmp <- d %>% filter(project == proj)
        p[[proj]] <- ggplot(tmp, aes(x=gtex.upreg.in, y=eval(val),color=null)) + 
          geom_pointrange(aes(ymin=lb,
                              ymax=ub), size = .5) +
          facet_grid(.~gtex.tis + tcga.set,
                   labeller =
                   labeller(tcga.set=tcga.labs,
                            gtex.tis = gtex.labs)) +
          ggtitle(proj) +
          ylab(ylab) + xlab("Sex of upregulated GTEx genes") +
              theme(axis.text.x=element_text(angle = 90, hjust = 0)) +

                   scale_color_discrete(name="",
                         labels=c("Obsv", "Null"))
                  

        print(p[[proj]])
    
    }
    dev.off()

    print(paste("Wrote plot to",file))
}
        

summary.de.surv.plot <- function(tissues, projects,
                                 ps, qs,
                                 logPs.gtex, logFC.gtex,
                                 logPs.tcga, logFC.tcga,
                                 geneInfo, PlotDir,
                                 niter = 500){
    
    no.x.ind.p <- which(ps[[1]]$chr %in% c(1:22))
    no.x.ind.gtex <- which(geneInfo.gtex$chr %in% c(1:22))
    no.x.ind.tcga <- which(geneInfo.tcga$chr %in% c(1:22))

    de.study <- "gtex"
    if (sum(tisues %in% projects)){
        de.study = "tcga"
    }
    
    for (quant in c(0.001, 0.01)){
        
        file <- paste0(PlotDir,"SurvPvals_",de.study,".tcga_",quant,".pdf")
        surv.summary[["X"]][[quant]][["nb"]] <-
            analyze.ranked.survival(tissues, projects,
                                    ps, qs,
                                    logPs.gtex, logFC.gtex,  
                                    logPs.tcga, logFC.tcga,
                                    quant=0.001,
                                    niter=niter, both=FALSE,
                                    file)
        
        ## no X
        file <- paste0(PlotDir,"SurvPvals_",de.study,".tcga_noX",quant,".pdf")
        surv.summary[["NoX"]][[quant]][["nb"]] <-
            analyze.ranked.survival(tissues, projects,
                                    ps, qs,
                                    logPs.gtex[no.x.ind.gtex,], logFC.gtex[no.x.ind.gtex,],  
                                    logPs.tcga[no.x.ind.tcga,], logFC.tcga[no.x.ind.tcga,],
                                    quant=0.001,
                                    niter=niter, both=FALSE,
                                    file)


        if (quant == .001 ){ next }

        
        file <- paste0(PlotDir,"SurvPvals_",de.study,".tcga_inboth_",quant,".pdf")
        surv.summary[["X"]][[quant]][["both"]] <-                     
            analyze.ranked.survival(tissues, projects,
                                    ps, qs,
                                    logPs.gtex, logFC.gtex,
                                    logPs.tcga, logFC.tcga,
                                    quant=quant,
                                    niter=niter, both=TRUE,
                                    file)
        
        ## NO X ##
        file <- paste0(PlotDir,"SurvPvals",de.study,".tcga_inboth",quant,".pdf")
        surv.summary[["NoX"]][[quant]][["both"]] <-
            analyze.ranked.survival(tissues, projects,
                                    ps, qs,
                                    logPs.gtex[no.x.ind.gtex,], logFC.gtex[no.x.ind.gtex,],
                                    logPs.tcga[no.x.ind.tcga,], logFC.tcga[no.x.ind.tcga,],
                                    quant=quant,
                                    niter=niter, both=TRUE,
                                    file)
    }

}
        

  
get.sample.sizes <- function(projects){

    data.dir <- "/gpfs/data/stranger-lab/askol/TCGA2/Data/Expression/Aracne/"
    Ns <- c()
    for (project in projects){
        
        file.expression.male <- paste0(data.dir, project,"_lcpm_auto_male.txt")
        file.expression.female <- paste0(data.dir,project,"_lcpm_auto_female.txt")    
        expression.files <- c(file.expression.male, file.expression.female)
        
        wcs <- c()
        for (i in 1:2){
            
            wc <- system(paste0("awk -F' ' '{print NF; exit}' ",
                                expression.files[i]), intern=T)
            wc <- as.numeric(wc)
            wcs[i] <- wc
        }

        Ns <- rbind(Ns, c(project, wcs))
    }

    Ns <- data.frame(project= Ns[,1], N.m = Ns[,2], N.f = Ns[,3])
    return(Ns)
}

get.gtex.sample.sizes <- function(tissues){
    
    OrigDataDir <- paste0("/gpfs/data/gtex-group/sex_biased_regulation_v8/",
                          "sexDE_v8_final/meri/data/")
    Ns <- c()
    for (tissue in tissues){
        
        
        count_file <- paste0(OrigDataDir, "Phenotypes/",tissue,
                             ".v8.gene_counts.txt.gz")
        covs_file <- paste0(OrigDataDir, "/Covariates/covs.basalcovs.",tissue,
                            ".txt")
        sex.info <- read.table(file = covs_file, as.is=T, header=F, quote="",
                               comment.char = "", sep="\t")
        colnames(sex.info) = c('SUBJID', 'SEX','SMTSISCH','SMRIN','AGE')
        samps <- read.table(file = count_file, nrow=1, header=F, as.is=T)
        samps <- gsub("\\.","-",samps)
        
        males <- samps[samps %in% sex.info$SUBJID[sex.info$SEX=="M"]]
        females <- samps[samps %in% sex.info$SUBJID[sex.info$SEX=="F"]]
        
        Ns <- rbind(Ns, c(tissue, length(males), length(females)))

    }
    
    Ns <- data.frame(project= Ns[,1], N.m = Ns[,2], N.f = Ns[,3]) 
    
    return(Ns)
}


plot.nostudies <- function(NoDE, file, trunc="", min.studies = 5){

    ## TRUNC IS USED TO CREATE A BAR PLOT OF NUMBER OF GENES 
    ## SHARING DIFFERENT NUMBERS OF STUDIES. TRUNC IS USED TO START THE
    ## PLOTING FROM A LARGER NUMBER OF STUDIES IN ORDER TO HAVE RESOLUTION
    ## TO SEE THE NUMBER OF GENES FOR VALUES THAT HAVE FEWER GENES

    ## MIN.STUDIES IS USED FOR ASSESSING THE NUMBER OF GENES WITH AT LEAST
    ## MIN.STUDIES STUDIES WITH THAT GENE DE.
    
    pdf(file = file, height=8, width=12)
    
    t <- as.data.frame(table(NoDE$n))   
    t <- data.frame(NoStudies=as.numeric(t[,1]), Count = as.numeric(t[,2]))

    p <- ggplot(data=t, aes(x = NoStudies, y=Count)) +
        geom_bar(stat="identity", fill="steelblue") +
            theme_minimal() +
                theme(text=element_text(size=14)) +
                    ylab("Number of Genes") + 
                        xlab("Number of Studies Sharing Gene")

    print(p)

    if (trunc != ""){

        ind <- t$NoStudies >= trunc
        t <- t[ind,]
        p <- ggplot(data=t, aes(x = NoStudies, y=Count)) +
            geom_bar(stat="identity", fill="steelblue") +
                theme_minimal() +
                    theme(text=element_text(size=14)) +
                        ylab("Number of Genes") + 
                            xlab("Number of Studies Sharing Gene")
        
        print(p)
    }
    

    ## PLOT DE BY CHROMSOME
    allN <- as.data.frame(table(NoDE$chr), stringsAsFactors=FALSE)
    names(allN) <- c("Chr", "CountAll")
    filtN <- as.data.frame(table(NoDE$chr[NoDE$n >= min.studies] ), stringsAsFactors=FALSE)
    names(filtN) <- c("Chr", "Count")
    N <- inner_join(allN, filtN, by="Chr")
    N$ExpCount <- N$CountAll/sum(N$CountAll)*sum(N$Count)
    N <- select(N, -CountAll)
    N <- melt(N, id = "Chr", variable.name = "CountType",value.name="Count")

    N$Chr <- factor(N$Chr, levels=c(1:22,"X","Y"))

    p <- ggplot(data=N, aes(x=Chr, y=Count, fill=factor(CountType))) +
                 geom_bar(stat="identity", position="dodge") +
                theme_minimal() +
                theme(text=element_text(size=14), legend.title=element_blank()) +
                ylab(paste0("Number of Genes DE in > ", min.studies-1)) + 
                    xlab("Chromosome") +
                        scale_fill_discrete(labels = c("Expect", "Observed"))

    print(p)
    dev.off()
    
    print(paste0("Wrote plot in ",file))
}

plot.de.on.x <- function(Qs, geneInfo, PlotDir, file.pre="tcga"){

    if (any(names(Qs)=="ensgene")){
        Qs <- rename(Qs ,gene = ensgene)        
        geneInfo <- rename(geneInfo, gene=ensgene)
        geneInfo <- select(geneInfo, gene,chr,start)
    }else {
        Qs <- rename(Qs,gene=SYMBOL)
        geneInfo <- rename(geneInfo, gene=SYMBOL)
        geneInfo <- select(geneInfo, gene,chr,start)
    }
    
    xinfo <- filter(geneInfo, chr=="X")
    Qs.x <- inner_join(xinfo, Qs, by="gene")
    Qs.x <- select(Qs.x, -gene, -chr)
    Qs.x <- melt(Qs.x, id = "start", variable.name = "study", value.name="q")
    Qs.x <- mutate(Qs.x, qcat = cut(q, breaks=c(0,.01,.05,.10,.25,.5,1)))
    Qs.x <- mutate(Qs.x, start = start/1000000)
    Qs.x <- mutate(Qs.x, study = gsub("TCGA-","", study))
    clrs <- colorspace::diverge_hsv(21)[21-c(0:5)]
    colfunc<-colorRampPalette(c("red","yellow"))
    clrs <- colfunc(6)
    file <- paste0(PlotDir,file.pre,"_DEonX.pdf")
    pdf(file = file, width=12, height=8)
    p <- ggplot(data=Qs.x, aes(x=start, y=study, color=qcat)) +
        geom_point(size=.5) +
            theme_minimal() +
                scale_color_manual(values=clrs, name="DE q-value range") +
                    xlab("Position on X (Mb)") + ylab("") +
                        theme(legend.position="bottom") 
    

    print(p)
    dev.off()
}


plot.de.corr <- function(Qs, names1, names2="", file) {

    xlims <- c(0,.7)
    lims <- c(0, .4)
    
    if (length(names2)==1){
        corr <- Qs %>% select(names1) %>%
            cor(use="pairwise.complete.obs", method="spearman")
        corr.x <- Qs %>% filter(chr == "X") %>% select(names1) %>% 
            cor(use="pairwise.complete.obs", method="spearman")
        dc <- dist(corr) %>% hclust(); dcx <- dist(corr.x) %>% hclust()
        
        nms1 <- nms2<-dc$labels[dc$order]
        nmsx1 <- nmsx2 <- dcx$labels[dcx$order]
        
    }else{
        corr <- Qs %>% select(c(names1, names2)) %>%
            cor(use="pairwise.complete.obs", method="spearman")
        dc <- dist(corr) %>% hclust()
        nms <- dc$labels[dc$order]
        nms1 <- nms[nms %in% names1]; nms2 <- nms[nms %in% names2]
        corr <- corr[rownames(corr) %in% names1,
                     colnames(corr) %in% names2]
        
        corr.x <- Qs %>% filter(chr == "X") %>% select(c(names1,names2)) %>% 
            cor(use="pairwise.complete.obs", method="spearman")       
        dcx <- dist(corr.x) %>% hclust()
        nmsx <- dcx$labels[dcx$order]
        nmsx1 <- nms[nms %in% names1] 
        nmsx2 <- nmsx[nmsx %in% names2]
        corr.x <- corr.x[rownames(corr.x) %in% names1,
                         colnames(corr.x) %in% names2 ]
        
    }
    
    corr <- corr %>% as.data.frame() %>% mutate(study=rownames(corr)) %>%
        melt(id = "study", variable.name = "study2",value.name="cor") %>%
            mutate(cor = ifelse(cor < lims[1], lims[1]+.001, cor)) %>%
                mutate(cor = ifelse(cor > lims[2], lims[2]-.001, cor))
    
    corr$study <- factor(corr$study, levels=nms1)
    corr$study2 <- factor(corr$study2, levels=nms2)

    
    corr.x <- corr.x %>% as.data.frame() %>% mutate(study=rownames(corr.x)) %>%
        melt(id = "study", variable.name = "study2",value.name="cor") %>%
            mutate(cor = ifelse(cor < xlims[1], xlims[1]+.001, cor)) %>%
                 mutate(cor = ifelse(cor > xlims[2], xlims[2]-.001, cor))
    corr.x$study <- factor(corr.x$study, levels=nmsx1)
    corr.x$study2 <- factor(corr.x$study2, levels=nmsx2)

    pdf(file=file)
    p <- ggplot(data=as.data.frame(corr), aes(x=study, y=study2, fill=cor)) +
        geom_tile() +
            scale_fill_gradient2("Spearman Correlation", high = "blue",
                                 low="white", limits=lims)
    p <- p + xlab("") + ylab("") + ggtitle("Correlation using all genes") +
        theme(text=element_text(size=10),
              axis.text.x=element_text(angle=60, hjust=1)) +
                  theme(legend.position="bottom") 

    print(p)

    p <- ggplot(data=as.data.frame(corr.x), aes(x=study, y=study2, fill=cor)) +
        geom_tile() +
            scale_fill_gradient2("Spearman Correlation", high = "blue",
                                 low="white", limits=xlims)
    p <- p + xlab("") + ylab("") + ggtitle("Correlation using genes on X") +
        theme(text=element_text(size=10),
              axis.text.x=element_text(angle=60, hjust=1)) +
                   theme(legend.position="bottom") 
    print(p)
    
    dev.off()

    print(paste0("Wrote plots to ",file))

    return(list(corr = corr, corr.x = corr.x))
}
               
plot.cond.de.distribution <- function(N, file){

    N.max.gtex <- max(N$N.gtex); N.max.tcga <- max(N$N.tcga)
    N <- mutate(N, n.gtex.imp = round(n.gtex*N.max.gtex/N.gtex)) %>%
        mutate(n.tcga.imp = round(n.tcga*N.max.tcga/N.tcga))
    N.use <- subset(N, n.gtex.imp > 0 | n.tcga.imp > 0)
    pdf(file = file, width = 20, height = 10)
    ## CONDITIONAL ON GTEX ##
    p <- ggplot(data=N, aes(x=factor(n.gtex.imp), y=n.tcga.imp,
                    fill=cut(..count.., c(0,1,5,10,20,50,100,500,1000,Inf)))) +
                        geom_bin2d() +  ylab("No TCGA studies") + xlab("No GTEx tissues") +  
                            scale_fill_hue("count")

    print(p)

  

    dev.off()
    print(paste0("Wrote plots to ",file))
}

    
