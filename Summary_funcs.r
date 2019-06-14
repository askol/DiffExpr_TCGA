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
## FOR PCA ##
library(FactoMineR)
library(factoextra)
library(missMDA)
library(cluster)
library(gplots)
library(ztable)
library(missRanger) ## for imputation for clustering ##
library(metap) ## for combining p-values

## 
## qqplot code ##
source("/gpfs/data/stranger-lab/askol/Code/qq_unif_plot.r")

collect.results <- function(projects, ResultDir){
    
    logFC <- c()
    logPs <- c()
    logQs <- c()
    geneInfo <- c()
    for (project in projects){
        
        print(paste0("Working on project : ",project))
        rslt <- load.results(project, ResultDir)
        dupe.ind <- which(duplicated(rslt$SYMBOL))
        dupe.names <- unique(rslt$SYMBOL[dupe.ind])

        if (length(dupe.ind) > 0){
            rslt <- filter(rslt, SYMBOL %in% dupe.names == F)
        }
        print(paste0("Removing ",length(dupe.names), " duplicate genes"))
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



load.results <- function(project, ResultDir){

    file <- paste0(ResultDir, project, "_DE_voom.rslts")
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

        ## DE_limma was performed using voom ##
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
        dupe.ind <- duplicated(rslt$ensgene)
        dupe.genes <- unique(rslt$ensgene[dupe.ind])

        if (length(dupe.genes)>0){
            rslt <- filter(rslt, ensgene %in% dupe.genes)
        }
        
        logFCtmp <- logPstmp <- c()
        if (any(grepl("Est", names(rslt)))){
            logFCtmp <- rslt[,c("ensgene", "Est")]
            logFCtmp <- rename(logFCtmp, logFC = Est)
            logPstmp <- rslt[,c("ensgene", "P.Value")]
        }else{
            logFCtmp <- rslt[,c("ensgene", "logFC")]
            logPstmp <- rslt[,c("ensgene", "P.Value")]
        }
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
    trash <- plotNDE(Nsig.DE.tcga, file=file, max.count = 500)
    
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
    Ns$N.min <- apply(Ns[,c("N.m","N.f")], 1, min)
    mult.fact <- max(c(Ns$N.m, Ns$N.f))/max(Ns$NSig)
    
    ind <- which(names(Ns) == "tissue")
    if (length(ind)>0){ names(Ns)[ind] = "project" }

    ord <- order(Ns$N.min, decreasing=T)
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
    print(paste0("Wrote plots to ", file))
    return(Ns)
}

get.noDE.cancers.for.gene <- function(Qs, qthresh, geneInfo,
                                      Ntrunc, min.studies, 
                                      ResultDir, PlotDir, file.pre){

    NoDE <- Qs
    if (any(names(NoDE)=="SYMBOL")){
        NoDE <- rename(NoDE, gene=SYMBOL)
        geneInfo$gene <- geneInfo$symbol
    }else {
        NoDE <- rename(NoDE, gene=ensgene)
        geneInfo$gene <- geneInfo$ensgene
    }

    n <-  rowSums(NoDE[,-1] <= qthresh, na.rm=T)
    N <- rowSums(!is.na(NoDE[,-1]))
    NoDE$n <- n
    NoDE$N <- N
    NoDE.pre <- select(NoDE, gene, n, N)
    NoDE.pre <- inner_join(geneInfo, NoDE.pre, by="gene") %>% arrange(desc(n))
    NoDE.post <- select(NoDE, -n, -N)
    ind <- which(names(NoDE.post) == "gene")
    NoDE.post[,-ind] <- 1*(NoDE.post[,-ind] <= qthresh)
    NoDE <- inner_join(NoDE.pre, NoDE.post, by="gene")
    NoDE <- select(NoDE ,-gene, -entrez)
    ## distribution of number of studies with same de gene ##
    file = paste0(ResultDir,file.pre,"_nodestudiespgene.txt")
    write.table(file = file, NoDE, quote=F, row.names=F, col.names=T)

    file = paste0(PlotDir,file.pre,"_NoDEStud.per.gene.pdf")
    tmp <- plot.nostudies(NoDE, file, trunc=Ntrunc, min.studies=min.studies)        
    
    return(list(NoDE = NoDE, t=tmp$t, t.nox=tmp$t.nox))
}

get.common.common.DE.genes <- function(N.gtex, N.tcga, Qs.tcga, Qs.gtex,
                                       thresh.gtex = 10, thresh.tcga=5,
                                       ResultDir, PlotDir){

    gi <- grch38[,c("symbol","ensgene","chr","start")]
    N <- select(N.gtex, ensgene, n, N)
    names(N) <- gsub("n$","n.gtex", names(N))
    names(N) <- gsub("N$","N.gtex", names(N))
    N <- left_join(N, gi, by = "ensgene")
    N <- inner_join(N, N.tcga[,c("ensgene","n","N")], by="ensgene")
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
    N <- mutate(N, n.gtex.imp = round(n.gtex*n.max.gtex/N.gtex)) %>%
        mutate(n.tcga.imp = round(n.tcga*n.max.tcga/N.tcga))

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

get.n.dist <- function(Ns){
    
    n.max.gtex <- max(Ns$n.gtex); n.max.tcga <- max(Ns$n.tcga)

    N.dist.g.g.t <- matrix(0,n.max.gtex+1, n.max.tcga+1)
    rownames(N.dist.g.g.t) <-  0:n.max.gtex
    colnames(N.dist.g.g.t) <- 0:n.max.tcga
    N.dist.t.g.g <- t(N.dist.g.g.t)
    
    for (i in unique(Ns$n.gtex.imp)){
        ## condition value in the columns (gtex)
        tmp <- filter(Ns, n.gtex.imp==i) %>% select(n.tcga.imp) %>% table()
        N.dist.t.g.g[as.numeric(names(tmp))+1,i+1] <- tmp/sum(tmp)
    }
    ## <= cum prob P(n.tcga.de <= obs | n.gtex.de = j)
    N.dist.t.g.g.le <- apply(N.dist.t.g.g,2,cumsum)
    ## >= cum prob P(n.tcga.de >= obs | n.gtex.de = j)
    N.dist.t.g.g.ge <- apply(N.dist.t.g.g,2, function(x) rev(cumsum(rev(x))))      

    for (i in unique(Ns$n.tcga.imp)){
        ## conditional value in the column (tcga)
        tmp <- filter(Ns, n.tcga.imp==i) %>% select(n.gtex.imp) %>% table()
        N.dist.g.g.t[as.numeric(names(tmp))+1, i+1] <- tmp/sum(tmp)
    }

    ## >= cum prob P(n.gtex.de <= obs | n.tcga.de = j)
    N.dist.g.g.t.le <- as.data.frame(apply(N.dist.g.g.t,2,cumsum))
    N.dist.g.g.t.le <- mutate( N.dist.g.g.t.le,
                              n.gtex.imp = as.character(rownames(N.dist.g.g.t.le))) %>%
                                  melt(id = "n.gtex.imp", variable.name = "n.tcga.imp"
                                       ,value.name="prob.g.g.t.le") %>%
                              mutate(n.tcga.imp = as.character(n.tcga.imp))

    ## <=cum prob P(n.gtex.de <= obs | n.tcga.de = j)
    N.dist.g.g.t.ge <- as.data.frame(apply(N.dist.g.g.t,2, function(x) rev(cumsum(rev(x)))))
    N.dist.g.g.t.ge <- mutate( N.dist.g.g.t.ge,
                      n.gtex.imp = as.character(rownames(N.dist.g.g.t.ge))) %>%
                          melt(id = "n.gtex.imp", variable.name = "n.tcga.imp",value.name="prob.g.g.t.ge") %>%
                              mutate(n.tcga.imp = as.character(n.tcga.imp))
                                          
    N.dist.t.g.g.le <- as.data.frame(apply(N.dist.t.g.g,2,cumsum))
    N.dist.t.g.g.le <- mutate( N.dist.t.g.g.le, n.tcga.imp = as.character(rownames(N.dist.t.g.g.le))) %>%
        melt(id = "n.tcga.imp", variable.name = "n.gtex.imp",value.name="prob.t.g.g.le") %>%
                              mutate(n.gtex.imp = as.character(n.gtex.imp))
       
    N.dist.t.g.g.ge <- as.data.frame(apply(N.dist.t.g.g,2, function(x) rev(cumsum(rev(x)))))
    N.dist.t.g.g.ge <- mutate( N.dist.t.g.g.ge, n.tcga.imp = as.character(rownames(N.dist.t.g.g.ge))) %>%
        melt(id = "n.tcga.imp", variable.name = "n.gtex.imp",value.name="prob.t.g.g.ge") %>%
                              mutate(n.gtex.imp = as.character(n.gtex.imp))

    N.dist <- N.dist.g.g.t.le
    N.dist <- inner_join(N.dist, N.dist.g.g.t.ge, by=c("n.gtex.imp","n.tcga.imp"))
    N.dist <- inner_join(N.dist, N.dist.t.g.g.ge, by=c("n.gtex.imp","n.tcga.imp"))
    N.dist <- inner_join(N.dist, N.dist.t.g.g.le, by=c("n.gtex.imp","n.tcga.imp"))

    Ns <- mutate(Ns, n.gtex.imp = as.character(n.gtex.imp)) %>%
        mutate(n.tcga.imp = as.character(n.tcga.imp))
    N.dist <- left_join(Ns, N.dist, by=c("n.gtex.imp", "n.tcga.imp"))

    N.dist <- mutate(N.dist, lab = ifelse(prob.g.g.t.le < 0.001 | prob.g.g.t.ge < 0.001 |
                                 prob.t.g.g.ge < 0.001 | prob.t.g.g.le < 0.001, symbol, ""))
    N.dist <- mutate(N.dist, n.gtex.imp = as.numeric(n.gtex.imp), n.tcga.imp=as.numeric(n.tcga.imp))


    ## PROBABILITY THRESHOLD FOR INTERESTING DEVIATION BETWEEN n.gtex and n.tcga ##
    int.thresh = 0.001
    
    ## Fewer DE genes in TCGA THAN EXPECTED (BASED ON GTEX)
    N.dist <- mutate(N.dist, interest.gt = ifelse(prob.t.g.g.ge <= int.thresh, 1, 0))
    N.dist <- mutate(N.dist, interest.gt2 = ifelse(prob.g.g.t.le <= int.thresh, 1, 0))
    ## MORE DE GENE IN TCGA THAN EXPECTED (BASED ON GTEX)
    N.dist <- mutate(N.dist, interest.lt = ifelse(prob.t.g.g.le <= int.thresh, 1, 0))
    N.dist <- mutate(N.dist, interest.lt2 = ifelse(prob.g.g.t.ge <= int.thresh, 1, 0))
    return(N.dist)
      
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
    NoDE <- mutate(NoDE, chr = ifelse(chr %in% c(1:22,"X","Y")==FALSE, "Oth", chr))
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

    ## CHISQ TEST FOR INDEPENDENCE ##
    t <- chisq.test(cbind(N$Count, N$CountAll-N$Count), simulate.p.value=TRUE, B=1000000)
    print("Chisq test using all chromosomes")
    print(chisq.test(cbind(N$Count, N$CountAll-N$Count)))
    
    ind <- which(N$Chr %in% c("X","Y"))
    t.nox <- chisq.test(cbind(N$Count[-ind], N$CountAll[-ind]-N$Count[-ind]),
                        simulate.p.value=TRUE, B=1000000)
    print("Chisq test excluding X and Y")
    print(chisq.test(cbind(N$Count[-ind], N$CountAll[-ind]-N$Count[-ind])))
    
    N <- select(N, -CountAll)
    N <- melt(N, id = "Chr", variable.name = "CountType",value.name="Count")

    N$Chr <- factor(N$Chr, levels=c(1:22,"X","Y"))

  
    p <- ggplot(data=N, aes(x=Chr, y=Count, fill=factor(CountType))) +
                 geom_bar(stat="identity", position="dodge") +
                theme_minimal() +
                theme(text=element_text(size=14), legend.title=element_blank()) +
                ylab(paste0("Number of Genes DE in > ", min.studies-1)) + 
                    xlab("Chromosome") +
                        scale_fill_discrete(labels = c("Observed", "Expected"))

    print(p)
    dev.off()
    
    print(paste0("Wrote plot in ",file))
    return(list(t=t, t.nox=t.nox))
}

plot.de.on.x <- function(Ps, FC, geneInfo, PlotDir, file.pre="tcga"){

    if (any(names(Ps)=="ensgene")){
        Ps <- rename(Ps ,gene = ensgene)
        FC <- rename(FC, gene = ensgene)
         FC[,-1] = -1*FC[,-1]
        geneInfo <- rename(geneInfo, gene=ensgene)
        geneInfo <- select(geneInfo, gene,chr,start)
    }else {
        Ps <- rename(Ps,gene=SYMBOL)
        FC <- rename(FC, gene=SYMBOL)
       
        geneInfo <- rename(geneInfo, gene=SYMBOL)
        geneInfo <- select(geneInfo, gene,chr,start)
    }

    xinfo <- filter(geneInfo, chr=="X")
    Ps.x <- inner_join(xinfo, Ps, by="gene")
    FC.x <- inner_join(xinfo, FC, by="gene")
    Ps.x <- select(Ps.x, -gene, -chr)
    FC.x <- select(FC.x, -gene, -chr)
    Ps.x <- melt(Ps.x, id = "start", variable.name = "study", value.name="p")
    Ps.x <- mutate(Ps.x, p = 10^(-1*p))
    FC.x <- melt(FC.x, id = "start", variable.name = "study", value.name="FC")

    Ps.x <- inner_join(Ps.x, FC.x, by=c("start","study"))
    Ps.x <- mutate(Ps.x, study = gsub("TCGA-","", study))
    Ps.x <- filter(Ps.x, !is.na(p) & !is.na(FC))
    
    Ps.x <- group_by(Ps.x, study) %>% mutate(r = rank(p, ties.method="random"))
    Ps.x <- group_by(Ps.x, study) %>%
        mutate(rnk = cut(rank(p)/sum(!is.na(p)), breaks=c(0,.01,.05,.1,.25,1)))
    Ps.x <- mutate(Ps.x, up.down =  as.factor(1*(sign(FC) == 1)))
    
    Ps.x <- mutate(Ps.x, start = start/1000000)
    cols <- c("#e7298a", "#d95f02","#fdbf6f", "#a6cee3", "#d9d9d9")
    szs <- c(3, 1.5, .5, .1, .1)
    shps <- c(16, 18)
    file <- paste0(PlotDir,file.pre,"_DEonX.pdf")
    pdf(file = file, width=12, height=8)
    p <- ggplot(data=Ps.x, aes(x=start, y=study, color=rnk, size=rnk)) +
        geom_point(aes(shape=up.down)) +
            theme_minimal() +
                scale_color_manual(values=cols, name="Rank range")+
                   scale_size_manual(values=szs, name="") +
                   scale_shape_manual(values=shps) +
                   xlab("Position on X (Mb)") + ylab("") +
                   theme(legend.position="bottom") 
    
    print(p)
    dev.off()
    print(paste0("Wrote plot to ",file))
    return(Ps.x)
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
    N.use <- subset(N, n.gtex.imp > 0 | n.tcga.imp > 0)

    ind <- which(N$n.gtex.imp == 44)
    pdf(file = file, width = 20, height = 10)    
    ## CONDITIONAL ON GTEX ##
    p <- ggplot(data=N, aes(x=n.gtex.imp, y=n.tcga.imp,
                    fill=cut(..count.., c(0,1,5,10,20,50,100,500,1000,Inf)))) +
                        geom_bin2d() +  ylab("No TCGA studies") + xlab("No GTEx tissues") +  
                            scale_fill_hue("count")

    print(p)

    dev.off()
    print(paste0("Wrote plots to ",file))
}

removeMultiMapGenes <- function(geneInfo){

    if (any(names(geneInfo) %in% "SYMBOL")){
        geneInfo <- mutate(geneInfo, symbol = SYMBOL)
    }
    dupeGenes <- unique(geneInfo$symbol[duplicated(geneInfo$symbol)])

    for (gene in dupeGenes){

        t <- filter(geneInfo, symbol == gene) %>% select(chr) %>% table()
        r <- filter(geneInfo, symbol == gene) %>% select(start) %>% range(na.rm=T)
        r <- abs(diff(r))
        ## remove gene if on more than one chromosome or different start sites are
        ## > 500k apart ##
        if (length(t) > 1 | r > 500000){
            geneInfo <- filter(geneInfo, symbol != gene)
        }
    }
    return(geneInfo)
}    

plot.PCA <- function(logPs, file, geneInfo, exclX=F){

    ## LETS LOOK VIA PCA ##
    ## PCA and Corplot from here
    ## http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

    if (exclX == F){

        xgenes <- geneInfo$ensgene[geneInfo$chr %in% c("X","Y")]
        xgenes <- c(xgenes,
                    geneInfo$SYMBOL[geneInfo$chr %in% c("X","Y")])
        if (any(grepl("SYMBOL", names(logPs)))){
            logPs <- rename(logPs, ensgene = SYMBOL)
        }
        logPs <- logPs[logPs$ensgene %in% xgenes == F,]
        file <- gsub(".pdf","_exclX.pdf", file)
    }else{
        file <- gsub(".pdf", "_inclX.pdf", file)
    }
    train <- logPs[,-1]
    ## MAKE 30 THE MAXIMUM -LOG10 P-VALUE ##
    train <- sapply(train, function(x) ifelse (x>30, 30, x))
   
    ## REMOVE GENES MISSING IN MORE THAN 2 STUDIES ##
    ind <- which(rowSums(is.na(train)) > 2)
    genes <- logPs[-ind,1]
    train <- train[-ind,]
    train <- t(scale(train))

    ## impute missing values
    train <- imputePCA(train)
    
    train <- data.frame(train)
    train <- sapply(train, function(x) as.numeric(as.character(x)))
    train <- data.frame(tissue = names(logPs)[-1], train)
    train$tissue<-as.factor(train$tissue)
    
    ## for plotting
    colors = rainbow(length(unique(train$tissue)))
    names(colors) = unique(train$tissue)
      
    
    rownames(train) <- train$tissue
    res.pca <- PCA(train[,-1], scale.unit=FALSE, graph = FALSE)
    
    pdf(file)

    print(fviz_pca_ind(res.pca, col.ind = "cos2", axes=c(1,2),
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE, # Avoid text overlapping (slow if many points),
                 labelsize = 2))
    
    print(fviz_pca_ind(res.pca, col.ind = "cos2", axes=c(2,3), 
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE ,
                 labelsize=2))
    
    print(fviz_pca_ind(res.pca, col.ind = "cos2", axes=c(3,4), 
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE ,
                 labelsize=2))
    
    print(fviz_pca_ind(res.pca, col.ind = "cos2", axes=c(4,5),
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE ,
                 labelsize=2))
    
    dev.off()
    print(paste0("Wrote plots to file ",file))
}


get.GSEA.results <- function(projects, GSEADir){

    gsea <- c()
    for (project in projects){

        file <- paste0(GSEADir,gsub("-","_",project),"/",project,
                       "_msigdb.v6.2.symbols.gmt_GSEA_summary.txt")
        if (file.exists(file)){
            if (file.size(file) > 0){
                g <- read.table(file = file, as.is=T, header=T, sep="\t")
                g <- mutate(g, nox = 0)
                file <- paste0(GSEADir,gsub("-","_",project),"/",project,
                               "_msigdb.v6.2.symbols.gmt_nox_GSEA_summary.txt")
                if (file.size(file) > 0){
                    gnox <- read.table(file = file, as.is=T, header=T, sep="\t")            
                    gnox <- mutate(gnox, nox = 1)
                
                    g <- rbind(g, gnox)
                }
                g <- mutate(g, proj = project)
            }else{
                print(paste0("No GSEA results found for ",project))
                next
            }
        }
        gsea <- rbind(gsea, g)
    }

    return(gsea)
}
            

get.fgsea.results <- function(projects, GSEADir){

    gsea <- c()
    for (project in projects){

        file <- paste0(GSEADir,gsub("-","_",project),"/",project,
                       "_fGSEA_summary.txt")
        if (file.exists(file)){
            if (file.size(file) > 0){
                g <- read.table(file = file, as.is=T, header=T)
                g <- mutate(g, nox = 0)
                g <- select(g, -ES, -log2err)
                file <- paste0(GSEADir,gsub("-","_",project),"/",project,
                               "_fGSEA_noX_summary.txt")
                if (file.size(file) > 0){
                    gnox <- read.table(file = file, as.is=T, header=T)            
                    gnox <- mutate(gnox, nox = 1)
                    gnox <- select(gnox, -ES, -log2err, LEgenes)
                    g <- rbind(g, gnox)
                }
                g <- mutate(g, proj = project)
            }else{
                print(paste0("No GSEA results found for ",project))
                next
            }
        }
        gsea <- rbind(gsea, g)
    }

    return(gsea)
}


get.median.p.value <- function(g, file){
    
    tissues <- unique(g$proj)
    quants <- c()
    pdf(file=file)
    for (tissue in tissues){
        print(paste0("Plotting qq plots for ",tissue))
        for (NoX in c(0,1)){
            tit <- paste0(tissue, " : NoX = ",NoX)
            ps <- g %>% filter(nox == NoX, proj == tissue) %>% select(pval) %>% unlist()
            p <- qqunif.plot(ps, main=tit)
            print(p)
            quant <- quantile(ps)
            quants <- rbind(quants, c(tissue, NoX, quant))
        }
    }
    dev.off()
    print(paste0("Wrote plot ", file))
    return(quants)
}


        
        
        
summarizie.gsea.by.cancer <- function(g, qthr = 0.05) {

    ## examine the distribution of the tag and size of the enriched gene sets
    ## tag = # le genes / N gene set
    ## DETERMINE IF THE IS A PATTERN TO THE DISTRIBUTION THAT SHOULD BE TAKEN INTO CONSIDERATION
    
    g <- mutate(g, signif = (padj<=qthr)*1)

    ## look at the distribution of size and Ncontrib in significatn
    ## versus nonsignigicant

    projects <- unique(g$proj)
    
    for (NoX in c(0,1)){

        tmp <- g %>% filter(nox == NoX)
        
        file <- paste0(PlotDir,"GSEA_summary_NoX",NoX,".pdf")
        pdf(file = file)
        
        p<-ggplot(tmp, aes(x=size))+
            geom_histogram(color="black", fill="white")+
                facet_grid(signif ~ .) + ggtitle("Gene set size")
        
        q <- ggplot(tmp, aes(x=Ncontrib))+
            geom_histogram(color="black", fill="white")+
                facet_grid(signif ~ .) +
                    ggtitle("Number of contributing genes")
        print(p)
        print(q)
    
        for (project in projects){

            tmp <- g %>% filter(proj == project, signif==1, nox==NoX)
            p<-ggplot(tmp, aes(x=size))+
                geom_histogram(color="black", fill="white")+
                    facet_grid(signif ~ .) +
                        ggtitle(paste0(project, " : Gene set size"))
            
            q <- ggplot(tmp, aes(x=Ncontrib))+
                geom_histogram(color="black", fill="white")+
                    facet_grid(signif ~ .) +
                        ggtitle(paste0(project, " : Number of contributing genes"))
            print(p)
            print(q)
        }
        
        dev.off()
    }
    ## REPORT THE TOP TEN GENE SETS FOR EACH CANCER
    

    ## EXAMINE THE OVERLAP OF THE GENES IN THE SIGNIFICANT GENE SETS

}

gsea.OL <- function(g, NoX = 0, topN = 250, plot.file, LEcutoff = 0){

    projects <- unique(g$proj)
    gseaOL <- c()
    genesContrib <- c()
    
    pdf(file = plot.file, width=16, height = 10)
    for (project in projects){

        print(paste0("Working on project ",project))
        
        g.proj <- filter(g, proj == project, nox == NoX, Ncontrib > LEcutoff)
        g.proj <- arrange(g.proj, pval)
        g.proj <- g.proj[1:topN,]
        
        ## REMOVE GENESETS WITH SMALL NUMBERS OF LEADING EDGE GENES ##
        
        all.genes <- g.proj %>% select(LEgenes) %>% unlist()
        all.gene.sets <- sapply(all.genes, strsplit, split=";")
        all.genes <- do.call(c, all.gene.sets)
        Ngenes.wdupes <- length(all.genes)
        t <- table(all.genes)
        qgenes <- quantile(t)
        top <- sort(t, decreasing=T)[1:topN]
        all.genes <- unique(all.genes)
        Nunique.genes <- length(all.genes)

        genesContrib <- rbind(genesContrib,
                              c(project, Ngenes.wdupes, Nunique.genes, top))
        
        genes.common <- lapply(all.gene.sets, function(x) 1*(all.genes %in% x))

        genes.common <- do.call(rbind, genes.common)
        rownames(genes.common) <- g.proj$pathway
        colnames(genes.common) <- all.genes
        
        genes.common.cnt <- genes.common %*% t(genes.common)
        ns <- diag(genes.common.cnt)
        Ns <- diag(1/ns)
        U <- Ns %*% genes.common.cnt
        diag(U) <- 0
        ord1 <- apply(U, 1, function(x) order(x, decreasing=T)[1])
        ord2 <- apply(U, 2, function(x) order(x, decreasing=T)[1])
        ol1 <- apply(U, 1, function(x) sort(x, decreasing=T)[1])
        ol2 <- apply(U, 2, function(x) sort(x, decreasing=T)[1])

        n1 <- ns*ol1; n2 = ns[ord2]*ol2
        N1 <- ns; N2 = ns[ord2]

        best.gs.n1 <- colnames(U)[ord1]
        best.gs.n2 <- colnames(U)[ord2]
        out <- data.frame(project = project, set = colnames(U),
                          best.ol.n1 = ol1, n1 = n1, N1 = N1,
                          best.ol.n2 = ol2, n2 = n2, N2 = N2,
                          best.gs.n1 = best.gs.n1, best.gs.n2 = best.gs.n2)
        gseaOL <- rbind(gseaOL, out)

        ## CREATE PLOT OF GENE SET OVERLAP BY: X-AXIS GENES, Y-AXIS GENESETS ##
        ## FIGURE OUT HOW TO SORT GENES AND GENE SETS SO THAT PLOT IS PRETTY
        d <- melt(genes.common, varnames=c("set", "gene"), value.name = "isin")

        dist.bin <- function(x) dist(x, method="binary")
        hclust.ave <- function(x) hclust(x, method="average")
        cols <- c("gray95","black")
        heatmap.2(genes.common, trace="none", distfun=dist.bin,
                  hclustfun=hclust.ave,scale="row", col = cols,
                  key=FALSE, cexRow=.3, cexCol=.3, main=project,
                  xlab = "Gene")
                
    }
    
    dev.off()
    print(paste0("Created plot in ",plot.file))
    genesContrib <- data.frame(project = genesContrib[,1],
                               apply(genesContrib[,-1], 2, as.numeric))
    return(list(gseaOL = gseaOL, genesContrib = genesContrib))    
}

summarize.OL <- function(gc){

    ## MEAN NUMBER OF CONTRIBUTING GENES IN TOP 100 GENES SETS WITH AND WITHOUT DUPES 
    contr.genes.mean <- colMeans(gc[,2:3])
    print("Mean unique genes in top 100 gene sets / Total number of genes in top 100 GS")
    print(    contr.genes.mean[2]/contr.genes.mean[1])
    top100contrib <- rowSums(gc[,-c(1:3)])
    top100contrib <- cbind(top100contrib, top100contrib/gc[,2])
    top100contrib <- cbind(top100contrib, top100contrib[,1]/100)
    top100contrib <- as.data.frame(top100contrib)
    names(top100contrib) <- c("NgsTop100", "PropGenesRepTop100", "MeanGSgeneIn")
    ## mean number of genesets the 100 mostly hightly represented genes are in ##
    print(paste0("Proportion of genes represented by 100 most common genes in genesets / ",
          "Total number of genes across all gene sets (including dupes)"))
    print(mean(top100contrib$PropGenesRepTop100))
    print("Mean number of gene sets 100 most common genes are found in")
    print(mean(top100contrib$MeanGSgeneIn))
    return(top100contrib)        
}

reportTop10 <- function(g, NoX=0, TopRank = 10, maxlogP = 39.99, file){

    
    g <- filter(g, nox == NoX)
    g <- g %>% group_by(proj) %>% mutate(rnk = rank(pval))
    g.nofilt <- g
    g <- filter(g, rnk <= TopRank)
    g <- select(g, -rnk, -nox)
    g <- select(g, proj, everything())
    
    write.table(file = file, g, quote=F, col.names = T, row.names=F)
    print(paste0("Wrote table to ",file))
    
    paths <- unique(g$pathway)
    h <- filter(g.nofilt, pathway %in% paths)
    h <- select(h, pathway, proj, padj)
    h <- mutate(h, padj = -log10(padj))
    h <- mutate(h, padj = ifelse(padj >= maxlogP, maxlogP, padj))
    h <- dcast(h, pathway ~ proj)
    h <- mutate(h, pathway = substr(pathway, start=1, stop=60))
  
    path.ord <- names(sort(table(g$pathway), decreasing=T))
    ord <- match(path.ord, h$pathway)
    h <- h[ord,]

    tbl.file <- gsub("\\.txt", "_pretty.html", file)
    options(ztable.type="html")
    z <- ztable(h) %>% makeHeatmap(palette="Blues")
    sink(tbl.file)
    print(z)
    sink()
    
    print(paste0("Wrote ", tbl.file))

    return(h)
}
    
make.cytoscape.out <- function(g, outDir){

    projects <- unique(p$proj)

    
    for (project in projects){

        for (NoX in 0:1){

            file <- paste0(outDir, "Cyto_", project, "_nox_",NoX, ".txt")
            GMTFile <- paste0(outDir, project,"_",NoX,"_LEGeneSets.gmt")

            ## MAKE GMT FILE
            gs <- filter(g, proj == project, nox == NoX)
            gs <- select(gs, pathway, LEgenes) %>% mutate(pathway = paste(pathway, pathway, LEgenes,sep=";"))
            write.table(file = GMTFile, gs$pathway, quote=F, row.names=F, col.names = F)
            cmd <- paste("sed -i  \'s/;/\\t/g\' ",GMTFile)
            system(cmd)                               

            ## MAKE CYTOSCAPE INPUT FILE
            out <- filter(g, proj == project, nox == NoX) %>% select(pathway, pval, padj, NES) %>%
                mutate(Phenotype = ifelse(sign(NES) > 0, "+1", "-1"))
            
            out <- mutate(out, ID = pathway) %>% rename(Description = pathway, p.Val = pval,
                                   FDR = padj)
            out <- select(out, ID, Description, p.Val, FDR, Phenotype)
            write.table(file = file, out, row.names=F, col.names=F, quote=F, sep="\t")
        }
    }
}
            
genesetOL <- function(g, NoLE = 5, qThresh = 0.05,
                      maxGS = 1000, ResultDir){

    for (NoX in c(0,1)){

        gs <- filter(g, nox == NoX,  Ncontrib >= NoLE)
        gs <- group_by(gs, proj) %>% mutate(rnk = rank(pval, ties.method="random")) %>% ungroup()
        gs <- filter(gs, padj <= qThresh & rnk <= maxGS)
        gs <- gs %>% group_by(pathway) %>%
            mutate(Npath = length(proj)) %>% ungroup()
        gs <- arrange(gs, desc(Npath), pathway, desc(Ncontrib))
        
        t <- select(gs, pathway) %>% table()

        file <- paste0(ResultDir, "GenesetSharing_NoX_",NoX,
                       ".txt")
        write.table(file = file, gs, quote=F, row.names=F,
                    col.names=T)
        print(paste0("Wrote table to ",file))
        
        file <- paste0(ResultDir,"Plots/ShareGeneSets_NoX",NoX,".pdf")
        pdf(file = file)
        ## PLOT GENESETS SHARED BY 10 (6) OR MORE CANCERS
        for (Np in c(10, 6)){

            cytofilePre <- paste0(ResultDir,"Tables/CytoSharedBy",Np,"_NoX",NoX)
            gs %>% filter(Npath >= Np) %>% makeCytoFiles(filePre = cytofilePre)
            gsp <- filter(gs, Npath >= Np) %>%
                select(pathway, proj)
            gsp$here <- 1
            gsp <- dcast(gsp, pathway ~ proj, value.var="here" )
            gsp[is.na(gsp)] <- 0
            rownames(gsp) <- gsp$pathway
            colnames(gsp) <- gsub("TCGA-","",colnames(gsp))
            gsp <- select(gsp, -pathway)
            gsp <- as.matrix(gsp)
            dist.bin <- function(x) dist(x, method="binary")
            hclust.ave <- function(x) hclust(x, method="average")
            cols <- c("gray95","black")
            mn <- paste0("Gene sets in > ",Np-1, " projects, ExclX=",NoX)
            labs <- heatmap.2(gsp, scale="none", trace="none",
                              distfun=dist.bin,
                              hclustfun=hclust.ave, col = cols,
                              key=FALSE, cexRow=.3, cexCol=.6, main=mn,
                              xlab = "Project")
            file.tbl <- paste0(ResultDir,"Tables/ShareGeneSets_NoX",NoX,
                               "_NoSharedSets_",Np,".txt")
            write.table(file = file.tbl,
                        gsp[rev(labs$rowInd), labs$colInd],
                        quote=F, row.names=T, col.names=T)
            
        }
        
        dev.off()
        print(paste0("Wrote plots to ",file))
                
    }
}
    
makeCytoFiles <- function(gs, filePre){

    ## 
    file <- paste0(filePre, ".txt")
    GMTFile <- paste0(filePre, ".gmt")

    ## MAKE GMT FILE
    ## ## EACH GENEST MAY HAVE MULTIPLE GENE SETS. CONCATINATE ALL LE GENES
    ## ## INTO SINGLE GENE SET
    gs.proj <- select(gs, pathway, proj) %>% mutate(val = 1)
    gs.proj <- dcast(gs.proj, pathway~proj, value.var = "val")
    gs.proj[is.na(gs.proj)] <- 0
    gs.proj$projCode <- apply(gs.proj[,-1], 1, paste, collapse="")
    gs.proj <- select(gs.proj, pathway, projCode)
    
    gs.concat <- group_by(gs, pathway) %>% mutate(Genes = paste(LEgenes, collapse=";"))
    gs.concat <- group_by(gs.concat, pathway) %>%
        mutate(Phenotype = mean(sign(NES))) %>%
            mutate(Phenotype = ifelse(Phenotype > 0, "+1","-1")) %>%
                mutate(N = length(NES)) %>% ungroup()
        
    
    gs.concat <- gs.concat[-which(duplicated(gs.concat$pathway)),]
    
    out <- mutate(gs.concat, out = paste(pathway, pathway, Genes,sep=";"))
    write.table(file = GMTFile, out$out,
                quote=F, row.names=F, col.names = F)
    cmd <- paste("sed -i  \'s/;/\\t/g\' ",GMTFile)
    system(cmd)                               

    print(paste0("Wrote ",GMTFile))
    ## MAKE CYTOSCAPE INPUT FILE    
    out <- gs.concat %>% select(pathway, N, Phenotype) %>%
        mutate(ID = pathway, Description = pathway, p.Val = 1/N, FDR = N,
               )
    out <- inner_join(out, gs.proj, by="pathway")
    
    out <- select(out, ID, Description, p.Val, FDR, Phenotype, projCode)
    write.table(file = file, out, row.names=F, col.names=F, quote=F, sep="\t")
    print(paste("Wrote ",file))
}


clust.genes <- function(gsclust){

    clusters <- sort(unique(gsclust$Cluster))
    clust.genes <- c()
    clust.gs <- c()
    for (cluster in clusters){
        if (cluster == 1){ next }
        ind <- gsclust$Cluster == cluster
        genes <- gsclust$EnrichmentMap..Genes[ind]
        gss <- gsclust$EnrichmentMap..Formatted_name[ind]
        gss <- gsub("\n","", gss)
        common.genes <- get.common.genes(genes)
        clust.genes <- rbind(clust.genes, cbind(cluster, common.genes))
        clust.gs <- rbind(clust.gs, cbind(cluster, gss))
    }    
    return(list(cluster.genes = clust.genes, cluster.genesets = clust.gs))
}


get.common.genes <- function(genes, reportMin = 2){

    ## EXPECTING A VECTOR WITH GENES SEPERATED BY | ##
    N <- length(genes)
    genes <- paste(genes, collapse="|")
    genes <- strsplit(genes, split="\\|")[[1]]
    genes <- table(genes)
    genes <- genes[genes >= reportMin]
    genes <- as.data.frame(genes)
    names(genes) <- c("gene","count")
    genes$NGeneSets <- N
    return(arrange(genes, desc(count)))
}

gene.dist.test <- function(gene.gs.dist, cl.genes){

    gene.dist <- cl.genes
    gene.dist <- filter(gene.dist, cluster != 2)
    gene.dist <- gene.dist %>% group_by(gene) %>% mutate(sum.count = sum(count)) %>% ungroup()
    gene.dist <- gene.dist[!duplicated(gene.dist$gene),] %>% select(gene, sum.count)
    
    gene.dist <- inner_join(gene.dist, gene.gs.dist[,c("gene","gs.count")])
    gene.dist <- gene.dist[rowMeans(gene.dist[,c("sum.count","gs.count")]) >= 5,]
    chtest <- chisq.test(gene.dist[,2:3])
    gene.dist$stres <- chtest$stdres[,1]
    gene.dist <- gene.dist %>% arrange(desc(stres))

    return(gene.dist)
}

gs.gene.concat <- function(gs){

    ## concatinate all genes in all genesets ##
    gs.concat <- group_by(gs, pathway) %>% mutate(Genes = paste(LEgenes, collapse=";")) %>% ungroup()
    
    gs.concat <- gs.concat[-which(duplicated(gs.concat$pathway)),]
    genes <- sapply(gs.concat$Genes, function(x){
        unique(strsplit(x, split=";")[[1]])})
    names(genes) <- c()
    genes <- do.call(c, genes)

    gene.freqs <- as.data.frame(table(genes))
    names(gene.freqs) <- c("gene","gs.count")
    gene.freqs$freq <- gene.freqs$gs.count / length(genes)
    gene.freqs <- arrange(gene.freqs, desc(gs.count))
    return(gene.freqs)
}


## CREATE A BINARY HEAT MAP FOR EACH HIGH-LEVEL CLUSTER WITH
## PROJECT ON THE X, AND GENE SET ON THE Y
plot.cluster.proj.geneset <- function(genesets, gs, file){


     cluster.labs <- rbind(c(0,"Mitochondria/Phosphorylation/TCA"),
                          c(2,"Transcription/Translation/RNA Metabolism"),
                          c(3,"Cancer Gene Neighborhood (MORF)"),
                          c(4,"DNA packaging"),
                          c(5,"Mitosis/Cell cycle/Cancer/Chromosome Segregation"),
                          c(6,"Immunity"),
                          c(7,"Cytokines/Inflammation"),
                          c(8,"Cell adhesion and Immunity"),
                          c(9,"Monocytes"),
                          c(10, "Immune Response"),
                          c(11, "Cancer Modules"),
                          c(12, "mRNA processing/splicing"),
                          c(13, "Mircroglia"),
                          c(15, "Extracellular Matrix"),
                          c(16, "Cancer/Immunity/Inflammation"),
                          c(17, "Rapamycin"),
                          c(18, "Cancer Modules2"))         
                            
    dist.bin <- function(x) dist(x, method="binary")
    hclust.ave <- function(x) hclust(x, method="average")
    cols <- c("gray95","black")
    
    ## gs = gsea
    ## genesets = genesets belonging to high-level functions ##
    genesets <- as.data.frame(genesets, stringsAsFactors=FALSE)
    clusters <- unique(genesets$cluster)
    gs <- filter(gs, nox == 1)
    gs <- gs %>% group_by(proj) %>% mutate(quant = rank(pval)/length(pval)) %>% ungroup()

    pdf(file = file)
    ENR <- c()
    for (clust in cluster.labs[,1]){

        cluster.name <- cluster.labs[cluster.labs[,1] == clust,2]
        gs.clust <- filter(genesets, cluster==clust) %>% select(gss)

        gs.clust <- gs %>% filter(pathway %in% gs.clust$gss) %>%
            mutate(plot.padj = ifelse(padj < 0.05, 1, 0)) %>%
                mutate(plot.quant = ifelse(quant < .01, 1, 0))

        
        
        mn <- paste0(cluster.name, " (q < 0.05)")
        gsp <- dcast(gs.clust, pathway ~ proj, value.var="plot.padj" )
        rownames(gsp) <- gsp$pathway
        colnames(gsp) <- gsub("TCGA-","",colnames(gsp))
        gsp <- select(gsp, -pathway)
        gsp <- as.matrix(gsp)
        labs <- heatmap.2(gsp, scale="none", trace="none",
                          distfun=dist.bin,
                          hclustfun=hclust.ave, col = cols,
                          key=FALSE, cexRow=.6, cexCol=.6,
                          xlab = "Project", dendrogram="none",
                          lhei=c(.1,1), lwid=c(.1,1),
                          margins=c(6,8), srtRow = -45)
         title(mn, cex.main=1)

        mn <- paste0(cluster.name, " (quant < 0.01)")
        gsp <- dcast(gs.clust, pathway ~ proj, value.var="plot.quant" )
        rownames(gsp) <- gsp$pathway
        colnames(gsp) <- gsub("TCGA-","",colnames(gsp))
        gsp <- select(gsp, -pathway)
        gsp <- as.matrix(gsp)
        heatmap.2(gsp, scale="none", trace="none",
                  distfun=dist.bin,
                  hclustfun=hclust.ave, col = cols,
                  key=FALSE, cexRow=.6, cexCol=.6, 
                  xlab = "Project", dendrogram="none",
                   lhei=c(.1,1), lwid=c(.1,1),
                  margins=c(6,8), srtRow = -45)
        title(mn, cex.main=1)

        ## DETERMINE WHICH CANCERS ARE THE MOST "ENRICHED" ##
        enrPadj <- gs.clust %>% group_by(proj) %>% summarize( fish.p.padj=sumlog(padj)$p) %>% 
            mutate(fish.p.padj = ifelse(fish.p.padj < 10^-100, 10^-100, fish.p.padj)) %>%
            mutate(enrPadj = (-log10(fish.p.padj) - min(-log10(fish.p.padj))) /
                   (max(-log10(fish.p.padj)) - min(-log10(fish.p.padj))))
        
        enrQuant <- gs.clust %>% group_by(proj) %>% summarize(fish.p.quant=sumlog(quant)$p) %>%
            mutate(fish.p.quant = ifelse(fish.p.quant < 10^-100, 10^-100, fish.p.quant)) %>%
            mutate(enrQuant = (-log10(fish.p.quant) - min(-log10(fish.p.quant))) /
                   (max(-log10(fish.p.quant)) - min(-log10(fish.p.quant))))
        enr <- inner_join(enrPadj, enrQuant, by = "proj") %>%
            mutate(cluster=cluster.name)
        ENR <- rbind(ENR, enr)
    }

     dev.off()
     print(paste0("Wrote plots to ",file))
     
     ENR <- mutate(ENR, proj = as.factor(gsub("TCGA-","",proj))) %>%
         mutate(cluster = as.factor(cluster))
     
     ENR <- mutate(ENR, enrPadj.size = as.factor(1 + 1*(enrPadj > .4) + 1*(enrPadj > .8))) %>%
         mutate(enrQuant.size = as.factor(1 + 1*(enrQuant > .4) + 1*(enrQuant > .8)))

     file <- gsub(".pdf", "_EnrichScore.pdf", file)
     pdf(file = file, width = 20, height = 10)

     cols <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2")[c(4,6:8)],
               brewer.pal(9,"Set1")[c(4,6,8,9)])
     sizes = c(2,4,7)

     for (enr in c("enrPadj", "enrQuant")){

         enr.size = paste0(enr,".size")
         ENR <- mutate(ENR, enr = !!sym(enr)) %>%  mutate(enr.size = !!sym(enr.size))
                           
         p <- ggplot(ENR, aes(x=proj, y=enr, color = cluster, size=enr.size)) +
         scale_color_manual(values=cols, name="High-level function") +
             xlab("Project") + ylab("Enrichment Score") +
                 scale_size_manual(values=sizes,
                                   labels = c("0 - 0.4","0.4 - 0.8", "0.8 - 1.0"),
                                   name = "Enrichment Score") +
                                       scale_y_continuous(breaks = seq(0,1,by=.2))
         p <- p + ggtitle(paste0("Enrichment Score (",gsub("enr","",enr),"based)")) + 
             geom_jitter(width=0, height = .05) +
                 theme(legend.text = element_text(size=12, face="bold"),
                       legend.title=element_text(size=14, face="bold"),
                       legend.background = element_rect(size=4),
                       legend.key = element_rect(fill = "white"),
                       axis.text=element_text(size=12),
                       axis.title=element_text(size=14,face="bold")) +
                           guides(color = guide_legend(override.aes = list(size=12)))
         print(p)
     }
     

     for (enr in c("enrPadj", "enrQuant")){

         a = ENR[,c("proj","cluster",enr)]        
         a <- rename(a, enr = !!sym(enr)) 
         a$cluster <- as.character(a$cluster)
         a$proj <- as.character(a$proj)
         
         E <- dcast(a, proj~cluster, mean, value.var="enr")
         rownames(E) <- E$proj
         E <- select(E, -proj) %>% as.matrix()
         E <- round(100*E)/100 %>% signif(digits=1)
         brks <- c(0,.4,.8,1)
         cols <- brewer.pal(3,"Blues")
         heatmap.2(E, scale="none", trace="none",Rowv = FALSE, Colv = FALSE,
                   key=FALSE, cexRow=1.2, 
                   xlab = "Project", dendrogram="none",               
                   lhei=c(.1,1), lwid=c(.1,1),
                   margins=c(6,8), srtCol = 15, cexCol=1.2,
                   col = cols, breaks = brks,
                   cellnote=E, notecex=1.0,
                   notecol="black", main = paste0("Enrichment Using ",gsub("enr","",enr)))
     }
     dev.off()

     
     print(paste0("Wrote plots to ",file))
     
}


plot.cluster.proj.genes.in.gs <-function(genes, logPs, file){

    make.palette <- function(x, max.arg=7){

        min.val <- min(x)
        max.val <- max(x)

        ll <- 0.07
        if (min.val < ll){ ll <- min.val-.001 }

        ul <- ceiling(min(max.arg, max.val))
        ## max ul x.3 instead of (x+1) if appropriate
        if (ul - max.val > 0.7){ ul = ul - 1 + 0.3}
        breaks <- c(ll, 0.7, seq(1.3, ul, by = 1), seq(2, ul, by=1), ul)
        breaks <- unique(breaks) ## gets rid of extra 10 when max is 10
        breaks <- sort(breaks)
        pal <- colorRampPalette(c("lightskyblue1", "purple4"))        
        cols <- pal(length(breaks)-1)
        return(list(breaks = breaks, cols = cols))
    }
    
    cluster.labs <- rbind(c(0,"Mitochondria/Phosphorylation/TCA"),
                          c(2,"Transcription/Translation/RNA Metabolism"),
                          c(3,"Cancer Gene Neighborhood (MORF)"),
                          c(4,"DNA packaging"),
                          c(5,"Mitosis/Cell cycle/Cancer/Chromosome Segregation"),
                          c(6,"Immunity"),
                          c(7,"Cytokines/Inflammation"),
                          c(8,"Cell adhesion and Immunity"),
                          c(9,"Monocytes"),
                          c(10, "Immune Response"),
                          c(11, "Cancer Modules"),
                          c(12, "mRNA processing/splicing"),
                          c(13, "Mircroglia"),
                          c(15, "Extracellular Matrix"),
                          c(16, "Cancer/Immunity/Inflammation"),
                          c(17, "Rapamycin"),
                          c(18, "Cancer Modules2"))                          
                            
##    dist.bin <- function(x) dist(x, method="binary")
##    hclust.ave <- function(x) hclust(x, method="average")
    clusters <- unique(genes$cluster)
    names(logPs) <- gsub("TCGA-","",names(logPs))
                         
    Q <- logPs
    rownames(Q) <- Q$SYMBOL; Q <- select(Q, -SYMBOL)
    Q <- apply(Q, 2, function(x){-log10(rank(-x, na.last="keep", ties.method = "random")/
                                     sum(!is.na(x)))})
    Q <- as.data.frame(Q)
    
    pdf(file = file)

    for (clust in clusters){
        
        cluster.name <- cluster.labs[cluster.labs[,1] == clust,2]

        g <- genes %>% filter(cluster == clust) %>% select(gene)
        g <- as.character(g$gene)

        ps <- logPs %>% filter(SYMBOL %in% g)
        qs <- Q[rownames(Q) %in% g,]
        
        rownames(ps) <- ps$SYMBOL; ps <- select(ps, -SYMBOL)
        tps <- as.data.frame(ps)
        qs <- as.data.frame(qs)

        file.pre <- gsub("\\.pdf","",file); file.pre= gsub("Plots/","",file.pre)
        file.out = paste0(file.pre,"_",gsub(" ","",
            gsub("/","_",cluster.labs[cluster.labs[,1] == clust,2])),"_logps.txt")
        write.table(file = file.out, ps, quote=F, row.names = TRUE, col.names=T)        
        print(paste0("Wrote :",file.out))
        
        file.out = gsub("logps","logqs",file.out)
        write.table(file = file.out, qs, quote=F, row.names = TRUE, col.names=T)
        print(paste0("Wrote :", file.out))

        imp <- missRanger(tps, pmm.k = 10, num.trees = 100)
        ps <- as.matrix(imp)

        imp <- missRanger(qs,  pmm.k = 10, num.trees = 100)
        qs <- as.matrix(imp)

        for (i in 1:2){
            tit = "P-value based"
            d = ps; if (i == 2){ d = qs; tit="Quantile based" }

            for (j in 1:2){
                if (i == 2 & j == 2){ next }
                if (i == 1 & j == 1){
                    ## MDS ##
                    ## SCALE ##
                    mds <- d %>% scale() %>% t() %>% dist() %>% cmdscale(k = 4) %>% as_tibble()
                    tit = "Scaled P-value based"
                }else{
                    ## DONT SCALE ##
                    if (i == 1){  tit = "P-value based" }
                    mds <- d %>% t() %>% dist() %>% cmdscale(k = 4) %>% as_tibble()
                }

                for (dim1 in 1:3){
                    for (dim2 in (dim1+1):4){
                        d1 = paste0("Dim.",dim1); d2 = paste0("Dim.",dim2)
                        colnames(mds)[c(dim1,dim2)] <- c(d1,d2)
                        ## Plot MDS
                        print(ggscatter(mds, x = d1, y = d2,
                                        label = colnames(d),
                                        size = 1,
                                        repel = TRUE,
                                        title=paste0(cluster.name," (",tit,")")) )
                    }
                }
            }
        }
        
        mn <- paste0(cluster.name)

        col.info <- make.palette(ps, max.arg=7)
        cols <- col.info$cols
        breaks <- col.info$breaks
        
        heatmap.2(ps, scale="none", trace="none",
                  ##                  distfun=dist.bin,
                  ##                  hclustfun=hclust.ave,
                  key=TRUE, cexRow=.25, cexCol=.8,
                  xlab = "Project", dendrogram="both",
                  breaks=breaks, col=cols, adjCol = c(1,.5),
                  key.xtickfun = function() {
                      breaks = pretty(parent.frame()$breaks)
                      breaks = breaks[c(1,length(breaks))]
                      list(at = parent.frame()$scale01(breaks),
                           labels = breaks)
                  })
        title(mn, cex.main=.8, adj=1)

        ## REPEAT BUT FOR QUANTILES ##
        col.info <- make.palette(qs)
        cols <- col.info$cols
        breaks <- col.info$breaks
        
        mn <- paste0(mn, "\n(Quantiles of DE p-values)")
        heatmap.2(qs, scale="none", trace="none",        
                  key=TRUE, cexRow=.25, cexCol=.8,
                  xlab = "Project", dendrogram="both",
                  breaks=breaks, col=cols, adjCol = c(1,.5),
                  key.xtickfun = function() {
                      breaks = pretty(parent.frame()$breaks)
                      breaks = breaks[c(1,length(breaks))]
                      list(at = parent.frame()$scale01(breaks),
                           labels = breaks)
                  })
        title(mn, cex.main=.8, adj=1)
    }
                     
    dev.off()
    print(paste0("Wrote plots to ",file))
}

    

get.most.ols <- function(d){

    gsets <- c()
    ngsets <- nrow(d)
    
    d$maxn <- apply(d[,c("n1","n2")], 1, max)
    d <- arrange(d, maxn)
}
    
gsea.OL2 <- function(g, qthresh=0.05){

    ## THE GOAL HERE IS TO DETERMINE WHICH GENE SETS ARE BEING FOUND IN MULTIPLE CANCERS

    ## perform for both autosomes only and all chromosomes
    for (NoX in c(0,1)){

        g.use <- filter(g, nox == NoX, padj<=qthresh)
    }
}
