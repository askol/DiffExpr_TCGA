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


collect.results.gtex <- function(tissues, ResultDirGtex){
    
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


load.results <- function(project, ResultDir){

    file <- paste0(ResultDir, project, "_DE_limma.rslts")

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
        




  
