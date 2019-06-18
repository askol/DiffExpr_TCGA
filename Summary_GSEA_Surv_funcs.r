
collect.de.results <- function(projects, ResultDir){
    
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


collect.surv.data <- function(SurvResultDir){

    patrn <- "TCGA_DESurv_Best\\.txt$|GTEx_DESurv_Best\\.txt$"
    files <- dir(ResultDir, pattern=patrn)
    data <- c()
    
    for (file in files){

        if (file.exists(paste0(ResultDir,file)) == FALSE){ next }
        print(paste0("Working on file: ",file))
        data.tmp <- read.table(file = paste0(ResultDir,file), header=T, as.is=T)        
        data.int <- filter(data.tmp, emppmn<.05 | emppmd<.05 | emppmin<.05) %>%
            filter(Psurv < 0.001)

        data <- rbind(data, data.int)

    }

    return(data)
}


get.data <- function(project, DataDir){

    file <- paste0(DataDir, project, ".RDS")
    data <- readRDS(file)
    gene.info <- data$genes[,c("SYMBOL","ENSEMBL","chr","start","biotype")]
    
    case.ids <- data$samples$ID
    
    covs <- names(data$samples)[grep("SV",names(data$samples))]
    
    ## ADD BIRTHYEAR IF LESS THAN 10% IS MISSING VALUES ##
    if (grep("BirthYear", names(data$samples))){

        if (mean(is.na(data$samples$BirthYear)) < 0.10){
            covs <- c(covs, "BirthYear")
        }
    }

    meta <- get.surv.data(project, case.ids)
    meta$gender <- as.factor(meta$gender)
    
    ## THREE POSSIBLE DATA TYPES ##
    ## 1) lcpm,
    ## 2) inverseNormal expression,
    ## 3) inverseNormal expression after removing batch effects 
    lcpm <- cpm(data, log=TRUE)

    nuniqu <- apply(lcpm, 1, function(x) length(unique(x)))
    ## require > 50% of unique values ##
    na.ind <- which(nuniqu < ncol(lcpm)/2)
    data <- data[-na.ind,]
    data <- calcNormFactors(data, method="TMM")
        
    lcpRmBat <- removeBatchEffect(lcpm, covariates=data$samples[,covs])
    zRmBat <-  inverse_quantile_normalization(lcpRmBat)
    zRmBat <- as.data.frame(t(zRmBat))
    zRmBat$ID <- rownames(zRmBat)
    zRmBat <- merge(zRmBat, meta, by.x="ID", by.y="case_id", all=T)
    
    form <- paste("~", paste(covs, collapse="+"))
    design <- model.matrix(eval(parse(text=form)), data=data$sample)
    v <- voom(data, design, plot=F)   
    
    ## ANALYSIZE AFTER INVERSE NORMALIZATION ##
    cpm.invnorm <- inverse_quantile_normalization(v$E)
    cpm.invnorm <- as.data.frame(t(cpm.invnorm))
    cpm.invnorm$ID <- rownames(cpm.invnorm)
    cpm.invnorm <- merge(cpm.invnorm, data$samples[,c("ID",covs)],
                         by = "ID", all=T)
    cpm.invnorm <- merge(cpm.invnorm, meta, by.x="ID", by.y="case_id", all=T)
    
    logCPM <- cpm(data, log=TRUE, prior.count=3)
    logCPM <- as.data.frame(t(logCPM))
  
    logCPM$ID <- rownames(logCPM)
    logCPM <- merge(logCPM, data$samples[, c("ID", covs)],
                    by="ID", all=T)
    logCPM <- merge(logCPM, meta, by.x="ID", by.y="case_id", all=T)
        
    return(list(zRmBat, cpm.invnorm, logCPM, gene.info))
}

   
get.surv.data <- function(project, case.ids){

    META.DIR <- "/gpfs/data/stranger-lab/askol/TCGA/TCGA_MetaData/"
    
    keep.col <- c("hits.cases.case_id",
                  "hits.cases.demographic.gender",
                  "hits.cases.diagnoses.days_to_death",
                  "hits.cases.diagnoses.days_to_last_follow_up",
                  "hits.cases.diagnoses.age_at_diagnosis",
                  "hits.cases.demographic.race",
                  "hits.cases.demographic.ethnicity")
    meta <- read.table(file <- paste0(META.DIR, project,".METAdata.txt"),
                       header=T, as.is=T, sep="\t", quote="")
    meta <- meta[,keep.col]
    names(meta) <- gsub("hits\\..+\\.","",names(meta))
    
    meta$case_id <- paste0("A",meta$case_id)
    
    meta <- subset(meta, case_id %in% case.ids)
    
    ## REMOVE DUPES ##
    ind <- which(duplicated(meta$case_id))
    if (length(ind) >0 ){
        meta <- meta[-ind,]
    }

    ## REMOVE DATA POINTS WITH NO TIME OF EVENT ##
    ind <- which(is.na(meta$days_to_last_follow_up) &
                 is.na(meta$days_to_death))
    if (length(ind)){
        meta <- meta[-ind,]
    }
    ## FOR THE SUBSET OF SAMPLES THAT HAVE BOTH DAYS TO DEATH AND TO LAST FOLLOW-UP
    ## ALL DEATH DAYS > FOLLOW-UP, SO WILL USE DEATH DAYS ##
    meta$event.days <- apply(meta[,c("days_to_last_follow_up", "days_to_death")], 1,
                             max, na.rm=T)
    cens.ind <- which(is.na(meta$days_to_death)==T &
                      is.na(meta$days_to_last_follow_up)==F)
    meta$censored = 2
    meta$censored[cens.ind] = 1    
    return(meta)
}
  

calc.surv <- function(data, project, OutDir, prefix, gene.info, use.covs=FALSE){
    
    genes <- names(data)[grep("ENSG", names(data))]
    out <- c()
    data.sub <- list()
    data.sub[["all"]] <- data
    data.sub[["male"]] <- subset(data, gender=="male")
    data.sub[["female"]] <- subset(data, gender=="female")
    sexes = c("all", "male", "female")

    ## CREATE A FILE TO CONTAIN COEF, SE, P ONE FOR MALES, FEMALES AND COMBINED ##
    if (use.covs){
        cov.names <- names(data)[grep("SV|Birth", names(data))]
        covs <- get.best.covs(data, cov.names)
    }
    gene.ind <- which(names(data) %in% genes)
            
    for (sex in sexes){
        print(sex)
        d <- data.sub[[sex]]

        no.events <- sum(d$censored == 2, na.rm=T)
        if (no.events < 10){ next }
    
        s <- Surv(d$event.days, d$censored)
        
        if (use.covs){
            
            ## MAKE SURE THERE AREN'T TOO MANY COVARIATES BEING FIT ##
            ## MAKE SURE THAT THE NUMBER OF COVARIATES IS < 10% OF THE
            ## TOTAL NUMBER OF EVENTS
            max.covs <- ceiling(no.events * 0.10)
            covs.use <- covs
            if (length(covs) > max.covs){                
                covs.use <- covs[1:max.covs]             
            }
            if (sex == "all"){
                covs.use <- c(covs.use, "gender")
            }
            cov.ind <- which(names(d) %in% covs.use)
            
            eval(parse(text = paste0("srv <- lapply(gene.ind, function(x) ",
                           "coxph(s ~ d[,x] + ",
                           paste( paste0("d$", covs.use), collapse="+"),"))")))
        }else{
            ifelse(sex != "all",
                   srv <- lapply(gene.ind, function(x) coxph(s ~ d[,x])),
                   srv <- lapply(gene.ind, function(x) coxph(s ~ d[,x] + d[,"gender"]))  )
                                 
        }
        
        srv <- lapply(srv, summary)
        coef <- t(sapply(srv, function(x) x$coefficients[,1]))
        se <- t(sapply(srv, function(x) x$coefficients[,3]))
        ps <- t(sapply(srv, function(x) x$coefficients[,5]))
        

        if (nrow(coef) == 1){ coef <- t(coef); se <- t(se); ps <- t(ps)}
        coef <- as.data.frame(coef)        
        coef <- data.frame(gene=genes, coef)
        names(coef)[2] <- c("gene.coef")
        names(coef) <- gsub("d\\.","",names(coef))
        coef <- merge(gene.info, coef, by.x="ENSEMBL", by.y="gene", all.x=F,
                      all.y=T)
        
        se <-as.data.frame(se)
        se <- data.frame(gene = genes, se)
        names(se)[2] <- c("gene.se")
        names(se) <- gsub("d\\.","", names(se))
        se <- merge(gene.info, se, by.x="ENSEMBL", by.y="gene", all.x=F,
                    all.y=T)
        
        ps <- as.data.frame(ps)
        ps <- data.frame(gene=genes, ps)
        names(ps)[2] <- c("gene.p")
        names(ps) <- gsub("d\\.","",names(ps))
        ps <- merge(gene.info, ps,  by.x="ENSEMBL", by.y="gene", all.x=F,
                    all.y=T)
        
        ## write files out ##
        
        file <- paste0(OutDir, project,"_",prefix,"_",sex,"_COEF.txt")
        write.table(file = file, coef, quote=F, row.names=F, col.names=T)
        print(paste0("Wrote ",file))
        
        file <- paste0(OutDir, project,"_",prefix,"_",sex,"_SE.txt")
        write.table(file = file, se, quote=F, row.names=F, col.names=T)
        print(paste0("Wrote ",file))
        
        file <- paste0(OutDir, project,"_",prefix,"_",sex,"_PS.txt")
        write.table(file = file, ps, quote=F, row.names=F, col.names=T)
        print(paste0("Wrote ",file))

        warnings()
        print(warnings())
    }
    
}
        
get.best.covs  <- function(d, cov.names, sig.thresh=0.05){

    s <- Surv(d$event.days, d$censored)
    cov.ind <- which(names(d) %in% cov.names)
    eval(parse(text = paste0("srv <- lapply(cov.ind, function(x) ",
                           "coxph(s ~ d[,x]))")))

    srv <- lapply(srv, summary)
    ps <- t(sapply(srv, function(x) x$coefficients[,5]))
    min.ord <- order(ps, decreasing=FALSE)[1]
    sig.covs <- c()
    if (ps[min.ord] <= sig.thresh){
        sig.covs <- cov.names[min.ord]
        sig.ind <- which(names(d) == sig.covs)
    }

    if (length(sig.covs) == 0){ return(c()) }

    search=1
    while(search){
        
        cov.ind <- which(names(d) %in% cov.names)
        cov.ind <- setdiff(cov.ind, sig.ind)
        
        eval(parse(text = paste0("srv <- lapply(cov.ind, function(x) ",
                       "coxph(s ~ d[,x] + ",
                       paste(paste0("d$",sig.covs), collapse="+"),"))")))
        
        srv <- lapply(srv, summary)
        ps <- t(sapply(srv, function(x) x$coefficients[,5]))
        min.ord <- order(ps[,1], decreasing=FALSE)[1]
        
        if (ps[min.ord,1] <= sig.thresh){
            sig.covs <- c(sig.covs, names(d)[cov.ind[min.ord]])
            sig.ind <- which(names(d) %in% sig.covs)
        }else{
            search <- 0
        }
    }
    return(sig.covs)
}
 
