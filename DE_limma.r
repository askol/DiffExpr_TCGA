DataDir <- "/gpfs/data/stranger-lab/askol/TCGA2/Data/Expression/"
ResultDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Results/"

library(limma)
library(edgeR)

run_limma <- function(project){

    data <- get.data(project, DataDir)
    
    count <- data[[3]]
    logCPM <- data[[2]]
    design <- data[[1]]
    design.nosvs <- data[[4]]

    ## LIMMA WITHOUT VOOM (INCLUDE SVS AND BIRTHYEAR) ##
    fit <- lmFit(logCPM, design)
    fit <- eBayes(fit, trend=TRUE)
    limma.rslt <- topTable(fit,coef=2, n=Inf, sort.by="none")

    limma.rslt <- data.frame(gene = count$genes[,1], count$genes[,-1], limma.rslt)
    coef <- as.data.frame(fit$coefficient)
    names(coef) <- paste("coef", gsub("\\(|\\)", "", names(coef)), sep=".")
    se <- as.data.frame(fit$s2.post*fit$stdev.unscaled)
    names(se) <- paste("se", gsub("\\(|\\)", "", names(se)), sep=".")
    df.resid <- as.data.frame(fit$df.resid)
    limma.rslt <- data.frame(limma.rslt, coef, se, df.resid)

    out.file <- paste0(ResultDir, project,"_DE_limma.rslts")
    write.table(file = out.file, limma.rslt[order(limma.rslt$P.Value),], quote=F,
                row.names=F, col.names=T, sep="\t")
    print(paste0("Wrote ", out.file))


    ## LIMMA WITH VOOM (INCLUDES SVS AND BIRTHYEAR)
    v <- voom(count, design, plot=F)
    vfit <- lmFit(v, design)    
    efit <- eBayes(vfit)
    voom.rslt <- topTable(efit, coef=2, n=Inf, sort.by="none")

    voom.rslt <- data.frame(gene = count$genes[,1], voom.rslt)
    coef <- as.data.frame(efit$coefficient)
    names(coef) <- paste("coef", gsub("\\(|\\)", "", names(coef)), sep=".")
    se <- as.data.frame(efit$s2.post * efit$stdev.unscaled)
    names(se) <- paste("se", gsub("\\(|\\)", "", names(se)), sep=".")
    df.resid <- as.data.frame(efit$df.resid)
    voom.rslt <- data.frame(voom.rslt, coef, se, df.resid)

    out.file <- paste0(ResultDir, project,"_DE_voom.rslts")
    write.table(file = out.file, voom.rslt[order(voom.rslt$P.Value),], quote=F,
                row.names=F, col.names=T, sep="\t")
    print(paste0("Wrote ", out.file))

    ## LIMMA WITH VOOM (INCLUDES BIRTHYEAR ONLY)
    v <- voom(count, design.nosvs, plot=F)
    vfit <- lmFit(v, design.nosvs)    
    efit <- eBayes(vfit)
    voom.rslt <- topTable(efit, coef=2, n=Inf, sort.by="none")

    voom.rslt <- data.frame(gene = count$genes[,1], voom.rslt)
    coef <- as.data.frame(efit$coefficient)
    names(coef) <- paste("coef", gsub("\\(|\\)", "", names(coef)), sep=".")
    se <- as.data.frame(efit$s2.post * efit$stdev.unscaled)
    names(se) <- paste("se", gsub("\\(|\\)", "", names(se)), sep=".")
    df.resid <- as.data.frame(efit$df.resid)
    voom.rslt <- data.frame(voom.rslt, coef, se, df.resid)

    out.file <- paste0(ResultDir, project,"_DE_voom_nosvs.rslts")
    write.table(file = out.file, voom.rslt[order(voom.rslt$P.Value),], quote=F,
                row.names=F, col.names=T, sep="\t")
    print(paste0("Wrote ", out.file))

    
    
    ## ANALYSIZE AFTER INVERSE NORMALIZATION ##
    cpm.invnorm <- inverse_quantile_normalization(v$E)
    txt <-  paste0("summary(lm(t(cpm.invnorm) ~ group + ", paste(
        names(count$sample)[grep("SV", names(count$sample))], collapse="+"),
        " ,data=count$samples))")
    out <- eval(parse(text = txt))
    out2 <- lapply(out, function(x) coefficients(x)[2,])
    lm.rslt <- as.data.frame(do.call(rbind, out2))
    names(lm.rslt) <- c("Est", "SE", "t.val","P.Value")
    lm.rslt <- data.frame(v$genes, lm.rslt)
    names(lm.rslt)[1] <- "ENSMBL"

    out.file <- paste0(ResultDir, project,"_LM_invnorm.rslts")
    write.table(file = out.file, lm.rslt[order(lm.rslt$P.Value),], quote=F,
                row.names=F, col.names=F, sep="\t")
    
    print(paste0("Wrote ", out.file))
    
}

get.data <- function(project, DataDir){

    file <- paste0(DataDir, project, ".RDS")
    data <- readRDS(file)

    ## ADD BIRTHYEAR IF LESS THAN 10% IS MISSING VALUES ##
    add.cov = ""
    if (grep("BirthYear", names(data$samples))){

        if (mean(is.na(data$samples$BirthYear)) < 0.10){
            add.cov <- " BirthYear + "
        }
    }
    txt <- paste0("model.matrix(~group + ",add.cov,
                  paste(names(data$samples)[grep("SV", names(data$samples))], collapse="+"),
                  ", data$samples)")
    design <- eval(parse(text = txt))
    design.nosvs <- model.matrix(~group + BirthYear, data$samples)
    
    logCPM <- cpm(data, log=TRUE, prior.count=3)

    return(list(design = design, logCPM = logCPM, count = data, design.nosvs = design.nosvs))
}

inverse_quantile_normalization <- function(gct) {
        gct = t(apply(gct, 1, rank, ties.method = "average"));
        gct = qnorm(gct / (ncol(gct)+1));
        return(gct)
}


args <- commandArgs(TRUE)
project = args[1]
run_limma(project)
print(date())

finish.file = paste0(project,".finished")
system(paste0("touch ",finish.file))
