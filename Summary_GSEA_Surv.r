ResultDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Results/"
SurvResultDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DE_Surv/Results/"
DataDir <- "/gpfs/data/stranger-lab/askol/TCGA2/Data/Expression/"

## READ IN DE RESULTS ##
projects = read.table(file = "/gpfs/data/stranger-lab/askol/TCGA/TCGA_Target_projects.txt",
    header=TRUE, as.is=TRUE, sep="\t")
projects = projects[grep("TCGA",projects[,1]),1]

skip.cancers <- c("TCGA-BRCA","TCGA-OV", "TCGA-CESC", "TCGA-PRAD", "TCGA-TGCT",
                  "TCGA-UCS", "TCGA-UCEC")
projects <- projects[!projects %in% skip.cancers] 

tmp <- collect.de.results(projects, ResultDir)

logFC.tcga <- tmp$logFC
logPs.tcga <- tmp$logPs
Qs.tcga <- tmp$Qs

## READ IN GSEA RESULTS ##
file <- paste0(ResultDir,"GSEA_Cluster_Genes_NoX_gt6shared.txt")
## LEADING EDGE GENES FROM GENE SETS ASSOCIATED WITH HIGH-LEVEL CLUSTER ASSIGNMENTS
## ONLY INCLUDES GENES THAT WERE PRESENT IN 25% OR MORE OF THE GENESETS BELONGING TO
## THE HL CLUSTER
cl.genes.25 <- read.table(file = file, header=T, as.is=T)

## CLUSTER DEFINITIONS
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

## READ IN SURVIVAL RESULTS ##
sdata <- collect.surv.data(SurvResultDir)
sdata <- mutate(sdata, DEPop = ifelse(grepl("TCGA",tissue), "TCGA" , "GTEx"))

## survival results used are thos with suffix: _lcpm.invnorm.covs_<sex>_COEF.txt ##


## PLOT SIGNED LOG P-VALUES VERSUS DE
## COLOR CODE FOR CANCER, SHAPE FOR ALL, FEMALE, MALE
## AND SIZE FOR SURVIVAL P-VALUE

## FOR EACH CANCER AND EACH HIGH LEVEL FUNCTION, USE A FISHER'S LOG SUM P-VALUE TO
## ASSESS HOW MUCH A HIGH LEVEL FUNCTION IS AFFECTING SURVIVAL

