ScriptDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Summary/Scripts/"
CodeDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Code/"
SummaryDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Summary/"

setwd(SummaryDir)

projects = read.table(file = "/gpfs/data/stranger-lab/askol/TCGA/TCGA_Target_projects.txt",
    header=TRUE, as.is=TRUE, sep="\t")
projects = projects[grep("TCGA",projects[,1]),1]

skip.cancers <- c("TCGA-BRCA","TCGA-OV", "TCGA-CESC", "TCGA-PRAD", "TCGA-TGCT",
                  "TCGA-UCS", "TCGA-UCEC")
projects <- projects[!projects %in% skip.cancers] 



if (dir.exists(ScriptDir)==FALSE){
    dir.create(ScriptDir)
}

create.pbs <- function(ScriptDir, project){
    
    file <- paste0(ScriptDir, project, "_summary.pbs")
    sh.txt <- rbind(
        "#!/bin/bash",
        "#PBS -l  nodes=1:ppn=1,mem=16gb",
        "#PBS -l walltime=96:00:00",
        paste0("#PBS -o ", ScriptDir),
        "#PBS -j oe",
        paste0("#PBS -N ",project),
        "module load R")
    
    write.table(file = file, sh.txt, quote=F, row.names=F, col.names=F)

    return(file=file)
}

for (project in projects){
    
    cmd <-  paste0("R CMD BATCH  '--args ",project,
                   "' ", CodeDir, 
                   "Compare_DE_results.r ",
                   ScriptDir,
                   project,".out")
    
    print(paste0("Working on project ",project))
    
    file <- create.pbs(ScriptDir, project)
    write.table(file = file, cmd, quote=F, row.names=F, col.names=F, append=T)

    cmd <- paste0("qsub ", file)
    system(cmd)
}

