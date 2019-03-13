job.files = c()

projects = read.table(file = "/gpfs/data/stranger-lab/askol/TCGA/TCGA_Target_projects.txt",
    header=TRUE, as.is=TRUE, sep="\t")
projects = projects[grep("TCGA",projects[,1]),1]

skip.cancers <- c("TCGA-BRCA","TCGA-OV", "TCGA-CESC", "TCGA-PRAD", "TCGA-TGCT",
                  "TCGA-UCS", "TCGA-UCEC")
projects <- projects[!projects %in% skip.cancers] 


job.files = c()

ScriptDir <- "/gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Scripts/"

if (!dir.exists(ScriptDir)){
    dir.create(ScriptDir)
}
setwd(ScriptDir)


for (project in projects){

    job.file = paste0(project,".pbs")
    job.files = c(job.files, job.file)

    write(file = job.file, "#!/bin/bash")
    write(file = job.file, paste("#PBS -l nodes=1:ppn=1,mem=8gb"), append=TRUE)
    write(file=job.file, "#PBS -l walltime=36:00:00", append=TRUE)
    write(file = job.file,
          paste0("#PBS -o ", ScriptDir,
                 project,"_pbs.out"), append=TRUE)
    write(file = job.file, "#PBS -j oe", append=TRUE)
    write(file = job.file, paste0("#PBS -N Crt_dat_",project), append=TRUE)
    write(file = job.file, paste0("R CMD BATCH  '--args ",project,
              "' /gpfs/data/stranger-lab/askol/TCGA2/DiffExpression/Code/DE_limma.r ",
              ScriptDir, project,".out"),
              append=TRUE)    
    system(paste("chmod +x", job.file))
}

for (job in job.files){
    system(paste("qsub",job))
}
