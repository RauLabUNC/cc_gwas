library(dplyr); library(tidyverse)

##BE CAREFUL WITH THE NAME
RunDescription="250602-VST_614"

#Create directory
output_dir <- file.path("Results/Consolidate_SNPs/GeneOwnChr",RunDescription)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


###CTRL ONLY--------------------------------------------------------------------
setwd(file.path("/proj/raulab/users/Anh/CC-miQTL_VST/Results/OneFiveLOD_all_SNPs/GeneOwnChr", RunDescription, "control/"))
MyFiles <- list.files(pattern="*.csv", all.files = T)
df = do.call(rbind, lapply(MyFiles, read.csv))

setwd(file.path("/proj/raulab/users/Anh/CC-miQTL_VST"))
write.csv(df, file=file.path(output_dir, "control_snps.csv"))


###ISO ONLY--------------------------------------------------------------------
setwd(file.path("/proj/raulab/users/Anh/CC-miQTL_VST/Results/OneFiveLOD_all_SNPs/GeneOwnChr", RunDescription, "iso/"))
MyFiles <- list.files(pattern="*.csv", all.files = T)
df = do.call(rbind, lapply(MyFiles, read.csv))

setwd(file.path("/proj/raulab/users/Anh/CC-miQTL_VST"))
write.csv(df, file=file.path(output_dir, "iso_snps.csv"))



