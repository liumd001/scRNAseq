

library(Seurat, quietly = TRUE)
library(tidyverse, quietly = TRUE)
require(DESeq2, quietly = TRUE)
require(GGally)
require(flextable)
require(ggplot2)
require(ggven)
require(VennDiagram)

datafolder = '/Volumes/home/Drive/learning_lmd/dataset/scRNAseq/'

#tuong data from prostate cancer atlas
tuong <- readRDS("/Volumes/home/Drive/learning_lmd/dataset/scRNAseq/prostate_portal_300921_tuong.RDS")

#Chen data, #chen, gse141445
list.files("/Volumes/home/Drive/learning_lmd/dataset/scRNAseq/GSE141445")
chen = read.table(gzfile("/Volumes/home/Drive/learning_lmd/dataset/scRNAseq/GSE141445/GSM4203181_data.matrix.txt.gz"))
chen_samples = getGEO(filename = file.path(datafolder,"GSE141445_series_matrix.txt.gz"))@phenoData@data
write.csv(chen_samples, file.path(datafolder,"chen_GSE141445_samples.csv"), row.names = FALSE)

#Song data
list.files(file.path(datafolder,"GSE176031"))
untar(file.path(datafolder,"GSE176031","GSE176031_RAW.tar"), exdir = file.path(datafolder,"GSE176031"))
song_files = list.files(file.path(datafolder,"GSE176031"), full.names = TRUE)
song_files = song_files[grepl("\\.gz", song_files)]
song = lapply(song_files,fread)
song = lapply(1:length(song), function(x){
  dta = song[[x]]
  colnames(dta)[2:ncol(dta)] <- paste0(colnames(dta)[2:ncol(dta)],"_",x)
  return(dta)
})
song <- Reduce(function(x,y) merge(x, y, by = "GENE", all = TRUE), song)
song_samples = GEOquery::getGEO(filename = file.path(datafolder,"GSE176031_series_matrix.txt"))@phenoData@data
write.csv(song, file.path(datafolder,"dta_song_renamecolumns.csv"), row.names = FALSE)
write.csv(song_samples, file.path(datafolder,"dta_song_samples.csv"), row.names = FALSE)

#Dong data
require(data.table)
require(openxlsx)
untar(file.path(datafolder,"GSE137829","GSE137829_RAW.tar"), exdir = file.path(datafolder,"GSE137829"))
dong_files = list.files(file.path(datafolder,"GSE137829"), full.names = TRUE)
dong_files = dong_files[grepl("\\.gz", dong_files)]
dong = lapply(dong_files,fread)
dong <- Reduce(function(x,y) merge(x, y, by = c("Gene_ID","Symbol"), all = TRUE), dong)
dong_samples = GEOquery::getGEO(filename = file.path(datafolder,"GSE137829","GSE137829_series_matrix.txt"))@phenoData@data
write.csv(dong,file.path(datafolder,"GSE137829","dta_dong.csv"), row.names = FALSE)
write.csv(dong_samples,file.path(datafolder,"GSE137829","dta_dong_samples.csv"), row.names = FALSE)









dongfolder = "C:\\Users\\mliu10\\Documents\\Amgen_documents\\prostate_scRNAseq\\GSE137829"
dongfiles = list.files(dongfolder)
dongfiles = dongfiles[grepl("^GSM", dongfiles)]
dong = lapply(file.path(dongfolder,dongfiles), fread)

dong_cluster = read.xlsx(file.path(dongfolder, "42003_2020_1476_MOESM4_ESM.xlsx"))
dong_cluster$sbj = gsub("patient #","P",dong_cluster$orig.ident)

names(dong) = as.data.frame(do.call(rbind, lapply(dongfiles, function(x) strsplit(x, split = "_")[[1]])))[,2]
dong1 = dong

for (i in 1:length(dong1)){
  dong1[[i]]$sbj = rep(names(dong)[i], nrow(dong1[[i]]))
  dong1[[i]] = dong1[[i]] %>%
    rename(gene = Gene_ID) %>%
    t() %>%
    as.data.frame() 
  colnames(dong1[[i]]) = dong1[[i]][1,]
  dong1[[i]] = dong1[[i]] %>%
    slice(-c(1,2)) %>%
    mutate(sbj = names(dong1)[i]) %>%
    rownames_to_column(var = "cell_id")
  
  dong1[[i]] = dong_cluster %>%
    select(-orig.ident) %>%
    filter(sbj == unique(dong1[[i]]$sbj)) %>%
    rename(cell_id = X1) %>%
    left_join(dong1[[i]], by = c("cell_id","sbj"))
  }

geneann = dong[[1]][,c(1,2)]

dta_dong = dong1[[1]] %>%
  full_join(dong1[[2]]) %>%
  full_join(dong1[[3]]) %>%
  full_join(dong1[[4]]) %>%
  full_join(dong1[[5]]) %>%
  full_join(dong1[[6]])

#check if FOLH1 is there
"ENSG00000086205" %in% colnames(dta_dong)

save(dta_dong, file = file.path(out_dir, "dta_dong.RData"))
save(geneann, file = file.path(out_dir,"geneann.RData"))

rm(dong)
rm(dong1)
rm(dta_dong)

save(tuong_row, tuong_column, chen_row, chen_column, song_row, song_column, file = file.path(out_dir, "row_column_summary.RData"))
