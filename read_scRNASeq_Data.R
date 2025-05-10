

library(Seurat, quietly = TRUE)
library(tidyverse, quietly = TRUE)
require(DESeq2, quietly = TRUE)
require(GGally)
require(flextable)
require(ggplot2)
require(ggven)
require(VennDiagram)


out_dir = "U:\\team_share\\tfls\\2022_12_31_scRNASeq"

tuong <- readRDS("W:/devtm/cbd/users/bdecato/Prostate Single Cell/PMID_34936871_fig1.rds")

tuong_row = apply(as.data.frame(t(as.matrix(tuong@assays$RNA@counts))), 1, function(x) sum(x == 0)) 
tuong_column = apply(as.data.frame(t(as.matrix(tuong@assays$RNA@counts))), 2, function(x) sum(x == 0)) 

dta_tuong = as.data.frame(t(as.matrix(tuong@assays$RNA@counts))) %>%
  rownames_to_column(var = "cell_id") %>%
  left_join(tuong@meta.data %>%
              rename(sbj = Sample.ID) %>%
              rownames_to_column(var = "cell_id"), by = "cell_id") %>%
  distinct() 
save(dta_tuong, file = file.path(out_dir, "dta_tuong.RData"))
rm(tuong)
rm(dta_tuong)

#Chen data
chen <- readRDS("W:/devtm/cbd/users/bdecato/Prostate Single Cell/GSE141445.rds")
chen_row = apply(as.data.frame(t(as.matrix(chen@assays$RNA@counts))), 1, function(x) sum(x == 0)) 
chen_column = apply(as.data.frame(t(as.matrix(chen@assays$RNA@counts))), 2, function(x) sum(x == 0)) 

dta_chen = as.data.frame(t(as.matrix(chen@assays$RNA@counts))) %>%
  rownames_to_column(var = "cell_id") %>%
  left_join(chen@meta.data %>%
              rownames_to_column(var = "cell_id"), by = "cell_id") %>%
  distinct() %>%
  rename(celltype1 = Author.s.cell.type,
         celltype2 = cell.type...subgroup..standardized...standardized.,
         sbj = Sample.ID)
save(dta_chen, file = file.path(out_dir, "dta_chen.RData"))
rm(chen)
rm(dta_chen)

#Song data
dge = readRDS("C:\\Users\\mliu10\\Documents\\Amgen_documents\\prostate_scRNAseq\\GSE176031\\Seurat_dataset\\dge_E.rds")
song_row = apply(as.data.frame(t(as.matrix(dge@assays$RNA@counts))), 1, function(x) sum(x == 0)) 
song_column = apply(as.data.frame(t(as.matrix(dge@assays$RNA@counts))), 2, function(x) sum(x == 0)) 


dta_song = dge@assays$RNA@counts %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "cell_id") %>%
  left_join(dge@meta.data %>%
              as.data.frame() %>%
              rownames_to_column(var = "cell_id"))
save(dta_song, file = file.path(out_dir, "dta_song.RData"))
rm(dge)
rm(dta_song)

#Dong data
require(data.table)
require(openxlsx)
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
