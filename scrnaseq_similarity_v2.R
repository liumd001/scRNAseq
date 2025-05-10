######################################################################
###
###
###   Analysis of prostate cancer scRNAseq data, version 2
###     
###   2023-05-22
### This version is made for the request that 100 random genes used as 
### the control of Jaccard similarity between selected genes and STEAP2
###         
######################################################################



#######################################
###   Data and settings
#######################################

require(plotly, quietly = TRUE)
require(GGally, quietly = TRUE)
require(pheatmap, quietly = TRUE)
require(compiler, quietly = TRUE)
require(boot, quietly = TRUE)
require(parallel, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(lme4, quietly = TRUE)
require(VennDiagram)
require(ggplot2, quietly = TRUE)
require(ggven, quietly = TRUE)
require(flextable, quietly = TRUE)
require(DESeq2, quietly = TRUE)
require(GGally, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(Seurat, quietly = TRUE)
require(ggvenn, quietly = TRUE)
require(cowplot, quietly = TRUE)
require(patchwork)
require(biomaRt)

set.seed(1234)
out_dir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq"
if (!exists(out_dir)) dir.create(out_dir)


select = dplyr::select
rename = dplyr::rename
recode = dplyr::recode

hk = c("ACTB", "PGK1","GAPDH","B2M","HPRT1","GUSB","SDHA","RPL13A")

emarkers = data.frame(fullname = c("Prostate cancer-associated protein 6/ P501S/ prostein",
                                   "Prostatic Acid Phosphatase",
                                   #"prostate-specific androgen-regulated transcription factor",
                                   "prostate specific antigen",
                                   "prostate cancer antigen 3",
                                   "alpha-methylacyl-CoA racemase"),
                      gene = c("SLC45A3","ACPP","KLK3","PCA3","AMACR"))

tgenes = c("STEAP1","FOLH1","KLK2","CD276","CD46","TMEFF2","TACSTD2","STEAP2")

scale_y_log2 = function (...){
  require(scales)
  scale_y_continuous(..., trans = log_trans(2))
}

scale_x_log2 = function (...){
  require(scales)
  scale_x_continuous(..., trans = log_trans(2))
}

##########################################################
##  Choose random genes whose expression levels are
##    comparable to STEAP2
##########################################################
if(FALSE){
  load("U:\\team_share\\tfls\\2022_12_31_scRNASeq/dta_dong.RData", verbose = T)
  sum_dong = apply(dta_dong[,10:29082], 2, function(x) sum(as.numeric(x), na.rm = TRUE))
  sum_dong = as.data.frame(sum_dong) %>%
    rownames_to_column(var = "ensembl_gene_id")
  write.csv(sum_dong, "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\data\\sumCounts_byGene_dong.csv", row.names = FALSE)
  sum_dong = read.csv("U:\\team_share\\tfls\\2023_01_27_scRNASeq\\data\\sumCounts_byGene_dong.csv")
  
  load("U:\\team_share\\tfls\\2022_12_31_scRNASeq/dta_tuong.RData", verbose = T)
  sum_tuong = apply(dta_tuong[,2:22394], 2, function(x) sum(as.numeric(x), na.rm = TRUE))
  sum_tuong = as.data.frame(sum_tuong) %>%
    rownames_to_column(var = "symbol")
  write.csv(sum_tuong, "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\data\\sumCounts_byGene_tuong.csv", row.names = FALSE)
  sum_tuong = read.csv("U:\\team_share\\tfls\\2023_01_27_scRNASeq\\data\\sumCounts_byGene_tuong.csv")
  
  load("U:\\team_share\\tfls\\2022_12_31_scRNASeq/dta_chen.RData", verbose = T)
  sum_chen = apply(dta_chen[,2:23555], 2, function(x) sum(as.numeric(x), na.rm = TRUE))
  sum_chen = as.data.frame(sum_chen) %>%
    rownames_to_column(var = "symbol")
  write.csv(sum_chen, "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\data\\sumCounts_byGene_chen.csv", row.names = FALSE)
  sum_chen = read.csv("U:\\team_share\\tfls\\2023_01_27_scRNASeq\\data\\sumCounts_byGene_chen.csv")
  
  load("U:\\team_share\\tfls\\2022_12_31_scRNASeq/dta_song.RData", verbose = T)
  sum_song = apply(dta_song[,2:21878], 2, function(x) sum(as.numeric(x), na.rm = TRUE))
  sum_song = as.data.frame(sum_song) %>%
    rownames_to_column(var = "symbol")
  write.csv(sum_song, "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\data\\sumCounts_byGene_song.csv", row.names = FALSE)
  sum_song = read.csv("U:\\team_share\\tfls\\2023_01_27_scRNASeq\\data\\sumCounts_byGene_song.csv")
}#ENDOFIFFALSE

sum_dong = read.csv("U:\\team_share\\tfls\\2023_01_27_scRNASeq\\data\\sumCounts_byGene_dong.csv") 
sum_chen = read.csv("U:\\team_share\\tfls\\2023_01_27_scRNASeq\\data\\sumCounts_byGene_chen.csv")
sum_song = read.csv("U:\\team_share\\tfls\\2023_01_27_scRNASeq\\data\\sumCounts_byGene_song.csv")
sum_tuong = read.csv("U:\\team_share\\tfls\\2023_01_27_scRNASeq\\data\\sumCounts_byGene_tuong.csv")


# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# query biomart
genelist <- getBM(attributes = c("ensembl_gene_id",  "hgnc_symbol"),
                  filters = "ensembl_gene_id", 
                  values = sum_dong$ensembl_gene_id,
                  mart = mart)

sum_dong  = sum_dong %>%
  left_join(genelist) %>%
  rename(symbol = hgnc_symbol, 
         sum = sum_dong)

sum_song = sum_song %>%
  rename(sum = sum_song)

sum_chen = sum_chen %>%
  rename(sum = sum_chen)

sum_tuong = sum_tuong %>%
  rename(sum = sum_tuong)

#get threshold, differ by 2 fold and centered on STEAP2
getThreshold = function(dta){
  th = dta %>%
    filter(symbol == "STEAP2") %>%
    pull(sum)
  th = c(th/2, th, th*2)
  return(th)
}

g1_dong = sum_dong %>%
  ggplot(aes(x = sum)) +
  geom_histogram() +
  scale_x_log2() + 
  geom_vline(xintercept = getThreshold(sum_dong), color = "blue") +
  theme_bw()

g1_song = sum_song %>%
  ggplot(aes(x = sum)) +
  geom_histogram() +
  scale_x_log2() + 
  geom_vline(xintercept = getThreshold(sum_song), color = "blue") +
  theme_bw()

g1_chen = sum_chen %>%
  ggplot(aes(x = sum)) +
  geom_histogram() +
  scale_x_log2() + 
  geom_vline(xintercept = getThreshold(sum_chen), color = "blue") +
  theme_bw()

g1_tuong = sum_tuong %>%
  ggplot(aes(x = sum)) +
  geom_histogram() +
  scale_x_log2() + 
  geom_vline(xintercept = getThreshold(sum_tuong), color = "blue") +
  theme_bw()

cowplot::plot_grid(g1_dong, g1_chen, g1_song, g1_tuong, ncol = 2)


#get genes between thresold
getGenes = function(dta){
  th = getThreshold(dta)
  genes = dta$symbol[dta$sum >=th[1] & dta$sum <= th[3]]
  return(genes)
  
}

rgenes = unlist(lapply(list(sum_chen,sum_song, sum_dong, sum_tuong), getGenes))
rgenes = as.data.frame(table(rgenes)) %>%
  rownames_to_column(var = "symbol") %>%
  arrange(desc(Freq)) %>%
  filter(Freq == 3) %>%
  pull(rgenes)
rgenes = as.character(rgenes)

geneann = data.frame(symbol = c(hk,emarkers$gene, tgenes, rgenes),
                     genecat = c(rep("HK", length(hk)), 
                                 rep("PRAD_marker", length(emarkers$gene)),
                                 rep("prostate",length(tgenes)),
                                 rep("random",length(rgenes))))

geneann = geneann %>%
  left_join(genelist, by = c("symbol" = "hgnc_symbol")) %>%
  filter(!is.na(ensembl_gene_id)) %>%
  filter(!(symbol == "TACSTD2"& genecat == "random")) %>%
  filter(!(symbol == "FOLH1"& genecat == "random")) 

save(geneann, file = "U:\\team_share\\tfls\\2022_12_31_scRNASeq/data/geneann.RData")


###########################################################
##  subsetting data using selected genes
############################################################

load("U:\\team_share\\tfls\\2022_12_31_scRNASeq/data/geneann.RData", verbose = TRUE)
sampledgenes = geneann

if (FALSE){
  load("U:\\team_share\\tfls\\2022_12_31_scRNASeq/dta_dong.RData", verbose = T)
  all(sampledgenes$ensembl_gene_id %in% colnames(dta_dong))
  
  dong1 = dta_dong[,c(colnames(dta_dong)[1:9],sampledgenes$ensembl_gene_id)]
  dong1 = dong1 %>% 
    mutate(sbj = paste0("P",sbj)) %>%
    pivot_longer(cols = 10:ncol(dong1), names_to = "ensembl_gene_id", values_to = "count") %>%
    distinct() %>%
    mutate(count = as.numeric(count)) %>%
    left_join(geneann %>%
                select(ensembl_gene_id, symbol)) %>%
    select(-ensembl_gene_id) %>%
    mutate(count = as.numeric(count)) %>%
    pivot_wider(names_from = "symbol", values_from = count)
  
  rm(dta_dong)
  
  load("U:\\team_share\\tfls\\2022_12_31_scRNASeq/dta_tuong.RData", verbose = T)
  all(sampledgenes$symbol %in% colnames(dta_tuong))
  
  tuong1 = dta_tuong[,c("cell_id",sampledgenes$symbol, colnames(dta_tuong)[22395:22449])]
  rm(dta_tuong)
  
  load("U:\\team_share\\tfls\\2022_12_31_scRNASeq/dta_chen.RData", verbose = T)
  all(sampledgenes$symbol %in% colnames(dta_chen))
  chen1 = dta_chen[,c("cell_id",sampledgenes$symbol, colnames(dta_chen)[23557:23609])]
  rm(dta_chen)
  
  load("U:\\team_share\\tfls\\2022_12_31_scRNASeq/dta_song.RData", verbose = T)
  all(sampledgenes$Symbol %in% colnames(dta_song))
  song1 = dta_song[,c("cell_id",sampledgenes$Symbol, colnames(dta_song)[21879:21889])]
  rm(dta_song)
  
  save(tuong1, chen1,song1, dong1, file = "U:\\team_share\\tfls\\2022_12_31_scRNASeq/subdta_sgenes_v2.RData")
}

#save(tuong1, chen1, song1,dong1, file = "U:\\team_share\\tfls\\2022_12_31_scRNASeq/subdta_sgenes.RData")
load( "U:\\team_share\\tfls\\2022_12_31_scRNASeq/subdta_sgenes_v2.RData", verbose = TRUE)

#filter for epithelial cells

#luminal epithelial cells KLK3/KLK4 for tuong1
tuong1a  = tuong1 %>%
  filter(Author.s.sub.cell.type %in% c("luminal epithelial - KLK3","luminal epithelial - KLK4"))

# luminal cell of prostate epithelium for chen1
chen1a = chen1 %>%
  filter(celltype2 == "luminal cell of prostate epithelium")

#ERGneg_Tumor/ERGpos_Tumor for song1
song1a = song1 %>%
  filter(ID %in% c("ERGneg_Tumor", "ERGpos_Tumor"))

#Epithelial cell for dong1
dong1a = dong1 %>%
  filter(CellType %in% c("Epithelial cell"))



#combinatioins of sgenes
glist = t(combn(sampledgenes$symbol, 2))
colnames(glist) = c("g1","g2")

glist_ensembl = glist %>%
  as.data.frame() %>%
  left_join(geneann %>%
              as.data.frame() %>%
              select(ensembl_gene_id, symbol), by = c("g1" = "symbol")) %>%
  rename(g1en = ensembl_gene_id) %>%
  left_join(geneann %>%
              as.data.frame() %>%
              select(ensembl_gene_id, symbol), by = c("g2" = "symbol")) %>%
  rename(g2en = ensembl_gene_id)  %>%
  select(g1en, g2en)

lowerfun <- function(data, mapping) {
  ggplot(data = data, mapping = mapping)+ 
    geom_point(alpha = .25) + 
    geom_smooth(method = "lm", formula = y ~ x, 
                fill = "blue", color = "red", size = 0.5)
}

getDropout = function(ddta, glist1, glist2, byvar){
  
  #check gene list
  if (length(glist1) != length(glist2)) warning("Gene List 1 and Gene List2 are not in same length!")
  
  dropout = data.frame(g1 = character(),
                       g2 = character(),
                       byvar = character(),
                       do2 = numeric(),
                       do1a = numeric(),
                       do1b = numeric(),
                       do0 = numeric())
  
  for (i in 1:length(glist1)){
    gene1 = glist1[i]
    gene2 = glist2[i]
    dta2 = ddta %>%
      as.data.frame() %>%
      dplyr::select(cell_id, !!gene1, all_of(glist2[i]), all_of(byvar)) %>%
      rename(byvar = !!byvar,
             g1 = !!glist1[i],
             g2 = !!glist2[i]) %>%
      mutate(g1 = g1 == 0,
             g2 = g2 == 0) %>%
      mutate(both = g1 + g2)
    
    tmp = dta2 %>%
      group_by(byvar) %>%
      summarise(do2 = sum(both == 2),
                do1a = sum(g1 & !g2),
                do1b = sum(!g1 & g2),
                do0 = sum(both == 0)) %>%
      mutate(g1 = glist1[i],
             g2 = glist2[i]) %>%
      relocate(g1,.before = byvar) %>%
      relocate(g2,.after = g1)
    
    dropout = bind_rows(dropout,tmp)
  }
  
  dropout = dropout %>%
    mutate(pct_do2 = 100*do2/(do2 +do1a + do1b + do0),
           pct_do1 = 100*(do1a +do1b)/(do2 + do1a + do1b + do0),
           pct_doa = 100*(do1a +do2)/(do1a + do2 + do0),
           pct_dob = 100*(do1b +do2)/(do1b + do2 + do0)) %>%
    rename(!!byvar:= byvar)
  
  return(dropout)
}


transformCounts = function(d){
  d = ifelse(is.na(d), NA, ifelse(d == 0,"NEG","POS"))
  return(d)
}

getVennDiagram = function(data, subs,filter, complist,ds, outdir, title){
  
  #subsetting
  if (!is.null(subs)) data = data %>% filter(eval(parse(text = subs)))
  
  #filtering
  if (!is.null(filter)) data  = data %>% filter(eval(parse(text = filter)))
  
  if (length(complist) == 3){
    p1 = data %>%
      ggplot(aes_string(A = complist[1],
                        B = complist[2],
                        C = complist[3])) +
      geom_venn(text_size = 2, 
                set_name_size = 3,
                stroke_size = 0.4,
                fill_alpha = 0.2,
                text_color = "blue",
                auto_scale = FALSE) +
      theme_bw() +
      theme(strip.text = element_text(size = 8),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank()) +
      coord_fixed() +
      ggtitle(paste0("Venn Diagram with ",ds," cells")) +
      facet_wrap(~ Subject, ncol = 5) 
    
    temp_file = tempfile(fileext = ".png")
    ggsave(temp_file,p1, width = 12, height = 10, unit = "in")
    file.copy(temp_file,file.path(outdir,paste0("VennDiagram_",ds,"_", length(complist),"genes.png")), overwrite = TRUE)
    
    for (id in unique(data$Subject)){
      p2 = data %>%
        filter(Subject == id) %>%
        ggplot(aes_string(A = complist[1],
                          B = complist[2],
                          C = complist[3])) +
        geom_venn(text_size = 6, 
                  set_name_size = 7,
                  stroke_size = 0.4,
                  fill_alpha = 0.2,
                  text_color = "blue",
                  auto_scale = FALSE) +
        theme_bw() +
        theme(strip.text = element_text(size = 8),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank()) +
        coord_fixed() +
        ggtitle(paste0("Venn Diagram with ",ds," cells, Subject: ",id))
      
      temp_file = tempfile(fileext = ".png")
      ggsave(temp_file,p2, width = 12, height = 10, unit = "in")
      file.copy(temp_file,file.path(outdir,paste0("VennDiagram_",ds,"_", length(complist),"genes_",id,".png")), overwrite = TRUE)
    }
    
  }else if (length(complist) == 2){
    p1 = data %>%
      ggplot(aes_string(A = complist[1],
                        B = complist[2])) +
      geom_venn(text_size = 2, 
                set_name_size = 3,
                stroke_size = 0.4,
                fill_alpha = 0.2,
                text_color = "blue",
                auto_scale = FALSE) +
      theme_bw() +
      theme(strip.text = element_text(size = 8),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank()) +
      coord_fixed() +
      ggtitle(paste0("Venn Diagram with ",ds," cells")) +
      facet_wrap(~ Subject, ncol = 5) 
    
    temp_file = tempfile(fileext = ".png")
    ggsave(temp_file,p1, width = 12, height = 10, unit = "in")
    file.copy(temp_file,file.path(outdir,paste0("VennDiagram_",ds,"_", length(complist),"genes.png")), overwrite = TRUE)
    
    for (id in unique(data$Subject)){
      p2 = data %>%
        filter(Subject == id) %>%
        ggplot(aes_string(A = complist[1],
                          B = complist[2])) +
        geom_venn(text_size = 6, 
                  set_name_size = 7,
                  stroke_size = 0.4,
                  fill_alpha = 0.2,
                  text_color = "blue",
                  auto_scale = FALSE) +
        theme_bw() +
        theme(strip.text = element_text(size = 8),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank()) +
        coord_fixed() #+
      #ggtitle(paste0("Venn Diagram with ",ds," cells, Subject: ",id))
      
      #add hk genes
      p3 = data %>%
        filter(Subject == id) %>%
        ggplot(aes(A = GAPDH,
                   B = RPL13A)) +
        geom_venn(text_size = 6, 
                  set_name_size = 7,
                  stroke_size = 0.4,
                  fill_alpha = 0.2,
                  text_color = "blue",
                  auto_scale = FALSE) +
        theme_bw() +
        theme(strip.text = element_text(size = 8),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank()) +
        coord_fixed() 
      
      p4 = plot_grid(p2, p3, ncol = 2,labels = paste0("Venn Diagram with ",ds," cells, Subject: ",id))
      
      
      temp_file = tempfile(fileext = ".png")
      ggsave(temp_file,p4, width = 12, height = 10, unit = "in")
      file.copy(temp_file,file.path(outdir,paste0("VennDiagram_",ds,"_", length(complist),"genes_",id,".png")), overwrite = TRUE)
    }
  }
}

 
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

#######################################################################
##    Jaccard similarity by data sources
#######################################################################
dta_similarity = lapply(1:length(sampledgenes$symbol), function(idx) {
  var1 = sampledgenes$symbol[idx]
  t1 = lapply(1:length(sampledgenes$symbol), function(j) {
    var2 = sampledgenes$symbol[j]
    if (var1 != var2){
      if (all(c(var1, var2) %in% colnames(tuong1a))){
        dta1 = tuong1a %>%
          rename(g1 = all_of(var1),
                 g2 = all_of(var2)) %>%
          select(cell_id, g1, g2, sbj) %>%
          mutate(g1 = transformCounts(g1),
                 g2 = transformCounts(g2)) %>%
          mutate(ds = "Tuong") 
        }
      
      if (all(c(var1, var2) %in% colnames(chen1a))){
        dta1 = dta1 %>%
          bind_rows(chen1a %>%
                    rename(g1 = all_of(var1),
                           g2 = all_of(var2)) %>%
                    select(cell_id, g1, g2,sbj) %>%
                    mutate(g1 = transformCounts(g1),
                           g2 = transformCounts(g2),
                           sbj = as.character(sbj)) %>%
                    mutate(ds = "Chen") ) }
      
      if (all(c(var1, var2) %in% colnames(song1a))){
        dta1 = dta1 %>%
          bind_rows(song1a %>%
                      rename(g1 = all_of(var1),
                             g2 = all_of(var2),
                             sbj = sample) %>%
                      select(cell_id, g1, g2,sbj) %>%
                      mutate(g1 = transformCounts(g1),
                             g2 = transformCounts(g2)) %>%
                      mutate(ds = "Song") )
        }
      
      if (all(c(var1, var2) %in% colnames(song1a))){
        dta1 = dta1 %>%
          bind_rows(dong1a %>%
                      rename(g1 = all_of(var1),
                             g2 = all_of(var2)) %>%
                      select(cell_id, g1, g2,sbj) %>%
                      mutate(g1 = transformCounts(g1),
                             g2 = transformCounts(g2)) %>%
                      mutate(ds = "Dong") )
      }
      
      tbl = dta1 %>%
        group_by(ds) %>%
        summarise(dual_neg = sum(g1 == "NEG" & g2 == "NEG"),
                  single_pos = sum(g1 == "NEG" & g2 == "POS") + sum(g1 == "POS" & g2 == "NEG") ,
                  dual_pos = sum(g1 == "POS" & g2 == "POS"),
                  ds = unique(ds)) %>%
        ungroup() %>%
        mutate(jaccard = dual_pos/(single_pos + dual_pos),
               concord_overall = (dual_pos + dual_neg)/(dual_pos + dual_neg + single_pos),
               concord_pos = dual_pos /(dual_pos + dual_neg + single_pos),
               concord_neg = dual_neg/(dual_pos + dual_neg + single_pos)) %>%
        mutate(g1 = var1, g2 = var2)
      
      return(tbl)}
  })
  t1 = do.call(rbind, t1)
  return(t1)
})

dta_similarity = do.call(rbind, dta_similarity) %>%
  filter(!g1 %in% c("AMACR","PCA3")) %>%
  filter(!g2 %in% c("AMACR","PCA3")) 
  
write.csv(dta_similarity,file.path(out_dir,'data',"JaccardMatrix_byDS_v2.csv"), row.names = FALSE)

dta3 = dta_similarity %>%
  filter(g1 %in% c(emarkers$gene, hk, tgenes)) %>%
  #filter(g2 %in% c(emarkers$gene, hk, tgenes)) %>%
  left_join(geneann %>%
              select(symbol, genecat), by = c("g2" = "symbol")) %>%
  distinct() %>%
  left_join(geneann %>%
              select(symbol, genecat) %>%
              rename(genecat2 = genecat), by = c("g1" = "symbol")) %>%
  distinct() %>%
  mutate(genecat2 = recode(genecat2, "HK" = "red","PRAD_marker" = "green","random" = "grey", "prostate" = "blue")) 

a = dta3 %>% 
  select(g1, genecat2) %>%
  arrange(g1)  %>%
  distinct() %>%
  pull(genecat2)
 

dta4 = dta3 %>%
  filter(g2 %in% c(hk, tgenes,emarkers$gene))

#jaccard
g1 = dta3 %>% 
  filter(!g2 %in% c(hk, tgenes,emarkers$gene)) %>%
  ggplot(aes(x = g1, y = jaccard, color = genecat)) +
  geom_point(position = position_jitter(0.2),alpha = 0.2) +
  geom_point(data = dta4, aes(x = g1, y = jaccard, color = genecat), position = position_jitter(0.2)) +
  theme_bw() +
  #scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.2, hjust = 0, color = a),
        legend.position = "bottom") +
  labs(x = "Gene 1", y = "Jaccard Similarity", color = "Gene 2 Category") + 
  facet_wrap(~ ds)
g1
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g1, width = 12, height = 10, units = "in")
file.copy(temp_file, file.path(out_dir, "Scatterplot_Jaccard_overview_v2.png"), overwrite = TRUE)

# overall concordance
g2 = dta3 %>% 
  filter(!g2 %in% c(hk, tgenes,emarkers$gene)) %>%
  ggplot(aes(x = g1, y = jaccard, color = genecat)) +
  geom_point(position = position_jitter(0.2),alpha = 0.3) +
  geom_point(data = dta4, aes(x = g1, y = concord_overall, color = genecat), position = position_jitter(0.2)) +
  theme_bw() +
  #scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.2, hjust = 0, color = a),
        legend.position = "bottom") +
  labs(x = "Gene 1", y = "Overall Concordance", color = "Gene 2 Category") + 
  facet_wrap(~ ds)
g2
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g2, width = 12, height = 10, units = "in")
file.copy(temp_file, file.path(out_dir, "Scatterplot_OverallConcordance_overview_v2.png"), overwrite = TRUE)

# positive concordance
g3 = dta3 %>% 
  filter(!g2 %in% c(hk, tgenes,emarkers$gene)) %>%
  ggplot(aes(x = g1, y = jaccard, color = genecat)) +
  geom_point(position = position_jitter(0.2),alpha = 0.3) +
  geom_point(data = dta4, aes(x = g1, y = concord_pos, color = genecat), position = position_jitter(0.2)) +
  theme_bw() +
  #scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.2, hjust = 0, color = a),
        legend.position = "bottom") +
  labs(x = "Gene 1", y = "Positive Concordance", color = "Gene 2 Category") + 
  facet_wrap(~ ds)
g3
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g3, width = 12, height = 10, units = "in")
file.copy(temp_file, file.path(out_dir, "Scatterplot_PositiveConcordance_overview_v2.png"), overwrite = TRUE)


## dropout by genes
dta5 = apply(tuong1a[,2:121],2, function(x) round(100*sum(x == 0)/nrow(tuong1),1)) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  rename(Tuong = ".") %>%
  left_join(apply(chen1a[,2:121],2, function(x) round(100*sum(x == 0)/nrow(chen1),1)) %>%
              as.data.frame() %>%
              rownames_to_column(var = "gene") %>%
              rename(Chen = "."), by = "gene") %>%
  left_join(apply(song1a[,2:121],2, function(x) round(100*sum(x == 0)/nrow(song1),1)) %>%
              as.data.frame() %>%
              rownames_to_column(var = "gene") %>%
              rename(Song = "."), by = "gene") %>%
  left_join(apply(dong1a[,10:129],2, function(x) round(100*sum(x == 0)/nrow(dong1),1)) %>%
              as.data.frame() %>%
              rownames_to_column(var = "gene") %>%
              rename(Dong = "."), by = "gene")  %>%
  arrange(gene)

dta5 %>%
  filter(gene %in% c(tgenes, emarkers$gene)) %>%
  flextable() %>%
  #set_caption("Dropout Rate (%) by Data Source") %>%
  bg(bg = "white", part = "all") %>%
  save_as_image(path = file.path(out_dir, "Table_dropout.png"))

#######################################################################
##    Jaccard similarity by subject (add STEAP2)
#######################################################################
# Jaccard similarity already saved
if (FALSE){

  dta_similarity_bySubject = lapply(1:length(sampledgenes$symbol), function(idx) {
    var1 = sampledgenes$symbol[idx]
    
    t1 = lapply(1:length(sampledgenes$symbol), function(j) {
      var2 = sampledgenes$symbol[j]
      if (var1 != var2){
        if (all(c(var1, var2) %in% colnames(tuong1a))){
          dta1 = tuong1a %>%
            rename(g1 = all_of(var1),
                   g2 = all_of(var2)) %>%
            select(cell_id, g1, g2, sbj) %>%
            mutate(g1 = transformCounts(g1),
                   g2 = transformCounts(g2)) %>%
            mutate(ds = "Tuong") 
        }
        
        if (all(c(var1, var2) %in% colnames(chen1a))){
          dta1 = dta1 %>%
            bind_rows(chen1a %>%
                        rename(g1 = all_of(var1),
                               g2 = all_of(var2)) %>%
                        select(cell_id, g1, g2,sbj) %>%
                        mutate(g1 = transformCounts(g1),
                               g2 = transformCounts(g2),
                               sbj = as.character(sbj)) %>%
                        mutate(ds = "Chen") ) }
        
        if (all(c(var1, var2) %in% colnames(song1a))){
          dta1 = dta1 %>%
            bind_rows(song1a %>%
                        rename(g1 = all_of(var1),
                               g2 = all_of(var2),
                               sbj = sample) %>%
                        select(cell_id, g1, g2,sbj) %>%
                        mutate(g1 = transformCounts(g1),
                               g2 = transformCounts(g2)) %>%
                        mutate(ds = "Song") )
        }
        
        if (all(c(var1, var2) %in% colnames(song1a))){
          dta1 = dta1 %>%
            bind_rows(dong1a %>%
                        rename(g1 = all_of(var1),
                               g2 = all_of(var2)) %>%
                        select(cell_id, g1, g2,sbj) %>%
                        mutate(g1 = transformCounts(g1),
                               g2 = transformCounts(g2)) %>%
                        mutate(ds = "Dong") )
        }  
        
        tbl = dta1 %>%
          group_by(sbj) %>%
          summarise(dual_neg = sum(g1 == "NEG" & g2 == "NEG"),
                    single_pos = sum(g1 == "NEG" & g2 == "POS") + sum(g1 == "POS" & g2 == "NEG") ,
                    dual_pos = sum(g1 == "POS" & g2 == "POS"),
                    ds = unique(ds)) %>%
          ungroup() %>%
          mutate(jaccard = dual_pos/(single_pos + dual_pos),
                 concord_overall = (dual_pos + dual_neg)/(dual_pos + dual_neg + single_pos),
                 concord_pos = dual_pos /(dual_pos + dual_neg + single_pos),
                 concord_neg = dual_neg/(dual_pos + dual_neg + single_pos)) %>%
          mutate(g1 = var1, g2 = var2)
        
        return(tbl)}
    })
    t1 = do.call(rbind, t1)
    return(t1)
  })
    
  dta_similarity_bySubject = do.call(rbind, dta_similarity_bySubject) 
  write.csv(dta_similarity_bySubject, file.path(out_dir, "JaccardSimilarity_BySubject_v2.csv"), row.names = FALSE)
}
#ENDOFIFELSE

dta_similarity_bySubject = read.csv(file.path(out_dir, "JaccardSimilarity_BySubject_v2.csv"))

dta3 = dta_similarity_bySubject %>%
  filter(g1 == "STEAP1") %>%
  #filter(g2 %in% c(emarkers$gene, hk, tgenes)) %>%
  left_join(geneann %>%
              select(symbol, genecat), by = c("g2" = "symbol")) %>%
  distinct() %>%
  left_join(geneann %>%
              select(symbol, genecat) %>%
              rename(genecat2 = genecat), by = c("g1" = "symbol")) %>%
  distinct() %>%
  mutate(genecat3 = recode(genecat, "hk" = "red","prostate" = "green","random" = "grey",'PRAD_marker' = "blue")) %>%
  mutate(genecat4 = ifelse(genecat == "prostate","red","blue"))

dta4 = dta3 %>%
  filter(genecat != "hk") %>%
  mutate(g2 = ifelse(genecat == "random", "random", g2)) 

a = dta4 %>% 
  group_by(g2) %>%
  summarise(jaccard = median(jaccard,na.rm = TRUE)) %>%
  ungroup() %>%
  distinct() %>%
  arrange(jaccard)
acolor = ifelse(a != "random","blue",'red')

g2t = as.character(a$g2)

dta4 = dta4 %>%
  mutate(g2 = factor(g2, levels = g2t)) 


#jaccard
g1 = dta4 %>% 
  ggplot(aes(x = g2, y = jaccard)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = ds), position = position_jitter(0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.2, hjust = 0, color = acolor),
        legend.position = "bottom") +
  labs(x = "Selected Genes", y = "Jaccard Similarity", color = "Data source") +
  ggtitle("Jaccard Similarity Between STEAP1 and Selected Genes")
g1
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g1, width = 10, height = 6, units = "in")
file.copy(temp_file, file.path(out_dir, "Scatterplot_Jaccard_BySubject.png"), overwrite = TRUE)

#facet by ds
g2 = dta4 %>% 
  ggplot(aes(x = g2, y = jaccard)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(color = "red", alpha = 0.6,position = position_jitter(0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.2, hjust = 0, color = acolor),
        legend.position = "bottom") +
  labs(x = "Selected Genes", y = "Jaccard Similarity", color = "Data source") +
  ggtitle("Jaccard Similarity Between STEAP1 and Selected Genes") +
  facet_wrap(~ ds)
g2
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g2, width = 10, height = 6, units = "in")
file.copy(temp_file, file.path(out_dir, "Scatterplot_Jaccard_BySubject_facet.png"), overwrite = TRUE)


########################################################################
##    2023-02-02
##    Scatterplot
##    STEAP1 and STEAP2 positive rate by subject
########################################################################
outdir = out_dir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\20230202"
if (!dir.exists(outdir)) dir.create(outdir)

var1 = "STEAP1"
var2 = "STEAP2"
dta2 = tuong1a %>%
  rename(g1 = all_of(var1),
         g2 = all_of(var2)) %>%
  select(cell_id, g1, g2, sbj) %>%
  mutate(g1 = transformCounts(g1),
         g2 = transformCounts(g2)) %>%
  mutate(ds = "Tuong") %>%
  bind_rows(chen1a %>%
              rename(g1 = all_of(var1),
                     g2 = all_of(var2)) %>%
              select(cell_id, g1, g2,sbj) %>%
              mutate(g1 = transformCounts(g1),
                     g2 = transformCounts(g2),
                     sbj = as.character(sbj)) %>%
              mutate(ds = "Chen") ) %>%
  bind_rows(song1a %>%
              rename(g1 = all_of(var1),
                     g2 = all_of(var2),
                     sbj = orig.ident) %>%
              select(cell_id, g1, g2,sbj) %>%
              mutate(g1 = transformCounts(g1),
                     g2 = transformCounts(g2)) %>%
              mutate(ds = "Song") ) %>%
  bind_rows(dong1a %>%
              rename(g1 = all_of(var1),
                     g2 = all_of(var2)) %>%
              select(cell_id, g1, g2,sbj) %>%
              mutate(g1 = transformCounts(g1),
                     g2 = transformCounts(g2)) %>%
              mutate(ds = "Dong") ) %>%
  filter(!grepl("_N", sbj))

dta3 = dta2 %>%
  group_by(ds,sbj) %>%
  summarise(steap1 = 100*sum(g1 == "POS"),
            steap2 = 100*sum(g2 == "POS"),
            steap12_dual = 100*sum(g1 == "POS" & g2 == "POS"),
            n = n()) %>%
  ungroup() %>%
  distinct() %>%
  mutate(steap1_pos = steap1/n,
         steap2_pos = steap2/n,
         dual_pos = steap12_dual/n)

g3 = dta3 %>%
  ggplot(aes(x = steap1_pos, y = steap2_pos)) +
  geom_point(aes(color = dual_pos)) +
  theme_bw() +
  labs(x = "STEAP1 Positivity by Subjects (%)", y = "STEAP2 Positivity by Subjects", color = "Dual Postivity (%)") +
  facet_wrap(~ ds)

g3
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g3, width = 12, height = 10, units = "in")
file.copy(temp_file, file.path(outdir, "Scatterplot_STEAP_12_Positivity.png"), overwrite = TRUE)









##################################################################
##  VennDiagram: KLK3, STEAP1, PSMA
##################################################################
#tuong, 3 genes
getVennDiagram(data = tuong1 %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = Author.s.cell.type,
                           Subject = sbj) %>%
                 filter(!grepl("_N", Subject)),
               complist = c("PSMA","STEAP1","PSA"),
               filter = NULL,
               subs = "grepl('EPithelial',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_1_17_scRNASeq\\VennDiagram\\Tuong",
               ds = "Tuong_Epithelial")

#tuong, 2 genes
getVennDiagram(data = tuong1 %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = Author.s.cell.type,
                           Subject = sbj) %>%
                 filter(!grepl("_N", Subject)),
               complist = c("PSMA","STEAP1"),
               filter = "PSA",
               subs = "grepl('EPithelial',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_1_17_scRNASeq\\VennDiagram\\Tuong",
               ds = "Tuong_Epithelial_PSA+")

#Chen, 3 genes
getVennDiagram(data = chen1 %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = celltype1,
                           Subject = sbj),
               complist = c("PSMA","STEAP1","PSA"),
               filter = NULL,
               subs = "grepl('Luminal|Basal',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_1_17_scRNASeq\\VennDiagram\\Chen",
               ds = "Chen_Epithelial")

#Chen, 2 genes
getVennDiagram(data = chen1 %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = celltype1,
                           Subject = sbj),
               complist = c("PSMA","STEAP1"),
               filter = "PSA",
               subs = "grepl('Luminal|Basal',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_1_17_scRNASeq\\VennDiagram\\Chen",
               ds = "Chen_Epithelial_PSA+")

#Song, 3 genes
getVennDiagram(data = song1 %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = ID,
                           Subject = orig.ident) %>%
                 filter(!grepl("_N", Subject)),
               complist = c("PSMA","STEAP1","PSA"),
               filter = NULL,
               subs = "grepl('LE',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_1_17_scRNASeq\\VennDiagram\\Song",
               ds = "Song_Epithelial")

#song, 2 genes
getVennDiagram(data = song1 %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = ID,
                           Subject = orig.ident) %>%
                 filter(!grepl("_N", Subject)),
               complist = c("PSMA","STEAP1"),
               filter = "PSA",
               subs = "grepl('LE',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_1_17_scRNASeq\\VennDiagram\\Song",
               ds = "Song_Epithelial_PSA+")

#Dong, 3 genes
getVennDiagram(data = dong1 %>%
                 filter(!is.na(KLK3)) %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = CellType,
                           Subject = sbj) %>%
                 filter(!grepl("_N", Subject)),
               complist = c("PSMA","STEAP1","PSA"),
               filter = NULL,
               subs = "grepl('Epithelial',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_1_17_scRNASeq\\VennDiagram\\Dong",
               ds = "Dong_Epithelial")

#Dong, 2 genes
getVennDiagram(data = dong1 %>%
                 filter(!is.na(KLK3)) %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = CellType,
                           Subject = sbj) %>%
                 filter(!grepl("_N", Subject)),
               complist = c("PSMA","STEAP1"),
               filter = "PSA",
               subs = "grepl('Epithelial',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_1_17_scRNASeq\\VennDiagram\\Dong",
               ds = "Dong_Epithelial_PSA+")


##################################################################
##  VennDiagram: KLK3, STEAP1, STEAP2
##################################################################
#tuong, 3 genes
getVennDiagram(data = tuong1 %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           STEAP2 = STEAP2 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = Author.s.cell.type,
                           Subject = sbj) %>%
                 filter(!grepl("_N", Subject)),
               complist = c("STEAP2","STEAP1","PSA"),
               filter = NULL,
               subs = "grepl('EPithelial',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\VennDiagram\\Tuong",
               ds = "Tuong_Epithelial")

#tuong, 2 genes
getVennDiagram(data = tuong1 %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           RPL13A = RPL13A > 0,
                           STEAP2 = STEAP2 > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = Author.s.cell.type,
                           Subject = sbj) %>%
                 filter(!grepl("_N", Subject)),
               complist = c("STEAP2","STEAP1"),
               filter = "PSA",
               subs = "grepl('EPithelial',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\VennDiagram\\Tuong",
               ds = "Tuong_Epithelial_PSA+")

#Chen, 3 genes
getVennDiagram(data = chen1 %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           STEAP2 = STEAP2 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = celltype1,
                           Subject = sbj),
               complist = c("STEAP2","STEAP1","PSA"),
               filter = NULL,
               subs = "grepl('Luminal|Basal',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\VennDiagram\\Chen",
               ds = "Chen_Epithelial")

#Chen, 2 genes
getVennDiagram(data = chen1 %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           STEAP2 = STEAP2 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = celltype1,
                           Subject = sbj),
               complist = c("STEAP2","STEAP1"),
               filter = "PSA",
               subs = "grepl('Luminal|Basal',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\VennDiagram\\Chen",
               ds = "Chen_Epithelial_PSA+")

#Song, 3 genes
getVennDiagram(data = song1 %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           STEAP2 = STEAP2 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = ID,
                           Subject = orig.ident) %>%
                 filter(!grepl("_N", Subject)),
               complist = c("STEAP2","STEAP1","PSA"),
               filter = NULL,
               subs = "grepl('LE',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\VennDiagram\\Song",
               ds = "Song_Epithelial")

#song, 2 genes
getVennDiagram(data = song1 %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           STEAP2 = STEAP2 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = ID,
                           Subject = orig.ident) %>%
                 filter(!grepl("_N", Subject)),
               complist = c("STEAP2","STEAP1"),
               filter = "PSA",
               subs = "grepl('LE',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\VennDiagram\\Song",
               ds = "Song_Epithelial_PSA+")

#Dong, 3 genes
getVennDiagram(data = dong1 %>%
                 filter(!is.na(KLK3)) %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           STEAP2 = STEAP2 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = CellType,
                           Subject = sbj) %>%
                 filter(!grepl("_N", Subject)),
               complist = c("STEAP2","STEAP1","PSA"),
               filter = NULL,
               subs = "grepl('Epithelial',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\VennDiagram\\Dong",
               ds = "Dong_Epithelial")

#Dong, 2 genes
getVennDiagram(data = dong1 %>%
                 filter(!is.na(KLK3)) %>%
                 summarise(PSA = KLK3 >0 ,
                           PSMA = FOLH1 > 0,
                           STEAP2 = STEAP2 > 0,
                           RPL13A = RPL13A > 0,
                           GAPDH = GAPDH > 0,
                           STEAP1 = STEAP1 > 0,
                           cell_id = cell_id,
                           Celltype = CellType,
                           Subject = sbj) %>%
                 filter(!grepl("_N", Subject)),
               complist = c("STEAP2","STEAP1"),
               filter = "PSA",
               subs = "grepl('Epithelial',data$Celltype, ignore.case = TRUE)",
               outdir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\VennDiagram\\Dong",
               ds = "Dong_Epithelial_PSA+")



###########################################################################################
##    2023-02-13
##    Jaccard Similarity between StEAP2 and other genes
###########################################################################################
if (FALSE){
var1 = "STEAP2"
dta_similarity_bySubject = lapply(1:length(sampledgenes$symbol), function(j) {
  var2 = sampledgenes$symbol[j]
  if (var1 != var2){
    dta1 = tuong1a %>%
      rename(g1 = all_of(var1),
             g2 = all_of(var2)) %>%
      select(cell_id, g1, g2, sbj) %>%
      mutate(g1 = transformCounts(g1),
             g2 = transformCounts(g2)) %>%
      mutate(ds = "Tuong") %>%
      bind_rows(chen1a %>%
                  rename(g1 = all_of(var1),
                         g2 = all_of(var2)) %>%
                  select(cell_id, g1, g2,sbj) %>%
                  mutate(g1 = transformCounts(g1),
                         g2 = transformCounts(g2),
                         sbj = as.character(sbj)) %>%
                  mutate(ds = "Chen") ) %>%
      bind_rows(song1a %>%
                  rename(g1 = all_of(var1),
                         g2 = all_of(var2),
                         sbj = orig.ident) %>%
                  select(cell_id, g1, g2,sbj) %>%
                  mutate(g1 = transformCounts(g1),
                         g2 = transformCounts(g2)) %>%
                  mutate(ds = "Song") ) %>%
      bind_rows(dong1a %>%
                  rename(g1 = all_of(var1),
                         g2 = all_of(var2)) %>%
                  select(cell_id, g1, g2,sbj) %>%
                  mutate(g1 = transformCounts(g1),
                         g2 = transformCounts(g2)) %>%
                  mutate(ds = "Dong") ) %>%
      filter(!grepl("_N", sbj))
    
    tbl = dta1 %>%
      group_by(sbj) %>%
      summarise(dual_neg = sum(g1 == "NEG" & g2 == "NEG"),
                single_pos = sum(g1 == "NEG" & g2 == "POS") + sum(g1 == "POS" & g2 == "NEG") ,
                dual_pos = sum(g1 == "POS" & g2 == "POS"),
                ds = unique(ds)) %>%
      ungroup() %>%
      mutate(jaccard = dual_pos/(single_pos + dual_pos),
             concord_overall = (dual_pos + dual_neg)/(dual_pos + dual_neg + single_pos),
             concord_pos = dual_pos /(dual_pos + dual_neg + single_pos),
             concord_neg = dual_neg/(dual_pos + dual_neg + single_pos)) %>%
      mutate(g1 = var1, g2 = var2)
    
    return(tbl)}
})
dta_similarity_bySubject = do.call(rbind, dta_similarity_bySubject) 
write.csv(dta_similarity_bySubject, file.path(out_dir, "JaccardSimilarity_BySubject_STEAP2.csv"), row.names = FALSE)
} #ENDOFIFELSE



dta_similarity_bySubject =  dta_similarity_bySubject %>%
  filter(!g1 %in% c("AMACR","PCA3")) %>%
  filter(!g2 %in% c("AMACR","PCA3")) %>%
  filter(g1 == "STEAP2")

dta3 = dta_similarity_bySubject %>%
  #filter(g2 %in% c(emarkers$gene, hk, tgenes)) %>%
  left_join(geneann %>%
              select(symbol, genecat), by = c("g2" = "symbol")) %>%
  distinct() %>%
  left_join(geneann %>%
              select(symbol, genecat) %>%
              rename(genecat2 = genecat), by = c("g1" = "symbol")) %>%
  distinct() %>%
  mutate(genecat3 = recode(genecat, "hk" = "red","PRAD_marker" = "green","random" = "grey", "prostate" = "blue")) %>%
  mutate(genecat4 = ifelse(genecat == "prostate","red","blue"))

dta4 = dta3 %>%
  filter(genecat != "hk") %>%
  mutate(g2 = ifelse(genecat == "random", "random", g2)) 

a = dta4 %>% 
  group_by(g2) %>%
  summarise(jaccard = median(jaccard,na.rm = TRUE)) %>%
  ungroup() %>%
  distinct() %>%
  arrange(jaccard)
acolor = ifelse(a != "random","blue",'red')

g2t = as.character(a$g2)

dta4 = dta4 %>%
  mutate(g2 = factor(g2, levels = g2t)) 


#jaccard
g1 = dta4 %>% 
  #filter(ds != "Song") %>%
  ggplot(aes(x = g2, y = jaccard)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = ds), position = position_jitter(0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.2, hjust = 0, color = acolor),
        legend.position = "bottom") +
  labs(x = "Selected Genes", y = "Jaccard Similarity", color = "Data source") +
  ggtitle("Jaccard Similarity Between STEAP2 and Selected Genes") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14))
g1
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g1, width = 10, height = 6, units = "in")
file.copy(temp_file, file.path(out_dir, "Scatterplot_Jaccard_BySubject_STEAP2.png"), overwrite = TRUE)

#facet by ds
g2 = dta4 %>% 
  #filter(ds != "Song") %>%
  ggplot(aes(x = g2, y = jaccard)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(color = "red", alpha = 0.6,position = position_jitter(0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.2, hjust = 0, color = acolor),
        legend.position = "bottom") +
  labs(x = "Selected Genes", y = "Jaccard Similarity", color = "Data source") +
  ggtitle("Jaccard Similarity Between STEAP2 and Selected Genes") +
  facet_wrap(~ ds)
g2
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g2, width = 10, height = 6, units = "in")
file.copy(temp_file, file.path(out_dir, "Scatterplot_Jaccard_BySubject_facet_STEAP2.png"), overwrite = TRUE)


################################################################################
##    Statisitical test for co-expression
################################################################################
outdir = out_dir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\20230216"
if (!dir.exists(outdir)) dir.create(outdir)

#function to get Jaccard by sbj and cell types
getJaccardByCellTypeBySbj2 = function(gene1, gene2, celltype){
  cat(c(gene1, gene2))
  cat('\n')
  dta1 = tuong1 %>%
    filter(KLK3 > 0) %>% #filter for PSA+
    rename(g1 = all_of(gene1),
           g2 = all_of(gene2),
           celltype = "Author.s.cell.type") %>%
    select(cell_id, g1, g2, sbj, celltype) %>%
    mutate(g1 = transformCounts(g1),
           g2 = transformCounts(g2)) %>%
    mutate(ds = "Tuong") %>%
    bind_rows(chen1 %>%
                filter(KLK3 > 0) %>% #filter for PSA+
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2),
                       celltype = celltype1) %>%
                select(cell_id, g1, g2,sbj, celltype) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2),
                       sbj = as.character(sbj)) %>%
                mutate(ds = "Chen") ) %>%
    bind_rows(song1 %>%
                filter(KLK3 > 0) %>% #filter for PSA+
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2),
                       sbj = orig.ident,
                       celltype = ID) %>%
                select(cell_id, g1, g2,sbj, celltype) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2)) %>%
                mutate(ds = "Song") ) %>%
    bind_rows(dong1 %>%
                filter(KLK3 > 0) %>% #filter for PSA+
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2),
                       celltype = CellType) %>%
                select(cell_id, g1, g2,sbj, celltype) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2)) %>%
                mutate(ds = "Dong") ) %>%
    filter(!grepl("_N", sbj))
  
  
  
  tbl = dta1 %>%
    group_by(ds,sbj, celltype) %>%
    summarise(dual_neg = sum(g1 == "NEG" & g2 == "NEG"),
              single_pos = sum(g1 == "NEG" & g2 == "POS") + sum(g1 == "POS" & g2 == "NEG") ,
              dual_pos = sum(g1 == "POS" & g2 == "POS")) %>%
    ungroup() %>%
    mutate(jaccard = dual_pos/(single_pos + dual_pos),
           concord_overall = (dual_pos + dual_neg)/(dual_pos + dual_neg + single_pos),
           concord_pos = dual_pos /(dual_pos + dual_neg + single_pos),
           concord_neg = dual_neg/(dual_pos + dual_neg + single_pos)) %>%
    mutate(g1 = gene1, g2 = gene2) 
  return(tbl)
}


#choose genes that can be found in all 4 datasets
geneann1 = geneann %>%
  filter(Symbol %in% colnames(tuong1)) %>%
  filter(Symbol %in% colnames(chen1)) %>%
  filter(Symbol %in% colnames(dong1)) %>%
  filter(Symbol %in% colnames(song1))

genepair = as.data.frame(expand.grid(g1 = c("STEAP1","STEAP2","GAPDH","KLK3"), g2 = geneann1$Symbol)) %>%
  mutate(g1 = as.character(g1),
         g2 = as.character(g2)) %>%
  filter(g1 != g2)

#function to get jaccard
getJaccard = function(gene1, gene2){
  dta1 = tuong1a %>%
    rename(g1 = all_of(gene1),
           g2 = all_of(gene2)) %>%
    select(cell_id, g1, g2, sbj) %>%
    mutate(g1 = transformCounts(g1),
           g2 = transformCounts(g2)) %>%
    mutate(ds = "Tuong") %>%
    bind_rows(chen1a %>%
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2)) %>%
                select(cell_id, g1, g2,sbj) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2),
                       sbj = as.character(sbj)) %>%
                mutate(ds = "Chen") ) %>%
    bind_rows(song1a %>%
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2),
                       sbj = orig.ident) %>%
                select(cell_id, g1, g2,sbj) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2)) %>%
                mutate(ds = "Song") ) %>%
    bind_rows(dong1a %>%
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2)) %>%
                select(cell_id, g1, g2,sbj) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2)) %>%
                mutate(ds = "Dong") ) %>%
    filter(!grepl("_N", sbj))
  
  tbl = dta1 %>%
    group_by(ds,sbj) %>%
    summarise(dual_neg = sum(g1 == "NEG" & g2 == "NEG"),
              single_pos = sum(g1 == "NEG" & g2 == "POS") + sum(g1 == "POS" & g2 == "NEG") ,
              dual_pos = sum(g1 == "POS" & g2 == "POS"),
              ds = unique(ds)) %>%
    ungroup() %>%
    mutate(jaccard = dual_pos/(single_pos + dual_pos),
           concord_overall = (dual_pos + dual_neg)/(dual_pos + dual_neg + single_pos),
           concord_pos = dual_pos /(dual_pos + dual_neg + single_pos),
           concord_neg = dual_neg/(dual_pos + dual_neg + single_pos)) %>%
    mutate(g1 = gene1, g2 = gene2) 
  return(tbl)
}

dta2 = getJaccard(gene1 = "STEAP1", gene2 = "STEAP2") %>%
  select(ds, sbj, g1, g2,jaccard) %>%
  bind_rows(getJaccard(gene1 = "STEAP1", gene2 = "KLK3") %>%
              select(sbj, g1,g2, jaccard) ) %>%
  bind_rows(getJaccard(gene1 = "STEAP1", gene2 = "AMACR") %>%
              select(ds,sbj, g1,g2, jaccard) )


kruskal.test(jaccard ~ g2, data = dta2)
pairwise.wilcox.test(dta2$jaccard, dta2$g2)

#function to get Jaccard by sbj and cell types
getJaccardByCellTypeBySbj = function(gene1, gene2, celltype){
  cat(c(gene1, gene2))
  cat('\n')
  dta1 = tuong1 %>%
    rename(g1 = all_of(gene1),
           g2 = all_of(gene2),
           celltype = "Author.s.cell.type") %>%
    select(cell_id, g1, g2, sbj, celltype) %>%
    mutate(g1 = transformCounts(g1),
           g2 = transformCounts(g2)) %>%
    mutate(ds = "Tuong") %>%
    bind_rows(chen1 %>%
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2),
                       celltype = celltype1) %>%
                select(cell_id, g1, g2,sbj, celltype) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2),
                       sbj = as.character(sbj)) %>%
                mutate(ds = "Chen") ) %>%
    bind_rows(song1 %>%
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2),
                       sbj = orig.ident,
                       celltype = ID) %>%
                select(cell_id, g1, g2,sbj, celltype) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2)) %>%
                mutate(ds = "Song") ) %>%
    bind_rows(dong1 %>%
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2),
                       celltype = CellType) %>%
                select(cell_id, g1, g2,sbj, celltype) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2)) %>%
                mutate(ds = "Dong") ) %>%
    filter(!grepl("_N", sbj))
  
  tbl = dta1 %>%
    group_by(ds,sbj, celltype) %>%
    summarise(dual_neg = sum(g1 == "NEG" & g2 == "NEG"),
              single_pos = sum(g1 == "NEG" & g2 == "POS") + sum(g1 == "POS" & g2 == "NEG") ,
              dual_pos = sum(g1 == "POS" & g2 == "POS")) %>%
    ungroup() %>%
    mutate(jaccard = dual_pos/(single_pos + dual_pos),
           concord_overall = (dual_pos + dual_neg)/(dual_pos + dual_neg + single_pos),
           concord_pos = dual_pos /(dual_pos + dual_neg + single_pos),
           concord_neg = dual_neg/(dual_pos + dual_neg + single_pos)) %>%
    mutate(g1 = gene1, g2 = gene2) 
  return(tbl)
}
kruskal.test(jaccard ~ celltype,data = jdta1)
pairwise.wilcox.test(jdta1$jaccard, jdta1$celltype)


jdta = do.call(rbind, lapply(1:nrow(genepair), function(x) {
  cat(x)
  cat("\n")
  m = getJaccardByCellTypeBySbj(gene1 = genepair[x,1], gene2 = genepair[x,2])
  return(m)
}
  ))

jdta1 = jdta %>%
  filter(single_pos != 0 &dual_pos != 0) %>%
  left_join(geneann %>% 
              select(Symbol, genecat), by = c("g1" = "Symbol")) %>%
  rename(g1cat = genecat) %>%
  left_join(geneann %>%
              select(Symbol, genecat), by = c("g2" = "Symbol") ) %>%
  rename(g2cat = "genecat") %>%
  mutate(genecat = paste(g1cat, g2cat, sep = "_"))

#Wilcoxon test for STEAP2-random and STEAP2-hk, STEAP2-STEAP1, STEAP2-KLK3
doPairWiseWilcoxon = function(var1, var2, data) {
  data = subset
  m = pairwise.wilcox.test(.data$var2)
}

jdta2 = jdta1 %>%
  filter(g1 == "STEAP2") %>%
  filter(celltype %in% c("Luminal","Epithelial cell","ERGneg_Tumor","ERGpos_Tumor","Luminal epithelial - KLK3")) %>%
  mutate(genecat = ifelse(g2 == "STEAP1","STEAP2-STEAP1", genecat))%>%
  filter(genecat != "prostate_prostate") %>%
  mutate(genecat = recode(genecat, 
                          "prostate_random" = 'STEAP2-random',
                          "prostate_hk" = 'STEAP2-hk'
                          ))




g4 = jdta1 %>%
  ggplot(aes(x = celltype, y = jaccard, color = genecat)) +
  #stat_boxplot(geom ='errorbar', width=0.25, size=0.7, coef=4, position=position_dodge(0.85)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single", width = .2)) +
  #geom_boxplot(outlier.shape = NA, varwidth = TRUE) +
  geom_point(position = position_jitterdodge(0), size = 0.5, alpha = 0.1) +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = -45, vjust = 0.2, hjust = 0.2)) +
  facet_wrap( ~ ds, scales = "free")
g4
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g4, width = 10, height = 8, unit = "in")
file.copy(temp_file, file.path(out_dir,"Boxplot_Jaccard_allcelltype.png"), overwrite = TRUE)

g5 = jdta2 %>%
  ggplot(aes(x = genecat, y = jaccard)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  labs(x = "Gene Pairs", y = "Jaccard Similarity", color = "Data Sources") +
  geom_point(position = position_jitter(0.2), aes(color = ds)) +
  ylim(0,1.2)

myfit = lapply(unique(jdta2$ds), function(m) {
  myfit = as.data.frame(pairwise.wilcox.test(x = jdta2$jaccard[jdta2$ds == m], g = jdta2$genecat[jdta2$ds == m], method = "fdr")$p.value)
  myfit = myfit %>%
    rownames_to_column(var = "comp1")  %>%
    pivot_longer(cols = c("STEAP2-hk","STEAP2-random"), names_to = "comp2", values_to = "Jaccard") %>%
    filter(!is.na(Jaccard))
  myfit$ds = m
  return(myfit)
})

t1 = as.data.frame(pairwise.wilcox.test(x = jdta2$jaccard, g = jdta2$genecat, method = "fdr")$p.value)
myfit[[5]] = t1 %>%
  rownames_to_column(var = "comp1")  %>%
  pivot_longer(cols = c("STEAP2-hk","STEAP2-random"), names_to = "comp2", values_to = "Jaccard") %>%
  filter(!is.na(Jaccard))
myfit[[5]]$ds = "Overall"
myfit = do.call(rbind, myfit)

t2 = myfit %>%
  mutate(Jaccard = ifelse(Jaccard < 0.001,"< 0.001", as.character(round(Jaccard,3)))) %>%
  pivot_wider(names_from = "comp2", values_from = "Jaccard") %>%
  relocate(ds, .before = comp1) %>%
  arrange(comp1) %>%
  as_tibble()

require(ggpmisc)
g6 = g5 + geom_table_npc(data = t2, npcx = 0.65, npcy = 1,label = list(t2))
g6

temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g6, width = 10, height = 8, unit = "in")
file.copy(temp_file, file.path(out_dir,"Boxplot_STEAP2Jaccard.png"), overwrite = TRUE)


################################################################################
##    Redo: Statisitical test for co-expression
##          using PSA+ epithelial cells
################################################################################
outdir = out_dir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq\\20230217"
if (!dir.exists(outdir)) dir.create(outdir)


#choose genes that can be found in all 4 datasets
geneann1 = geneann %>%
  filter(Symbol %in% colnames(tuong1)) %>%
  filter(Symbol %in% colnames(chen1)) %>%
  filter(Symbol %in% colnames(dong1)) %>%
  filter(Symbol %in% colnames(song1))

genepair = as.data.frame(expand.grid(g1 = c("STEAP1","STEAP2","GAPDH","FOLH1"), g2 = geneann1$Symbol)) %>%
  mutate(g1 = as.character(g1),
         g2 = as.character(g2)) %>%
  filter(g1 != g2)


#function to get Jaccard by sbj and cell types
getJaccardByCellTypeBySbj2 = function(gene1, gene2, celltype){
  cat(c(gene1, gene2))
  cat('\n')
  dta1 = tuong1 %>%
    filter(KLK3 > 0) %>% #filter for PSA+
    rename(g1 = all_of(gene1),
           g2 = all_of(gene2),
           celltype = "Author.s.cell.type") %>%
    select(cell_id, g1, g2, sbj, celltype) %>%
    mutate(g1 = transformCounts(g1),
           g2 = transformCounts(g2)) %>%
    mutate(ds = "Tuong") %>%
    bind_rows(chen1 %>%
                filter(KLK3 > 0) %>% #filter for PSA+
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2),
                       celltype = celltype1) %>%
                select(cell_id, g1, g2,sbj, celltype) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2),
                       sbj = as.character(sbj)) %>%
                mutate(ds = "Chen") ) %>%
    bind_rows(song1 %>%
                filter(KLK3 > 0) %>% #filter for PSA+
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2),
                       sbj = orig.ident,
                       celltype = ID) %>%
                select(cell_id, g1, g2,sbj, celltype) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2)) %>%
                mutate(ds = "Song") ) %>%
    bind_rows(dong1 %>%
                filter(KLK3 > 0) %>% #filter for PSA+
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2),
                       celltype = CellType) %>%
                select(cell_id, g1, g2,sbj, celltype) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2)) %>%
                mutate(ds = "Dong") ) %>%
    filter(!grepl("_N", sbj))
  
  
  
  tbl = dta1 %>%
    group_by(ds,sbj, celltype) %>%
    summarise(dual_neg = sum(g1 == "NEG" & g2 == "NEG"),
              single_pos = sum(g1 == "NEG" & g2 == "POS") + sum(g1 == "POS" & g2 == "NEG") ,
              dual_pos = sum(g1 == "POS" & g2 == "POS")) %>%
    ungroup() %>%
    mutate(jaccard = dual_pos/(single_pos + dual_pos),
           concord_overall = (dual_pos + dual_neg)/(dual_pos + dual_neg + single_pos),
           concord_pos = dual_pos /(dual_pos + dual_neg + single_pos),
           concord_neg = dual_neg/(dual_pos + dual_neg + single_pos)) %>%
    mutate(g1 = gene1, g2 = gene2) 
  return(tbl)
}

jdta = do.call(rbind, lapply(1:nrow(genepair), function(x) {
  cat(x)
  cat("\n")
  m = getJaccardByCellTypeBySbj2(gene1 = genepair[x,1], gene2 = genepair[x,2])
  return(m)
}))

jdta1 = jdta %>%
  filter(single_pos != 0 &dual_pos != 0) %>%
  left_join(geneann %>% 
              select(Symbol, genecat), by = c("g1" = "Symbol")) %>%
  rename(g1cat = genecat) %>%
  left_join(geneann %>%
              select(Symbol, genecat), by = c("g2" = "Symbol") ) %>%
  rename(g2cat = "genecat") %>%
  mutate(genecat = paste(g1cat, g2cat, sep = "_"))


#function to specify gene1
getPlotTest = function(targetgene){
  jdta2 = jdta1 %>%
    filter(g1 == all_of(targetgene)) %>%
    filter(celltype %in% c("Luminal","Epithelial cell","ERGneg_Tumor","ERGpos_Tumor","Luminal epithelial - KLK3")) %>%
    mutate(genecat = ifelse(g2 == "STEAP1",paste0(targetgene,"-STEAP1"), genecat)) %>%
    mutate(genecat = ifelse(g2 == "FOLH1", paste0(targetgene,"-FOLH1"), genecat)) %>%
    mutate(genecat = ifelse(g2 == "STEAP2", paste0(targetgene,"-STEAP2"), genecat)) %>%
    filter(genecat != "prostate_prostate") %>%
    mutate(genecat = recode(genecat, 
                            "prostate_random" = all_of(paste0(targetgene,'-random')),
                            "prostate_hk" = all_of(paste0(targetgene,'-HK'))))
  g = jdta2 %>%
    ggplot(aes(x = genecat, y = jaccard)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    labs(x = "Gene Pairs", y = "Jaccard Similarity", color = "Data Sources") +
    geom_point(position = position_jitter(0.2), aes(color = ds)) +
    ylim(0,1.2) +
    ggtitle(paste0("Co-expression of ",tg," and selected genes in PSA+ epithelial cells"))
  
  myfit = lapply(unique(jdta2$ds), function(m) {
    myfit = as.data.frame(pairwise.wilcox.test(x = jdta2$jaccard[jdta2$ds == m], g = jdta2$genecat[jdta2$ds == m], method = "fdr")$p.value)
    myfit = myfit %>%
      rownames_to_column(var = "comp1")  %>%
      pivot_longer(cols = 2:4, names_to = "comp2", values_to = "Jaccard") %>%
      filter(!is.na(Jaccard))
    myfit$ds = m
    return(myfit)})
  
  t1 = as.data.frame(pairwise.wilcox.test(x = jdta2$jaccard, g = jdta2$genecat, method = "fdr")$p.value)
  myfit[[5]] = t1 %>%
    rownames_to_column(var = "comp1")  %>%
    pivot_longer(cols = 2:4, names_to = "comp2", values_to = "Jaccard") %>%
    filter(!is.na(Jaccard))
  myfit[[5]]$ds = "Overall"
  myfit = do.call(rbind, myfit)
  
  t2 = myfit %>%
    mutate(Jaccard = ifelse(Jaccard < 0.001,"< 0.001", as.character(round(Jaccard,3)))) %>%
    pivot_wider(names_from = "comp2", values_from = "Jaccard") %>%
    relocate(ds, .before = comp1) %>%
    arrange(comp1) %>%
    as_tibble() %>%
    pivot_longer(cols = 3:5, names_to = "comp2", values_to = "pvalue") %>%
    mutate(Comparison = paste0("(",comp1,") vs (", comp2,")")) %>%
    select(-comp1, -comp2) %>%
    filter(!is.na(pvalue)) %>%
    pivot_wider(names_from = Comparison, values_from = "pvalue")
  
  ft = flextable(t2)
    
  require(ggpmisc)
  g6 =  g + inset_element(gen_grob(ft, fit = "auto", scaling = "min"),
                    left = 0.05,
                    bottom = 0.8,
                    right = 0.99,
                    top = 0.99)
  return(g6)
}

for (tg in c("STEAP1","STEAP2","FOLH1")){
  g = getPlotTest(targetgene = tg)
  temp_file = tempfile(fileext = ".png")
  ggsave(temp_file, g, width = 10, height = 8, unit = "in")
  file.copy(temp_file, file.path(out_dir,paste0("BoxplotTest_",tg,"Jaccard.png")), overwrite = TRUE)
}


#function to specify gene1
#targetgene = "STEAP2"
getPlotTest2 = function(targetgene){
  jdta2 = jdta1 %>%
    filter(!g2 %in% c("ACPP","AMACR")) %>%
    filter(g1 == all_of(targetgene)) %>%
    filter(celltype %in% c("Luminal","Epithelial cell","ERGneg_Tumor","ERGpos_Tumor","Luminal epithelial - KLK3")) %>%
    mutate(genecat = ifelse(g2cat == "prostate",g2, genecat)) %>%
    mutate(genecat = ifelse(g2cat == "random", "random", genecat)) %>%
    #mutate(genecat = ifelse(g2 == "FOLH1", paste0(targetgene,"-FOLH1"), genecat)) %>%
    #mutate(genecat = ifelse(g2 == "STEAP2", paste0(targetgene,"-STEAP2"), genecat)) %>%
    #mutate(genecat = recode(genecat, 
    #                        "prostate_random" = all_of(paste0(targetgene,'-random')),
    #                        "prostate_hk" = all_of(paste0(targetgene,'-HK')))) %>%
    filter(g2cat != "hk") %>%
    mutate(genecat = as.factor(genecat)) %>%
    mutate(genecat = relevel(genecat,ref = "random"))
  
  g = jdta2 %>%
    ggplot(aes(x = genecat, y = jaccard)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    labs(x = "Gene Pairs", y = "Jaccard Similarity", color = "Data Sources") +
    geom_point(position = position_jitter(0.2), aes(color = ds)) +
    ylim(0,1.2) +
    ggtitle(paste0("Co-expression of ",targetgene," and selected genes in PSA+ epithelial cells"))
  
  myfit = lapply(unique(jdta2$ds), function(m) {
    myfit = as.data.frame(pairwise.wilcox.test(x = jdta2$jaccard[jdta2$ds == m], g = jdta2$genecat[jdta2$ds == m], method = "fdr")$p.value)
    myfit = myfit %>%
      rownames_to_column(var = "comp1")  %>%
      pivot_longer(cols = 2:11, names_to = "comp2", values_to = "Jaccard") %>%
      filter(!is.na(Jaccard))
    myfit$ds = m
    return(myfit)})
  
  t1 = as.data.frame(pairwise.wilcox.test(x = jdta2$jaccard, g = jdta2$genecat, method = "fdr")$p.value)
  myfit[[5]] = t1 %>%
    rownames_to_column(var = "comp1")  %>%
    pivot_longer(cols = 2:11, names_to = "comp2", values_to = "Jaccard") %>%
    filter(!is.na(Jaccard))
  myfit[[5]]$ds = "Overall"
  myfit = do.call(rbind, myfit)
  
  t2 = myfit %>%
    filter(comp2 == "random") %>%
    mutate(Jaccard = ifelse(Jaccard < 0.001,"< 0.001", as.character(round(Jaccard,3)))) %>%
    pivot_wider(names_from = "comp1", values_from = "Jaccard") %>%
    relocate(ds, .before = comp2) %>%
    #arrange(comp1) %>%
    as_tibble() %>%
    pivot_longer(cols = 3:12, names_to = "comp1", values_to = "pvalue") %>%
    #mutate(Comparison = comp2) %>%
    select(-comp2) %>%
    filter(!is.na(pvalue)) %>%
    pivot_wider(names_from = comp1, values_from = "pvalue")
  
  ft = flextable(t2)
  
  require(ggpmisc)
  g6 =  g + inset_element(gen_grob(ft, fit = "auto", scaling = "min"),
                          left = 0.05,
                          bottom = 0.8,
                          right = 0.99,
                          top = 0.99)
  return(g6)
}

g_steap2 = getPlotTest2(targetgene = "STEAP2")

temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g_steap2, width = 10, height = 8, unit = "in")
file.copy(temp_file, file.path(out_dir,paste0("BoxplotTest_","STEAP2","_Jaccard_B.png")), overwrite = TRUE)


########################################################
##    Create UMAP clustering Plots
##    2023-02-21
########################################################

library(Seurat, quietly = TRUE)
library(tidyverse, quietly = TRUE)
require(DESeq2, quietly = TRUE)
require(GGally)
require(flextable)
require(ggplot2)
require(ggven)
require(VennDiagram)
require(patchwork)
require(Matrix)


out_dir = "U:\\team_share\\tfls\\2023_01_27_scRNASeq/20230202"

#tuong
tuong <- readRDS("W:/devtm/cbd/users/bdecato/Prostate Single Cell/PMID_34936871_fig1.rds")

#normalize
tuong = NormalizeData(tuong, normalization.method = "LogNormalize", scale.factor = 10000)

#identifying highly variable features
tuong = FindVariableFeatures(tuong, selection.method = "vst", nfeatures = 200)

#scaling data
all.genes = rownames(tuong)
tuong = ScaleData(tuong, features = all.genes)

#dimensional reduction
tuong = RunPCA(tuong, features = VariableFeatures(object = tuong))

#determine dimensionality
tuong = JackStraw(tuong, num.replicate = 100)
tuong = ScoreJackStraw(tuong, dims = 1:20)

#clustering
tuong = FindNeighbors(tuong, dims = 1:10)
tuong = FindClusters(tuong, resolution = 0.5)
head(Idents(tuong),5)

#UMAP
tuong = RunUMAP(tuong, dims = 1:14)

g_umap_tuong  = tuong[['umap']]@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = "cellid") %>%
  left_join(tuong@meta.data %>%
              select(Author.s.cell.type) %>%
              rownames_to_column(var = "cellid"), by = "cellid") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = Author.s.cell.type)) +
  geom_point(size = 0.5, alpha = 0.5) +
  #geom_text(x = -Inf, y = Inf, label = "Tuong", vjust = 2,hjust = -0.5, color = "black") +
  theme_bw() 

#save(g_umap_tuong, file = file.path(out_dir,"CompositePlot.RData"))

#chen
chen <- readRDS("W:/devtm/cbd/users/bdecato/Prostate Single Cell/GSE141445.rds")

#normalize
chen = NormalizeData(chen, normalization.method = "LogNormalize", scale.factor = 10000)

#identifying highly variable features
chen = FindVariableFeatures(chen, selection.method = "vst", nfeatures = 200)

#scaling data
all.genes = rownames(chen)
chen = ScaleData(chen, features = all.genes)

#dimensional reduction
chen = RunPCA(chen, features = VariableFeatures(object = chen))

#determine dimensionality
chen = JackStraw(chen, num.replicate = 100)
chen = ScoreJackStraw(chen, dims = 1:20)

#clustering
chen = FindNeighbors(chen, dims = 1:10)
chen = FindClusters(chen, resolution = 0.5)

#UMAP
chen = RunUMAP(chen, dims = 1:14)

g_umap_chen  = chen[['umap']]@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = "cellid") %>%
  left_join(chen@meta.data %>%
              select(Author.s.cell.type) %>%
              rownames_to_column(var = "cellid"), by = "cellid") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = Author.s.cell.type)) +
  geom_point(size = 0.5, alpha = 0.5) +
  #geom_text(x = -Inf, y = Inf, label = "Chen", vjust = 2,hjust = -0.5, color = "black") +
  theme_bw() 

g_umap_chen
#save(g_umap_tuong, g_umap_chen, file = file.path(out_dir,"CompositePlot.RData"))

#Song
dge = readRDS("C:\\Users\\mliu10\\Documents\\Amgen_documents\\prostate_scRNAseq\\GSE176031\\Seurat_dataset\\dge_E.rds")
g_umap_song  = dge[['umap']]@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column(var = "cellid") %>%
  left_join(dge@meta.data %>%
              select(ID) %>%
              rownames_to_column(var = "cellid"), by = "cellid") %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = ID)) +
  labs(color = "Cell Type") +
  geom_point(size = 0.5, alpha = 0.5) +
  #geom_text(x = -Inf, y = Inf, label = "Song", vjust = 2,hjust = -0.5, color = "black") +
  theme_bw() 

g_umap_song
#save(g_umap_tuong, g_umap_chen, g_umap_song,file = file.path(out_dir,"CompositePlot.RData"))


######################################
##scatterplot
var1 = "STEAP1"
var2 = "FOLH1"
dta2 = tuong1a %>%
  rename(g1 = all_of(var1),
         g2 = all_of(var2)) %>%
  select(cell_id, g1, g2, sbj) %>%
  mutate(g1 = transformCounts(g1),
         g2 = transformCounts(g2)) %>%
  mutate(ds = "Tuong") %>%
  bind_rows(chen1a %>%
              rename(g1 = all_of(var1),
                     g2 = all_of(var2)) %>%
              select(cell_id, g1, g2,sbj) %>%
              mutate(g1 = transformCounts(g1),
                     g2 = transformCounts(g2),
                     sbj = as.character(sbj)) %>%
              mutate(ds = "Chen") ) %>%
  bind_rows(song1a %>%
              rename(g1 = all_of(var1),
                     g2 = all_of(var2),
                     sbj = orig.ident) %>%
              select(cell_id, g1, g2,sbj) %>%
              mutate(g1 = transformCounts(g1),
                     g2 = transformCounts(g2)) %>%
              mutate(ds = "Song") ) %>%
  bind_rows(dong1a %>%
              rename(g1 = all_of(var1),
                     g2 = all_of(var2)) %>%
              select(cell_id, g1, g2,sbj) %>%
              mutate(g1 = transformCounts(g1),
                     g2 = transformCounts(g2)) %>%
              mutate(ds = "Dong") ) %>%
  filter(!grepl("_N", sbj))

dta3 = dta2 %>%
  group_by(ds,sbj) %>%
  summarise(steap1 = 100*sum(g1 == "POS"),
            steap2 = 100*sum(g2 == "POS"),
            steap12_dual = 100*sum(g1 == "POS" & g2 == "POS"),
            n = n()) %>%
  ungroup() %>%
  distinct() %>%
  mutate(steap1_pos = steap1/n,
         steap2_pos = steap2/n,
         dual_pos = steap12_dual/n)

scatterplot_steap1 = dta3 %>%
  filter(ds!= "Dong") %>%
  ggplot(aes(x = steap1_pos, y = steap2_pos)) +
  geom_point(aes(color = dual_pos)) + 
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "STEAP1 Positivity by Subjects (%)", y = "PSMA Positivity by Subjects", color = "Dual Postivity (%)") +
  facet_wrap(~ ds, ncol = 3)

scatterplot_steap1


#jaccard
#choose genes that can be found in all 4 datasets
geneann1 = geneann %>%
  filter(Symbol %in% colnames(tuong1)) %>%
  filter(Symbol %in% colnames(chen1)) %>%
  filter(Symbol %in% colnames(dong1)) %>%
  filter(Symbol %in% colnames(song1))

genepair = as.data.frame(expand.grid(g1 = c("STEAP1","STEAP2","GAPDH","FOLH1"), g2 = geneann1$Symbol)) %>%
  mutate(g1 = as.character(g1),
         g2 = as.character(g2)) %>%
  filter(g1 != g2)


#function to get Jaccard by sbj and cell types
getJaccardByCellTypeBySbj2 = function(gene1, gene2, celltype){
  cat(c(gene1, gene2))
  cat('\n')
  dta1 = tuong1 %>%
    filter(KLK3 > 0) %>% #filter for PSA+
    rename(g1 = all_of(gene1),
           g2 = all_of(gene2),
           celltype = "Author.s.cell.type") %>%
    select(cell_id, g1, g2, sbj, celltype) %>%
    mutate(g1 = transformCounts(g1),
           g2 = transformCounts(g2)) %>%
    mutate(ds = "Tuong") %>%
    bind_rows(chen1 %>%
                filter(KLK3 > 0) %>% #filter for PSA+
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2),
                       celltype = celltype1) %>%
                select(cell_id, g1, g2,sbj, celltype) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2),
                       sbj = as.character(sbj)) %>%
                mutate(ds = "Chen") ) %>%
    bind_rows(song1 %>%
                filter(KLK3 > 0) %>% #filter for PSA+
                rename(g1 = all_of(gene1),
                       g2 = all_of(gene2),
                       sbj = orig.ident,
                       celltype = ID) %>%
                select(cell_id, g1, g2,sbj, celltype) %>%
                mutate(g1 = transformCounts(g1),
                       g2 = transformCounts(g2)) %>%
                mutate(ds = "Song") ) %>%
    filter(!grepl("_N", sbj))
  
  
  
  tbl = dta1 %>%
    group_by(ds,sbj, celltype) %>%
    summarise(dual_neg = sum(g1 == "NEG" & g2 == "NEG"),
              single_pos = sum(g1 == "NEG" & g2 == "POS") + sum(g1 == "POS" & g2 == "NEG") ,
              dual_pos = sum(g1 == "POS" & g2 == "POS")) %>%
    ungroup() %>%
    mutate(jaccard = dual_pos/(single_pos + dual_pos),
           concord_overall = (dual_pos + dual_neg)/(dual_pos + dual_neg + single_pos),
           concord_pos = dual_pos /(dual_pos + dual_neg + single_pos),
           concord_neg = dual_neg/(dual_pos + dual_neg + single_pos)) %>%
    mutate(g1 = gene1, g2 = gene2) 
  return(tbl)
}

if (FALSE){
  jdta = do.call(rbind, lapply(1:nrow(genepair), function(x) {
    cat(x)
    cat("\n")
    m = getJaccardByCellTypeBySbj2(gene1 = genepair[x,1], gene2 = genepair[x,2])
    return(m)
  }))
  
  jdta1 = jdta %>%
    filter(single_pos != 0 &dual_pos != 0) %>%
    left_join(geneann %>% 
                select(Symbol, genecat), by = c("g1" = "Symbol")) %>%
    rename(g1cat = genecat) %>%
    left_join(geneann %>%
                select(Symbol, genecat), by = c("g2" = "Symbol") ) %>%
    rename(g2cat = "genecat") %>%
    mutate(genecat = paste(g1cat, g2cat, sep = "_"))
  
  write.csv(jdta1, file.path(out_dir, "jaccardmatrix.csv"), row.names = FALSE)
}

jdta1 = read.csv(file.path("U:\\team_share\\tfls\\2023_01_27_scRNASeq/20230202", "jaccardmatrix.csv")) 
jdta1 = jdta1 %>%
  filter(ds != "Dong")

#function to specify gene1
#targetgene = "STEAP1"
getPlotTest2 = function(targetgene){
  jdta2 = jdta1 %>%
    filter(!g2 %in% c("ACPP","AMACR")) %>%
    filter(g1 == all_of(targetgene)) %>%
    filter(celltype %in% c("Luminal","Epithelial cell","ERGneg_Tumor","ERGpos_Tumor","Luminal epithelial - KLK3")) %>%
    mutate(genecat = ifelse(g2cat == "prostate",g2, genecat)) %>%
    mutate(genecat = ifelse(g2cat == "random", "random", genecat)) %>%
    #mutate(genecat = ifelse(g2 == "FOLH1", paste0(targetgene,"-FOLH1"), genecat)) %>%
    #mutate(genecat = ifelse(g2 == "STEAP2", paste0(targetgene,"-STEAP2"), genecat)) %>%
    #mutate(genecat = recode(genecat, 
    #                        "prostate_random" = all_of(paste0(targetgene,'-random')),
    #                        "prostate_hk" = all_of(paste0(targetgene,'-HK')))) %>%
    filter(g2cat != "hk") %>%
    mutate(genecat = as.factor(genecat)) %>%
    mutate(genecat = relevel(genecat,ref = "random"))
  
  genecat_level = jdta2 %>%
    group_by(genecat) %>%
    summarise(median = median(jaccard)) %>%
    ungroup() %>%
    arrange(median) %>%
    pull(genecat)
  
  g = jdta2 %>%
    mutate(genecat = factor(genecat, levels = genecat_level)) %>%
    arrange(genecat) %>%
    ggplot(aes(x = genecat, y = jaccard)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw() +
    labs(x = "Gene Pairs", y = "Jaccard Similarity", color = "Data Sources") +
    geom_point(position = position_jitter(0.2), aes(color = ds)) +
    ylim(0,1.2) +
    ggtitle(paste0("Co-expression of ",targetgene," and selected genes in PSA+ epithelial cells"))
  
  myfit = lapply(unique(jdta2$ds), function(m) {
    myfit = as.data.frame(pairwise.wilcox.test(x = jdta2$jaccard[jdta2$ds == m], g = jdta2$genecat[jdta2$ds == m], method = "fdr")$p.value)
    myfit = myfit %>%
      rownames_to_column(var = "comp1")  %>%
      pivot_longer(cols = 2:11, names_to = "comp2", values_to = "Jaccard") %>%
      filter(!is.na(Jaccard))
    myfit$ds = m
    return(myfit)})
  
  t1 = as.data.frame(pairwise.wilcox.test(x = jdta2$jaccard, g = jdta2$genecat, method = "fdr")$p.value)
  myfit[[5]] = t1 %>%
    rownames_to_column(var = "comp1")  %>%
    pivot_longer(cols = 2:11, names_to = "comp2", values_to = "Jaccard") %>%
    filter(!is.na(Jaccard))
  myfit[[5]]$ds = "Overall"
  myfit = do.call(rbind, myfit)
  
  t2 = myfit %>%
    filter(comp2 == "random") %>%
    mutate(Jaccard = ifelse(Jaccard < 0.001,"< 0.001", as.character(round(Jaccard,3)))) %>%
    pivot_wider(names_from = "comp1", values_from = "Jaccard") %>%
    relocate(ds, .before = comp2) %>%
    #arrange(comp1) %>%
    as_tibble() %>%
    pivot_longer(cols = 3:12, names_to = "comp1", values_to = "pvalue") %>%
    #mutate(Comparison = comp2) %>%
    select(-comp2) %>%
    filter(!is.na(pvalue)) %>%
    pivot_wider(names_from = comp1, values_from = "pvalue")
  
  ft = flextable(t2)
  
  require(ggpmisc)
  g6 =  g + inset_element(gen_grob(ft, fit = "auto",  scaling = "min"),
                          left = 0.05,
                          bottom = 0.7,
                          right = 0.99,
                          top = 0.99)
  return(list(g, ft))
}

boxplot_list = getPlotTest2(targetgene = "STEAP1")


save(g_umap_tuong, 
     g_umap_chen, 
     g_umap_song,
     scatterplot_steap1,
     boxplot_list,
     file = file.path(out_dir,"CompositePlot.RData"))


load(file.path(out_dir,"CompositePlot.RData"),verbose = TRUE)

require(cowplot)
umap_legend_fontsize = 8
umap_legend_space = 0.01
legend_title_fontsize = 10
key_height = 0.6
axis_title_fontsize = 8
axis_text_fontsize = 6
g_umap_tuong1 = g_umap_tuong +
  theme(legend.text = element_text(size = umap_legend_fontsize),
        legend.spacing.y = unit(umap_legend_space, "cm"),
        legend.title = element_blank(),
        axis.title = element_text(size = axis_title_fontsize),
        axis.text = element_text(size = axis_text_fontsize)) +
  geom_text(x = -Inf, y = Inf, label = "Tuong", color = "black", size = 3, vjust = 1.2, hjust = -0.3) +
  guides(color = guide_legend(byrow = TRUE, keyheight = key_height))
g_umap_chen1 = g_umap_chen +
  theme(legend.text = element_text(size = umap_legend_fontsize),
        legend.spacing.y = unit(umap_legend_space, "cm"),
        legend.title = element_blank(),
        axis.title = element_text(size = axis_title_fontsize),
        axis.text = element_text(size = axis_text_fontsize)) +
  geom_text(x = -Inf, y = Inf, label = "Chen", color = "black", size = 3, vjust = 1.2, hjust = -0.3) +
  guides(color = guide_legend(byrow = TRUE, keyheight = key_height))
g_umap_song1 = g_umap_song +
  theme(legend.text = element_text(size = umap_legend_fontsize),
        legend.spacing.y = unit(umap_legend_space, "cm"),
        legend.title = element_blank(),
        axis.title = element_text(size = axis_title_fontsize),
        axis.text = element_text(size = axis_text_fontsize)) +
  geom_text(x = -Inf, y = Inf, label = "Song", color = "black", size = 3, vjust = 1.2, hjust = -0.3) +
  guides(color = guide_legend(byrow = TRUE, keyheight = key_height))

scatterplot_steap11 = scatterplot_steap1 +
  theme(legend.text = element_text(size = umap_legend_fontsize),
        legend.spacing.y = unit(umap_legend_space, "cm"),
        legend.title = element_text(size = 8),
        axis.title = element_text(size = axis_title_fontsize),
        axis.text = element_text(size = axis_text_fontsize))

boxplot_steap11 = boxplot_list[[1]] +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(size = 10),
        #legend.position = 'none',
        #legend.spacing.y = unit(umap_legend_space, "cm"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8)) +
  inset_element(gen_grob(boxplot_list[[2]], fit = "auto",  scaling = "min"),
                left = 0.05,
                bottom = 0.7,
                right = 0.99,
                top = 0.99)


g1 = plot_grid(g_umap_chen1, g_umap_song1, g_umap_tuong1, ncol = 3, rel_widths = c(2.7,2.7,3))
g2 = plot_grid(g1, scatterplot_steap11, boxplot_steap11, nrow = 3, ncol = 1, rel_heights = c(3,3,7), rel_widths = c(10,8,10))
g2
temp_file = tempfile(fileext = ".png")
ggsave(temp_file, g2, width = 12, height = 12, unit = "in")
file.copy(temp_file, file.path(out_dir, "CompositePlot.png"), overwrite = TRUE)
