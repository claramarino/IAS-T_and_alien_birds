# Functional space in n dim using pcoa

rm(list=ls())

library(tidyverse)
library(Hmisc)

#### Load data ####

traits <- readRDS("Data/00_Traits_alien_iast_birds")


#### Clean trait data ####

data_ft <- traits %>%
  select(binomial, Hand.Wing.Index:insul_level) %>%
  distinct() %>%
  column_to_rownames("binomial") %>%
  filter(!is.na(hab_sum)) %>%
  mutate(hab_sum=as.character(hab_sum)) %>%
  mutate(hab_sum=if_else(!(hab_sum %in% c("1","2","3","4")),"5+",hab_sum)) %>%
  mutate(Habitat=factor(hab_sum, ordered = T,
                        levels = c("1","2","3","4","5+"))) %>%
  select(-c(hab_sum, volant, insular_endemic))


str(data_ft)


#### Run PCoA on traits ####

# 1. compute the dissimilarity matrix between species 
# (Gower distance for non numeric var) & compute PCoA on dis matrix
# 2. use a pca function which can include mixed data

# Mouillot et al 2020
mat_dis <- cluster::daisy(data_ft, metric = "gower")
mat_pcoa <- ape::pcoa(mat_dis)

nbdim = 10
nbdim<-min(nbdim,ncol(mat_pcoa$vectors) )
# keeping species coordinates on the 'nbdim' axes
mat_coord<-mat_pcoa$vectors[,1:nbdim]
row.names(mat_coord)<-row.names(data_ft)
colnames(mat_coord)<-paste("PC",1:nbdim,sep="")

sum(mat_pcoa$values$Eigenvalues)
mat_pcoa$values$Eigenvalues[1]/sum(mat_pcoa$values$Eigenvalues)
mat_pcoa$values$Eigenvalues[2]/sum(mat_pcoa$values$Eigenvalues)
mat_pcoa$values$Eigenvalues[3]/sum(mat_pcoa$values$Eigenvalues)
mat_pcoa$values$Eigenvalues[4]/sum(mat_pcoa$values$Eigenvalues)



#### Add trait correlation with ecological space ####

# select coordinates from pcoa
coordinates = as.data.frame(mat_coord) 
colnames(coordinates) <- paste0("V", 1:ncol(coordinates))

# Shape traits as numeric or binary variables (with dummy)
# to compute correlation with species
str(data_ft)

data_ft_num <- data_ft %>%
  mutate(Habitat = as.numeric(Habitat),
         value = 1,
         Trophic.Level = paste0("Diet.", Trophic.Level)) %>% 
  pivot_wider(names_from = Trophic.Level, values_from = value, values_fill = 0) %>%
  mutate(value = 1,
         Primary.Lifestyle = paste0("ForN.", Primary.Lifestyle)) %>% 
  pivot_wider(names_from = Primary.Lifestyle, values_from = value, values_fill = 0) %>%
  mutate(value = 1,
         bioreg = paste0("Bioreg.", bioreg)) %>% 
  pivot_wider(names_from = bioreg, values_from = value, values_fill = 0)

colnames(data_ft_num)

# change to matrix
data_ft_num_mat <- as.matrix(data_ft_num)

# Compute correlations between species coordinates and trait values
cor_mat<-rcorr(as.matrix(cbind(data_ft_num_mat, 
                               coordinates %>% dplyr::select(V1:V4))), 
               type = 'spearman')
# Matrix of correlation coeficients
R <- as.data.frame(round(cor_mat$r, 2)) %>% 
  select(V1:V4) %>%
  rownames_to_column("Trait") %>%
  filter(!(Trait %in% c("V1", "V2", "V3", "V4")))
# Matrix of p-value 
p <- as.data.frame(round(cor_mat$P,3))  %>% 
  select(V1:V4) %>%
  rownames_to_column("Trait") %>%
  filter(!(Trait %in% c("V1", "V2", "V3", "V4")))

colnames(R) <- c("Trait", paste0("r_PC", 1:4))
colnames(p) <- c("Trait", paste0("p_PC", 1:4))

cor_tab <- left_join(R, p, by="Trait")
str(cor_tab)



##### save PCA data & correlation traits axes ####

pc_to_save <- left_join(
  traits %>% select(binomial:Species1),
  bind_cols(data_ft, coordinates) %>% rownames_to_column("binomial"))

saveRDS(pc_to_save, "Output/01_PCoA_traits_groups")

class(cor_tab)
write.table(x = cor_tab, 
            file = "Output/01_Corr_traits_axis_pcoa.csv",
            sep=";", row.names = F)
