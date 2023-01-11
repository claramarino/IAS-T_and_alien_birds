# Evaluation of status prediction for potential DD species

rm(list=ls())

library(tidyverse)
library(cluster)
options(dplyr.summarise.inform = FALSE)

#### Load data ####

tg <- readRDS("Data/00_Traits_alien_iast_birds")
all_groups <- tg %>% filter(!is.na(hab_sum))

#### Set up species groups ####

ei_no_dd <- all_groups %>% 
  filter(combi!="0_0_0") %>%
  mutate_if(is.character, as.factor)

#### Compute distance matrix ####

eicat_iast_ft <- tg %>%
  select(Species1, Hand.Wing.Index:insul_level) %>%
  distinct() %>%
  column_to_rownames("Species1") %>%
  filter(!is.na(hab_sum)) %>%
  mutate(hab_sum=as.character(hab_sum)) %>%
  mutate(hab_sum=if_else(!(hab_sum %in% c("1","2","3","4")),"5+",hab_sum)) %>%
  mutate(Habitat=factor(hab_sum, ordered = T,
                        levels = c("1","2","3","4","5+"))) %>%
  select(-c(hab_sum, volant, insular_endemic))


# distance matrix based on traits 
dis.gower <- daisy(eicat_iast_ft,metric='gower')
summary(dis.gower)
matrix = as.matrix(dis.gower)
dim(matrix)
matrix[1,2]

# select only no_DD species in rows & in columns
# we will evaluate the prediction using only species with a known status

no_dd <- unique(as.character(pull(ei_no_dd, Species1)))

df_dist <- as.data.frame(matrix)  %>%
  select(all_of(no_dd)) %>%
  rownames_to_column("no_DD_sp") %>%
  filter(no_DD_sp %in% no_dd) %>%
  column_to_rownames("no_DD_sp")

# create species to group correspondence
corresp <- ei_no_dd %>% distinct(Species1, group2) %>%
  rename(group = group2) 

#### Method 'Out of the bag' ####

# fill the  table for each of the 574 species (to remove in the predictors)

#------------- Initialize tables for metrics

# for min and mean distance
long_df_all_sp <- data.frame(
  no_DD_sp = character(0), group = character(0),
  mean_dist = numeric(0), min_dist = numeric(0))
# 10 closest neighbors
k_closest <- k_closest_sp <- data.frame(matrix(ncol = 10, nrow = 0)) %>%
  mutate_all(as.character)
colnames(k_closest) <- colnames(k_closest_sp) <- paste("closest_", 1:10, sep="")
# species in buffers 
d_buffer <- c("d1" = 0.21, 
              "d2" = mean(as.matrix(df_dist))) 
buffer_df_all <- data.frame(
  group = character(0), mean_buffer = numeric(0), n_b = numeric(0),
  prop_b = numeric(0), no_DD_sp = character(0),d_buffer=character(0))
#source function for calculating sp in buffer, mean dist & prop in buffer
source("R/find_sp_in_buffer.R")

#------------ Loop on all species

for (sp in no_dd){
  
  # test iteration
  #sp = no_dd[6]
  
  # select focus species & remove it from the predictors
  df_dist_out <- df_dist %>%
    select(-all_of(sp))%>%
    rownames_to_column("no_DD_sp") %>%
    filter(no_DD_sp == sp) %>%
    column_to_rownames("no_DD_sp")
  
  #___________________ Min & mean distance with groups
  
  long_df <- df_dist_out %>% rownames_to_column("no_DD_sp") %>%
    pivot_longer(!no_DD_sp, names_to = "sp2", values_to = "dist")
  
  long_df_sp <- left_join(long_df, corresp %>% rename(sp2 = Species1), 
                          by = "sp2") %>%
    group_by(no_DD_sp, group) %>%
    summarise(mean_dist = mean(dist), min_dist = min(dist))
  
  # implement in final long df with all no_dd species
  long_df_all_sp <- bind_rows(long_df_all_sp, long_df_sp)
  
  #____________________ 10 closest neighbours
  
  col1 <- as.data.frame(t(df_dist_out)) %>% rownames_to_column("sp")
  sorted <- col1[order(col1[,2]),]
  for (k in 1:10){
      k_closest_sp[1,k]<- sorted$sp[k]
  } 
  row.names(k_closest_sp) <- sp
  
  # implement in final k_closest dataset
  k_closest <- bind_rows(k_closest, k_closest_sp)
  
  #____________________ species in buffers
  
  buffer_d1 <- find_sp_in_buffer_role(df_dist_out, sp, d_buffer[1])
  buffer_d2 <- find_sp_in_buffer_role(df_dist_out, sp, d_buffer[2])
  
  buffer_d1[[1]] <- buffer_d1[[1]] %>% 
    mutate(no_DD_sp = sp, d_buffer = "d1")
  buffer_d2[[1]] <- buffer_d2[[1]] %>% 
    mutate(no_DD_sp = sp, d_buffer = "d2")
  
  buffer_df_all <- bind_rows(buffer_df_all, 
                             do.call(rbind.data.frame, buffer_d1),
                             do.call(rbind.data.frame, buffer_d2))
  
}

#___________________ Final metrics dataset

k_closest_lg <- k_closest %>%
  rownames_to_column("no_DD_sp") %>%
  pivot_longer(!no_DD_sp, names_to = "rank", values_to = "sp_2")

k_closest_gps <- left_join(k_closest_lg, corresp %>%
                             rename(sp_2 = Species1), 
                           by = "sp_2") %>% 
  group_by(no_DD_sp, group) %>% summarise(n_10_closest = n())

buffer_d1_df <- buffer_df_all %>%
  filter(d_buffer == "d1") %>%
  rename(mean_b1 = mean_buffer,
         n_b1 = n_b,
         prop_b1 = prop_b) %>%
  select(-d_buffer)
  
buffer_d2_df <- buffer_df_all %>%
  filter(d_buffer == "d2") %>%
  rename(mean_b2 = mean_buffer,
         n_b2 = n_b,
         prop_b2 = prop_b) %>%
  select(-d_buffer)
  
evaluate_no_dd <- full_join(
  # Add k closest neighbours
  full_join(long_df_all_sp, k_closest_gps, by=c("no_DD_sp","group")),
  # Add buffer metrics
  # join with d2 first because no NA
  full_join(buffer_d2_df, buffer_d1_df, by=c("no_DD_sp","group"))) %>%
  # replace NA by 0 for proportions (but nor for means)
  mutate_at(c("n_10_closest","n_b1","prop_b1", "n_b2","prop_b2"), 
            ~ ifelse(is.na(.), 0, .))

saveRDS(evaluate_no_dd, "Output/04_Metrics_evaluate_no_DD_role")


#### Attribute a status for each metric ####

# for min_dist & mean_dist, attribute status of minimal distance group
# idem for mean_b1 & mean_b2: minimal mean distance in each buffer

# for proportion of closest species:
# take group status for which the proportion is maximal when removing null_prop

# load metrics file
evaluate_no_dd <- readRDS("Output/04_Metrics_evaluate_no_DD_role") %>%
  mutate(group = as.character(group))

# calculate null proportions of group representation
null_prop <- c(table(as.character(ei_no_dd$group2))/nrow(ei_no_dd))

evaluate_no_dd <- left_join(evaluate_no_dd,
                            as.data.frame(null_prop) %>% 
                              rownames_to_column("group"),
                            by = "group")

evaluate_no_dd <- evaluate_no_dd %>%
  mutate(prop_10_closest = n_10_closest/10) %>%
  mutate(prop_b1_corr = prop_b1 - null_prop,
         prop_b2_corr = prop_b2 - null_prop,
         prop_10_closest_corr = prop_10_closest - null_prop) %>%
  mutate(prop_b1_corr = if_else(is.na(mean_b1), mean_b1, prop_b1_corr),
         prop_b2_corr = if_else(is.na(mean_b2), mean_b2, prop_b2_corr))

# initialize empty df to fill with predicted status
status_all_meth <- data.frame(matrix(ncol = 8, nrow = length(no_dd)))
names(status_all_meth) <- c("no_DD_sp","min_dist","mean_dist","mean_b1", 
                            "prop_b1_corr", "mean_b2", "prop_b2_corr",
                            "prop_10_closest_corr")
status_all_meth$no_DD_sp <- no_dd

# attribute status for each dd sp & each method
for(sp in no_dd){
  
  evaluate_no_dd_sp <- evaluate_no_dd %>% filter(no_DD_sp == sp)
  
  # vars for which to take minimal value
  vars_min <- c("min_dist","mean_dist","mean_b1","mean_b2")
  
  for(col in vars_min){
    # predict only if at least 1 group is rpzted in the metric
    if(sum(is.na(pull(evaluate_no_dd_sp,col)))!=3){
      status_all_meth[which(status_all_meth$no_DD_sp==sp),col] <-
        evaluate_no_dd_sp$group[which.min(pull(evaluate_no_dd_sp,col))]
    } else {
      status_all_meth[which(status_all_meth$no_DD_sp==sp),col] <- NA
    }
  }
  # vars for which to take maximal value after correction by null_prop
  vars_max <- c("prop_b1_corr", "prop_b2_corr", "prop_10_closest_corr")
  for(col2 in vars_max){
    if(sum(is.na(pull(evaluate_no_dd_sp,col2)))!=3){
      status_all_meth[which(status_all_meth$no_DD_sp==sp),col2] <-
        evaluate_no_dd_sp$group[which.max(pull(evaluate_no_dd_sp,col2))]
    } else {
      status_all_meth[which(status_all_meth$no_DD_sp==sp),col2] <- NA
    }
  }
  
}


##### Evaluate which metrics are the best predictors #####

compare <- left_join(corresp, 
                     status_all_meth %>% rename(Species1 = no_DD_sp))%>%
  mutate_if(is.factor, as.character)

# Make a contingency table per metric
metric_names <- c("min_dist","mean_dist","mean_b1", "prop_b1_corr",
                  "mean_b2","prop_b2_corr", "prop_10_closest_corr")
ct_list <- vector(mode="list", length(metric_names))
names(ct_list) <- metric_names
  
for(col in metric_names){
  ct_list[[col]] <- table(compare$group, pull(compare, col))
}


# percent of well predicted values (matrix diagonal)
lapply(ct_list, function(x){
  sum(diag(as.matrix(x)))/length(no_dd)
})

# Focus on EICAT_imp --> TP, TN, FP, FN => specificity, sensibility
calculate_eicat_imp_success <- function(mat){
  TP <- mat[1,1]
  TN <- sum(mat[2:4,2:4])
  FP <- sum(mat[2:4,1])
  FN <- sum(mat[1,2:4])
  index <- c("TP" = TP, "TN" = TN, "FP" = FP, "FN" = FN,
             "sensitivity" = TP/(TP+FN), "specificity" = TN/(TN+FP))
  return(index)
}

lapply(ct_list, function(x){
  mat <- as.matrix(x)
  index <- calculate_eicat_imp_success(mat)
  return(index)
  })

