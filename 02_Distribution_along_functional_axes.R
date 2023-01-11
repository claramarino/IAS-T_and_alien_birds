rm(list=ls())

########## Analyses in functional space ##########

source("R/R_rainclouds.R")
library(ggpubr)
library(tidyverse)


# load data 

ei_no_dd <- readRDS("Output/01_PCoA_traits_groups") %>% 
  filter(!is.na(V1)) %>%
  filter(combi!="0_0_0") %>%
  mutate(combi2 = case_when(
    combi=="0_0_1" ~ "Only indirect",
    combi=="0_1_1" ~ "Direct_Indirect",
    combi=="1_1_1" ~ "All mecha",
    combi=="0_1_0" ~ "Only direct",
    combi=="1_0_0" ~ "Only ecosystem",
    combi=="1_1_0" ~ "Ecosystem_Direct",
    combi=="1_0_1" ~ "Ecosystem_Indirect")) %>%
  mutate(group2 = case_when(
    group2=="EICAT_imp" ~ "Alien birds with impact",
    group2=="EICAT_no_imp" ~ "Alien birds without impact",
    group2=="IAS-T_imp" ~ "Native birds impacted by IAS",
    group2=="IAS-T_no_imp" ~ "Native birds not impacted by IAS"
  )) %>%
  mutate_if(is.character, as.factor)

# remove all groups with less than ndim+1 species ?
table(ei_no_dd$combi, ei_no_dd$group)
ei_no_dd_m <- ei_no_dd %>% filter(!(group=="EICAT" & combi =="1_1_0")) %>%
  filter(!(group=="IAS-T" & combi =="1_0_1")) %>%
  filter(combi!="1_1_1") # remove cate all because not usefull
table(ei_no_dd_m$combi, ei_no_dd_m$group)

group_color = c("Direct_Indirect" = "darkorange1", "Ecosystem_Direct" = 
                  "darkorchid3", "Ecosystem_Indirect" = "chartreuse3",
                "Only direct" = "firebrick2","Only ecosystem" = "dodgerblue3",
                "Only indirect" = "gold")

# keep only alone mecha
ei_no_dd_alone <- ei_no_dd_m %>%
  filter(combi2 %in% c("Only direct", "Only ecosystem", "Only indirect"))
table(ei_no_dd_alone$combi2, ei_no_dd_alone$group)



#### Fig 1 -- Impact groups along PCA axis ####

# IAS-T-imp, IAS-T-no-imp, EICAT-imp, EICAT-no-imp

p1i <- ggplot(ei_no_dd, aes(x= 1, y = V1, fill = group2)) +
  geom_flat_violin(aes(fill = group2),
                   position = position_nudge(x = .1, y = 0), 
                   adjust = 1, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = 1-.15, y = V1, colour = group2),
             position = position_jitter(width = .05), size = .25, shape = 20)+
  geom_boxplot(aes(x=1, y = V1, fill = group2),
               outlier.shape = NA, alpha = .5, width = .2, colour = "grey30")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab("PC1") +
  labs(fill="", colour="") +
  coord_flip()

p2i <- ggplot(ei_no_dd, aes(x= 1, y = V2, fill = group2)) +
  geom_flat_violin(aes(fill = group2),
                   position = position_nudge(x = .1, y = 0), 
                   adjust = 1, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = 1-.15, y = V2, colour = group2),
             position = position_jitter(width = .05), size = .25, shape = 20)+
  geom_boxplot(aes(x=1, y = V2, fill = group2),
               outlier.shape = NA, alpha = .5, width = .2, colour = "grey30")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("PC2") +
  labs(fill="", colour="") +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip()

p3i <- ggplot(ei_no_dd, aes(x= 1, y = V3, fill = group2)) +
  geom_flat_violin(aes(fill = group2),
                   position = position_nudge(x = .1, y = 0), 
                   adjust = 1, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = 1-.15, y = V3, colour = group2),
             position = position_jitter(width = .05), size = .25, shape = 20)+
  geom_boxplot(aes(x=1, y = V3, fill = group2),
               outlier.shape = NA, alpha = .5, width = .2, colour = "grey30")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("PC3") +
  labs(fill="", colour="") +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip()

p4i <- ggplot(ei_no_dd, aes(x= 1, y = V4, fill = group2)) +
  geom_flat_violin(aes(fill = group2),
                   position = position_nudge(x = .1, y = 0), 
                   adjust = 1, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = 1-.15, y = V4, colour = group2),
             position = position_jitter(width = .05), size = .25, shape = 20)+
  geom_boxplot(aes(x=1, y = V4, fill = group2),
               outlier.shape = NA, alpha = .5, width = .2, colour = "grey30")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ylab("PC4") +
  labs(fill="", colour="") +
  theme_classic()+
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), 
        axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_flip()

a = ggarrange(p1i, p2i, p3i, ncol = 3, legend = "top", common.legend = T)
a

# save plot
pdf(file = "Output/02_Fig1_impact_groups_pcoa.pdf", 8, 4)
print(a)
dev.off()



#### Statistical test for comparing distributions ####
# add mean comparison along PC axes

axis = c("V1", "V2", "V3", "V4")
duo_imp = combn(c("Alien birds with impact", "Alien birds without impact",
                  "Native birds impacted by IAS", "Native birds not impacted by IAS"), 2)

ks_wt_df_imp <- data.frame(
  Axis = rep(axis, each = 6),
  a = character(24),
  b = character(24))

for (i in 1:ncol(duo_imp)){
  ks_wt_df_imp[seq(i,24, by = 6),"a"] <- duo_imp[1,i]
  ks_wt_df_imp[seq(i,24, by = 6),"b"] <- duo_imp[2,i]
}

for(i in 1:nrow(ks_wt_df_imp)){
  gpa = pull(ei_no_dd %>% 
               filter(group2 == ks_wt_df_imp$a[i]), 
             ks_wt_df_imp$Axis[i])
  gpb = pull(ei_no_dd %>% 
               filter(group2 == ks_wt_df_imp$b[i]), 
             ks_wt_df_imp$Axis[i])
  k = ks.test(gpa, gpb)
  w = wilcox.test(gpa, gpb)
  ks_wt_df_imp$KS_D[i] = k$statistic
  ks_wt_df_imp$KS_p_val[i] = k$p.value
  ks_wt_df_imp$WT_w[i] = w$statistic
  ks_wt_df_imp$WT_p_val[i] = w$p.value
  
}

ks_wt_df_imp <- ks_wt_df_imp %>%
  mutate(KS_p_val_round = round(ks_wt_df_imp$KS_p_val, 4),
         WT_p_val_round = round(ks_wt_df_imp$WT_p_val, 4))

write.table(ks_wt_df_imp, "Output/02_KS_WT_test_impact.csv", sep = ";", row.names = F)




#### Graphs - Figure 4 ####

p1 <- ggplot(ei_no_dd_alone, aes(x = group, y = V1, fill = combi2)) +
  geom_flat_violin(aes(fill = combi2),
                   position = position_nudge(x = .15, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_boxplot(aes(x = group, y = V1, fill = combi2),
               outlier.shape = NA, alpha = .8, width = .3, colour = "grey30")+
  scale_fill_manual(values = group_color)+
  ylab("PC1") +
  labs(fill="", colour="") +
  theme_classic()+
  coord_flip()


p2 <- ggplot(ei_no_dd_alone, aes(x = group, y = V2, fill = combi2)) +
  geom_flat_violin(aes(fill = combi2),
                   position = position_nudge(x = .15, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_boxplot(aes(x = group, y = V2, fill = combi2),
               outlier.shape = NA, alpha = .8, width = .3, colour = "grey30")+
  scale_fill_manual(values = group_color)+ 
  ylab("PC2") +
  labs(fill="", colour="") +
  theme_classic()+
  coord_flip()

p3 <- ggplot(ei_no_dd_alone, aes(x = group, y = V3, fill = combi2)) +
  geom_flat_violin(position = position_nudge(x = 0.15, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_boxplot(aes(x = group, y = V3, fill = combi2),
               outlier.shape = NA, alpha = .8, width = .3, colour = "grey30")+
  scale_fill_manual(values = group_color)+
  ylab("PC3") +
  labs(fill="") +
  theme_classic()+
  coord_flip()

p4 <- ggplot(ei_no_dd_alone, aes(x = group, y = V4, fill = combi2)) +
  geom_flat_violin(position = position_nudge(x = 0.15, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_boxplot(aes(x = group, y = V4, fill = combi2),
               outlier.shape = NA, alpha = .8, width = .3, colour = "grey30")+
  scale_fill_manual(values = group_color)+
  ylab("PC4") +
  labs(fill="") +
  theme_classic()+
  coord_flip()

b = ggarrange(p1, p2, p3, ncol = 3, legend = "top", common.legend = T)
b

# save plot
pdf(file = "Output/02_Fig4_mecha_groups_pcoa.pdf", 10, 6)
print(b)
dev.off()

#### Statistical tests ####

axis = c("V1", "V2", "V3", "V4")
groups = c("EICAT","IAS-T")
duo = combn(c("1_0_0", "0_1_0", "0_0_1"), 2)

ks_wt_df <- data.frame(
  Axis = rep(axis, each = 6),
  Group = rep(groups, times = 4, each = 3),
  a = character(24),
  b = character(24))

for (i in 1:ncol(duo)){
  ks_wt_df[seq(i,24, by = 3),"a"] <- duo[1,i]
  ks_wt_df[seq(i,24, by = 3),"b"] <- duo[2,i]
}

for(i in 1:nrow(ks_wt_df)){
  gpa = pull(ei_no_dd %>% 
               filter(group==ks_wt_df$Group[i] & combi==ks_wt_df$a[i]), 
             ks_wt_df$Axis[i])
  gpb = pull(ei_no_dd %>% 
               filter(group==ks_wt_df$Group[i] & combi==ks_wt_df$b[i]), 
             ks_wt_df$Axis[i])
  k = ks.test(gpa, gpb)
  w = wilcox.test(gpa, gpb)
  ks_wt_df$KS_D[i] = k$statistic
  ks_wt_df$KS_p_val[i] = k$p.value
  ks_wt_df$WT_w[i] = w$statistic
  ks_wt_df$WT_p_val[i] = w$p.value
}


ks_wt_df <- ks_wt_df %>%
  mutate(KS_p_val_round = round(ks_wt_df$KS_p_val, 4),
         WT_p_val_round = round(ks_wt_df$WT_p_val, 4))

write.table(ks_wt_df, "Output/02_KS_test_mecha_alone.csv", sep = ";", row.names = F)
