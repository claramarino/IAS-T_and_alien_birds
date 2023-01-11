rm(list=ls())

##### Trait by trait analyses #####

library(tidyverse)
library(ggpubr)
library(car)
library(emmeans)

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
  mutate_if(is.character, as.factor)

# clean group data 
tg <- ei_no_dd %>% select(binomial, group2, Hand.Wing.Index:Habitat)

# load trait data
all_birds <- readRDS("Data/00_Traits_worldwide_birds") %>%
  mutate_if(is.character, as.factor) %>%
  # convert Mass, Beak & Clutch with log 
  mutate(ln.Mass = log(Mass),
         ln.Clutch = log(Clutch),
         ln.Beak.Depth = log(Beak.Depth),
         ln.Beak.Length_Nares = log(Beak.Length_Nares)) %>%
  # remove converted var + trophic niche
  select(-c(Mass, Clutch, Beak.Depth, Beak.Length_Nares, Trophic.Niche)) %>%
  mutate(group2 = "All_birds") %>%
  filter(!is.na(hab_sum)) %>%
  mutate(hab_sum=as.character(hab_sum)) %>%
  mutate(hab_sum=if_else(!(hab_sum %in% c("1","2","3","4")),"5+",hab_sum)) %>%
  mutate(Habitat=factor(hab_sum, ordered = T,
                        levels = c("1","2","3","4","5+"))) %>%
  select(-c(hab_sum))


tg_all <- bind_rows(tg, all_birds)
length(tg_all$binomial[duplicated(tg_all$binomial)])

group_order <- c("All_birds", "IAS-T_imp", "IAS-T_no_imp",
                 "EICAT_no_imp", "EICAT_imp")
tg_all <- tg_all %>%
  mutate(group2 = as.factor(group2)) %>%
  mutate(group2 = factor(group2, level = group_order))

str(tg_all)

#####################################################################

# Impact versus no impact

#---------------------- numeric var

# all morphological measures are corrected by ln.Mass

# HWI corrected by bodymass
a <- ggplot(tg_all, aes(x= group2, y= log(Hand.Wing.Index/ln.Mass + 1))) +
  geom_boxplot(alpha = 0.5)+
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) +
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1)+ 
  #xlab("") +
  theme_classic2()

# Tail Length corrected by mass
b <- ggplot(tg_all, aes(group2, Tail.Length/ln.Mass)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1) +
  xlab("") +
  theme_classic2()


# Beak Depth corrected by bodymass
c <- ggplot(tg_all, aes(group2, ln.Beak.Depth/ln.Mass)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1)+ 
  xlab("") +
  theme_classic2()

# Beak Length corrected by mass
d <- ggplot(tg_all, aes(group2, log(ln.Beak.Length_Nares/ln.Mass+1))) +
  geom_boxplot(alpha = 0.5) +
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1)+ 
  xlab("") +
  theme_classic2()


# Mass 
e <- ggplot(tg_all, aes(group2, ln.Mass)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1)+ 
  xlab("") +
  theme_classic2()

# Clutch 
f <- ggplot(tg_all, aes(group2, ln.Clutch)) +
  geom_boxplot(alpha = 0.5) +
  geom_point(data=tg_all %>% filter(group2 != "All_birds"),
             position = position_jitter(width = .2), 
             col = "darkcyan", alpha = 0.4) + 
  geom_boxplot(alpha = 0.5) +
  stat_summary(fun = "mean", color="violetred3", shape=18, size = 1)+ 
  xlab("") +
  theme_classic2()


pdf("Output/03_Fig3_Single_traits_num_impact.pdf", 10, 6)
ggarrange(a, b, e, c, d, f, ncol=3, nrow = 2 )
dev.off()


# statistical test with all birds

ls_all <- data.frame()

for(i in 1:100){
  
  tg_all_rand <- bind_rows(
    tg_all %>%
      filter(group2 != "All_birds"),
    tg_all %>%
      filter(group2 == "All_birds" & 
               binomial %in% sample(levels(tg_all$binomial), 150)))
  
  # hwi / ln.mass
  ls_all <- bind_rows(ls_all,
                      as.data.frame(lsmeans(lm(Hand.Wing.Index/ln.Mass ~ group2, data=tg_all_rand), 
                                            pairwise ~ group2, adjust = "bonferroni")$contrasts) %>%
                        mutate(Trait = "hwi_by_mass"),
                      as.data.frame(lsmeans(lm(Tail.Length/ln.Mass ~ group2, data=tg_all_rand), 
                                            pairwise ~ group2, adjust = "bonferroni")$contrasts) %>%
                        mutate(Trait = "tail_by_mass"),
                      as.data.frame(lsmeans(lm(ln.Beak.Depth/ln.Mass ~ group2, data=tg_all_rand), 
                                            pairwise ~ group2, adjust = "bonferroni")$contrasts) %>%
                        mutate(Trait = "beak_depth_by_mass"),
                      as.data.frame(lsmeans(lm(ln.Beak.Length_Nares/ln.Mass ~ group2, data=tg_all_rand), 
                                            pairwise ~ group2, adjust = "bonferroni")$contrasts) %>%
                        mutate(Trait = "beak_length_by_mass"), 
                      as.data.frame(lsmeans(lm(ln.Mass ~ group2, data=tg_all_rand), 
                                            pairwise ~ group2, adjust = "bonferroni")$contrasts) %>%
                        mutate(Trait = "mass"), 
                      as.data.frame(lsmeans(lm(ln.Clutch ~ group2, data=tg_all_rand), 
                                            pairwise ~ group2, adjust = "bonferroni")$contrasts) %>%
                        mutate(Trait = "clutch"))
  
  
  tg_all_ias <- tg_all_rand %>% 
    filter(group2 %in% c("IAS-T_imp", "All_birds", "IAS-T_no_imp")) %>%
    mutate_if(is.factor, as.character)
  
  
  tg_all_eicat <- tg_all_rand %>% 
    filter(group2 %in% c( "EICAT_imp", "EICAT_no_imp")) %>%
    mutate_if(is.factor, as.character)
  
  tbl <- table(tg_all_ias$group2, tg_all_ias$Habitat)
  chisq.test(tbl)
  
  tbl <- table(tg_all_eicat$group2, tg_all_eicat$Habitat)
  tbl
  chisq.test(tbl)
  
  a$p.value  
  
  
  
}

str(ls_all)

ls_summar <- ls_all %>% 
  group_by(contrast, Trait) %>%
  summarise(mean_est = mean(estimate),
            mean_pval = mean(p.value)) %>%
  mutate(mean_pval_round = round(mean_pval, 3)) %>%
  filter(mean_pval_round < 0.1)

write.table(ls_summar, "Output/03_ls_means_all_birds.csv", sep=";", row.names = F)



#---------------------- Categorial var


# trophic level
a1 <- ggplot(tg_all, aes(group2, fill=Trophic.Level)) +
  geom_bar(aes(y=..count..), position = "fill") +
  xlab("") +ylab("")+
  theme_classic()+
  theme(legend.position = "top")+
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  scale_color_viridis_d(option = "viridis", direction = -1)

chisq.test(tg_all$bioreg, tg_all$Trophic.Level)

# habitat
b1 <- ggplot(tg_all, 
             aes(group2, fill=Habitat)) +
  geom_bar(aes(y=..count..), position = "fill") +
  xlab("") +ylab("")+
  theme_classic()+
  theme(legend.position = "top")+
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  scale_color_viridis_d(option = "viridis", direction = -1)

chisq.test(tg_all$Habitat, tg_all$group2)

# primary lifestyle
c1 <- ggplot(tg_all,
             aes(group2, fill=Primary.Lifestyle)) +
  geom_bar(aes(y=..count..), position = "fill") +
  xlab("") +ylab("")+
  theme_classic()+
  theme(legend.position = "top")+
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  scale_color_viridis_d(option = "viridis", direction = -1)

chisq.test(tg_all$Primary.Lifestyle, tg_all$group2)


# insular level
d1 <- ggplot(tg_all %>% mutate(insul_level = as.factor(insul_level)),
             aes(group2, fill=insul_level)) +
  geom_bar(aes(y=..count..), position = "fill") +
  xlab("") +ylab("")+
  theme_classic()+
  theme(legend.position = "top")+
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  scale_color_viridis_d(option = "viridis", direction = -1)

chisq.test(tg_all$insul_level, tg_all$group2)


# bioregion
e1 <- ggplot(tg_all %>% filter(!is.na(bioreg)),
             aes(group2, fill=bioreg)) +
  geom_bar(aes(y=..count..), position = "fill") +
  xlab("") +ylab("")+
  theme_classic()+
  theme(legend.position = "top")+
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  scale_color_viridis_d(option = "viridis", direction = -1)

chisq.test(tg_all$bioreg, tg_all$group2)

pdf("Output/03_Fig2_Single_traits_cate_impact.pdf", 9, 7)
ggarrange(a1, b1, c1, d1, e1, ncol=3, nrow = 2)
dev.off()

