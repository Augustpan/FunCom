library(tidyverse)
library(vegan)
library(lmerTest)
library(ggpubr)
library(cowplot)
library(VennDiagram)

load("../data/full_data_set.RData")

guild = read_tsv("../data/taxonomy.guilds.tsv")
guild = select(
  guild,
  otu_name = "OTU ID",
  trophic_mode="Trophic Mode", 
  guild="Guild", 
  growth_morph="Growth Morphology"
)

pack_leaf = leaf %>% 
  left_join(smd, by="SMD_pid") %>% 
  na.omit()

otu_leaf_nat = filter(pack_leaf, SMD_origin=="nat") %>% select(starts_with("OTU"))
otu_leaf_int = filter(pack_leaf, SMD_origin=="int") %>% select(starts_with("OTU"))

patho_otu = filter(guild, trophic_mode=="Pathotroph")
patho_otu_nat = otu_leaf_nat[,colSums(otu_leaf_nat)>0] %>% select(any_of(patho_otu$otu_name)) %>% colnames()
patho_otu_int = otu_leaf_int[,colSums(otu_leaf_int)>0] %>% select(any_of(patho_otu$otu_name)) %>% colnames()
venn.diagram(list(nat=symbio_otu_nat, int=symbio_otu_int), filename="../result/patho_venn.jpg")

symbio_otu = filter(guild, trophic_mode=="Symbiotroph")
symbio_otu_nat = otu_leaf_nat[,colSums(otu_leaf_nat)>0] %>% select(any_of(symbio_otu$otu_name)) %>% colnames()
symbio_otu_int = otu_leaf_int[,colSums(otu_leaf_int)>0] %>% select(any_of(symbio_otu$otu_name)) %>% colnames()
venn.diagram(list(nat=symbio_otu_nat, int=symbio_otu_int), filename="../result/symbio_venn.jpg")

sapro_otu = filter(guild, trophic_mode=="Saprotroph")
sapro_otu_nat = otu_leaf_nat[,colSums(otu_leaf_nat)>0] %>% select(any_of(sapro_otu$otu_name)) %>% colnames()
sapro_otu_int = otu_leaf_int[,colSums(otu_leaf_int)>0] %>% select(any_of(sapro_otu$otu_name)) %>% colnames()
venn.diagram(list(nat=sapro_otu_nat, int=sapro_otu_int), filename="../result/sapro_venn.jpg")

pack_t_leaf = t_leaf %>% 
  left_join(smd, by="SMD_pid") %>% 
  na.omit()
tro_leaf = select(pack_t_leaf, starts_with("TM"))

fac_comp = as.factor(pack_t_leaf$compartment)
fac_comp = ordered(fac_comp, c(2,3,1,4,6,7,5))
fac = as.factor(pack_t_leaf$compartment)
pairwise.t.test(pack_t_leaf$TM_patho, pack_t_leaf$compartment)

ggplot(aes(x=SMD_lat, y=TM_patho), data=pack_t_leaf) +
  geom_point(aes(color=fac)) + 
  geom_smooth() + 
  facet_wrap(~SMD_origin)

# BOX PLOT
g1 = ggplot(aes(x=fac_comp, y=TM_patho), data=pack_t_leaf) +
  geom_boxplot()
g2 = ggplot(aes(x=SMD_origin, y=TM_patho), data=pack_t_leaf) +
  geom_boxplot() +
  stat_compare_means(method="t.test",comparisons=list(c("int","nat"))) +
  ylab("")
plot_grid(g1,g2,rel_widths=c(2.5,1))
ggsave("../result/patho_comp.jpg")

g1 = ggplot(aes(x=fac_comp, y=TM_sapro), data=pack_t_leaf) +
  geom_boxplot()
g2 = ggplot(aes(x=SMD_origin, y=TM_sapro), data=pack_t_leaf) +
  geom_boxplot() +
  stat_compare_means(method="t.test",comparisons=list(c("int","nat"))) +
  ylab("")
plot_grid(g1,g2,rel_widths=c(2.5,1))
ggsave("../result/sapro_comp.jpg")

g1 = ggplot(aes(x=fac_comp, y=TM_symbio), data=pack_t_leaf) +
  geom_boxplot()
g2 = ggplot(aes(x=SMD_origin, y=TM_symbio), data=pack_t_leaf) +
  geom_boxplot() +
  stat_compare_means(method="t.test",comparisons=list(c("int","nat"))) +
  ylab("")
plot_grid(g1,g2,rel_widths=c(2.5,1))
ggsave("../result/symbio_comp.jpg")

g1 = ggplot(aes(x=fac_comp, y=TM_patho_sapro), data=pack_t_leaf) +
  geom_boxplot()
g2 = ggplot(aes(x=SMD_origin, y=TM_patho_sapro), data=pack_t_leaf) +
  geom_boxplot() +
  stat_compare_means(method="t.test",comparisons=list(c("int","nat"))) +
  ylab("")
plot_grid(g1,g2,rel_widths=c(2.5,1))
ggsave("../result/patho_sapro_comp.jpg")

g1 = ggplot(aes(x=fac_comp, y=TM_patho_symbio), data=pack_t_leaf) +
  geom_boxplot()
g2 = ggplot(aes(x=SMD_origin, y=TM_patho_symbio), data=pack_t_leaf) +
  geom_boxplot() +
  stat_compare_means(method="t.test",comparisons=list(c("int","nat"))) +
  ylab("")
plot_grid(g1,g2,rel_widths=c(2.5,1))
ggsave("../result/patho_symbio_comp.jpg")

g1 = ggplot(aes(x=fac_comp, y=TM_sapro_symbio), data=pack_t_leaf) +
  geom_boxplot()
g2 = ggplot(aes(x=SMD_origin, y=TM_sapro_symbio), data=pack_t_leaf) +
  geom_boxplot() +
  stat_compare_means(method="t.test",comparisons=list(c("int","nat"))) +
  ylab("")
plot_grid(g1,g2,rel_widths=c(2.5,1))
ggsave("../result/sapro_symbio_comp.jpg")

g1 = ggplot(aes(x=fac_comp, y=TM_patho_sapro_symbio), data=pack_t_leaf) +
  geom_boxplot()
g2 = ggplot(aes(x=SMD_origin, y=TM_patho_sapro_symbio), data=pack_t_leaf) +
  geom_boxplot() +
  stat_compare_means(method="t.test",comparisons=list(c("int","nat"))) +
  ylab("")
plot_grid(g1,g2,rel_widths=c(2.5,1))
ggsave("../result/patho_sapro_symbio_comp.jpg")