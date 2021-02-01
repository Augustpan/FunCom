library(tidyverse)
library(vegan)
library(lmerTest)

load("../data/full_data_set.RData")

host_predictor = c(
  "TRA_leafherb", "TRA_leafbio", 
  "TRA_SSL", "TRA_SLA", "TRA_BI",
  "TRA_trichome", "TRA_saponins", "TRA_lignin", "TRA_CN")
clim_predictor = c("cpc1", "cpc2")
soil_predictor = c(
  "SOL_water", "SOL_PH", 
  "SOL_SOM", "SOL_SCN", 
  "SOL_NH", "SOL_NO", "SOL_AP")
all_predictor = c(clim_predictor, host_predictor, soil_predictor)

pack_leaf = leaf %>% 
  left_join(smd, by="SMD_pid") %>% 
  na.omit()
otu_leaf = select(pack_leaf, starts_with("OTU"))
otu_leaf = otu_leaf[rowSums(otu_leaf) > 5,] > 0

# CA
ca_leaf = cca(otu_leaf)
fit = envfit(ca_leaf, select(pack_leaf, all_of(all_predictor)), perm = 999)

ss = ca_leaf$CA$u
hc = cbind(ss[,1], ss[,2]) %>%
  dist() %>%
  hclust(method="ward.D")
cl = cutree(hc, 7)
ggplot() +
  geom_point(aes(x=ss[,1],y=ss[,2],color=as.factor(cl))) +
  theme_bw()
ggsave("../result/ca_otu.jpg")

jpeg("../result/ca_otu_env.jpg", width=720, height=720)
plot(ca_leaf)
plot(fit)
dev.off()

sformula = str_glue(
  "ss[,1] ~ SMD_origin + ",
  paste0(all_predictor, collapse=" + "),
  " + ",
  paste0(all_predictor, ":SMD_origin", collapse = " + "),
  " + (1|SMD_site)")
lmer.fit = lmer(as.formula(sformula),data=pack_leaf) %>% summary()
write.csv(lmer.fit$coefficients, "../result/ss1_otu.csv")

sformula = str_glue(
  "ss[,2] ~ SMD_origin + ",
  paste0(all_predictor, collapse=" + "),
  " + ",
  paste0(all_predictor, ":SMD_origin", collapse = " + "),
  " + (1|SMD_site)")
lmer.fit = lmer(as.formula(sformula),data=pack_leaf) %>% summary()
write.csv(lmer.fit$coefficients, "../result/ss2_otu.csv")

# CCA
sformula = str_glue(
  "otu_leaf ~ ",
  paste0(all_predictor, collapse="+")
)
cca_leaf = cca(as.formula(sformula), data=pack_leaf)
aov.cca = anova(cca_leaf,by="term", parallel=8)
site.cca = summary(cca_leaf)$sites
var.cca = summary(cca_leaf)$biplot * 2.1

ggplot() +
  geom_point(aes(x=site.cca[,1],y=site.cca[,2],color=as.factor(cl))) + 
  geom_segment(
    aes(x=0,y=0,xend=var.cca[,1],yend=var.cca[,2]), 
    arrow=arrow(angle=22.5, length=unit(0.2,"cm"), type="closed")) +
  geom_text(aes(x=var.cca[,1],y=var.cca[,2], label=rownames(var.cca)), hjust=0.2, vjust=1.5) +
  theme_bw()
ggsave("../result/cca_otu.jpg")
write.csv(aov.cca, "../result/cca_otu.csv")

# NMDS
nmds = metaMDS(otu_leaf)
ggplot() +
  geom_point(aes(x=nmds$points[,1],y=nmds$points[,2],color=as.factor(cl))) +
  theme_bw()
ggsave("../result/nmds_otu.jpg")

# RDA
otu_leaf = decostand(otu_leaf, method="hellinger")
sformula = str_glue(
  "otu_leaf ~ SMD_origin + ",
  paste0(all_predictor, collapse=" + "),
  " + ",
  paste0(all_predictor, ":SMD_origin", collapse = " + "))
rda_leaf = rda(as.formula(sformula), data=pack_leaf)
aov.rda = anova(rda_leaf, by="term", parallel=8)
write.csv(aov.rda, "../result/rda_otu.csv")

vp = varpart(
  otu_leaf,
  select(pack_leaf, all_of(clim_predictor),all_of(soil_predictor)),
  select(pack_leaf, all_of(host_predictor))
)
write.csv(vp$part$indfract, "../result/vp_otu.csv")

for (comp in c(1,2,3,4,5,6,7)) {
  pack_leaf_f = filter(pack_leaf, compartment==comp)
  otu_leaf_f = pack_leaf_f %>% 
    select(starts_with("OTU")) %>%
    decostand(method="hellinger")
  sformula = str_glue(
    "otu_leaf_f ~ ",
    paste0(all_predictor, collapse=" + "))
  rda_leaf = rda(as.formula(sformula), data=pack_leaf_f)
  aov.rda = anova(rda_leaf, by="term", parallel=8)
  write.csv(aov.rda, str_glue("../result/rda_otu_c",comp,".csv"))
  
  vp = varpart(
    otu_leaf_f,
    select(pack_leaf_f, all_of(clim_predictor),all_of(soil_predictor)),
    select(pack_leaf_f, all_of(host_predictor))
  )
  write.csv(vp$part$indfract, str_glue("../result/vp_otu_c",comp,".csv"))
}

