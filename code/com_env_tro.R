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

pack_t_leaf = t_leaf %>% 
  left_join(select(smd, all_of(all_predictor), starts_with("SMD"), starts_with("dbMEM"), compartment), by="SMD_pid") %>% 
  na.omit()
tro_leaf = select(pack_t_leaf, starts_with("TM"))

# CA
ca_t_leaf = cca(tro_leaf)
fit = envfit(ca_t_leaf, select(pack_t_leaf, all_of(all_predictor)), perm = 999)

ss = ca_t_leaf$CA$u
cl = as.factor(pack_t_leaf$compartment)
ggplot() +
  geom_point(aes(x=ss[,1],y=ss[,2],color=cl)) +
  theme_bw()
ggsave("../result/ca_tro.jpg", width=7.43, height=6.44)

jpeg("../result/ca_tro_env.jpg", width=720, height=720)
plot(ca_t_leaf)
plot(fit)
dev.off()

sformula = str_glue(
  "ss[,1] ~ SMD_origin + ",
  paste0(all_predictor, collapse=" + "),
  " + ",
  paste0(all_predictor, ":SMD_origin", collapse = " + "),
  " + (1|SMD_site)")
lmer.fit = lmer(as.formula(sformula),data=pack_t_leaf) %>% summary()
write.csv(lmer.fit$coefficients, "../result/ss1_tro.csv")

sformula = str_glue(
  "ss[,2] ~ SMD_origin + ",
  paste0(all_predictor, collapse=" + "),
  " + ",
  paste0(all_predictor, ":SMD_origin", collapse = " + "),
  " + (1|SMD_site)")
lmer.fit = lmer(as.formula(sformula),data=pack_t_leaf) %>% summary()
write.csv(lmer.fit$coefficients, "../result/ss2_tro.csv")

# CCA
sformula = str_glue(
  "tro_leaf ~ ",
  paste0(all_predictor, collapse="+")
)
cca_t_leaf = cca(as.formula(sformula), data=pack_t_leaf)
aov.cca = anova(cca_t_leaf,by="term", parallel=8)
site.cca = summary(cca_t_leaf)$sites
var.cca = summary(cca_t_leaf)$biplot * 2.1

ggplot() +
  geom_point(aes(x=site.cca[,1],y=site.cca[,2],color=as.factor(cl))) + 
  geom_segment(
    aes(x=0,y=0,xend=var.cca[,1],yend=var.cca[,2]), 
    arrow=arrow(angle=22.5, length=unit(0.2,"cm"), type="closed")) +
  geom_text(aes(x=var.cca[,1],y=var.cca[,2], label=rownames(var.cca)), hjust=0.2, vjust=1.5) +
  theme_bw()
ggsave("../result/cca_tro.jpg")
write.csv(aov.cca, "../result/cca_tro.csv")

# NMDS
nmds = metaMDS(tro_leaf)
ggplot() +
  geom_point(aes(x=nmds$points[,1],y=nmds$points[,2],color=as.factor(cl))) +
  theme_bw()
ggsave("../result/nmds_tro.jpg")

## RDA
#tro_leaf = decostand(tro_leaf, method="hellinger")
sformula = str_glue(
  "tro_leaf ~ SMD_origin + ",
  paste0(all_predictor, collapse=" + "),
  " + ",
  paste0(all_predictor, ":SMD_origin", collapse = " + "))
rda_t_leaf = rda(as.formula(sformula), data=pack_t_leaf)
aov.rda = anova(rda_t_leaf, by="term", parallel=8)
write.csv(aov.rda, "../result/rda_tro.csv")

vp = varpart(
  tro_leaf,
  select(pack_t_leaf, all_of(clim_predictor),all_of(soil_predictor)),
  select(pack_t_leaf, all_of(host_predictor))
)
write.csv(vp$part$indfract, "../result/vp_tro.csv")

for (comp in c(1,2,3,4,5,6,7)) {
  pack_t_leaf_f = filter(pack_t_leaf, compartment==comp)
  tro_leaf_f = pack_t_leaf_f %>% 
    select(starts_with("TM")) #%>%
    #decostand(method="hellinger")
  sformula = str_glue(
    "tro_leaf_f ~ ",
    paste0(all_predictor, collapse=" + "))
  rda_leaf = rda(as.formula(sformula), data=pack_t_leaf_f)
  aov.rda = anova(rda_leaf, by="term", parallel=8)
  write.csv(aov.rda, str_glue("../result/rda_tro_c",comp,".csv"))
  
  vp = varpart(
    tro_leaf_f,
    select(pack_t_leaf_f, all_of(clim_predictor),all_of(soil_predictor)),
    select(pack_t_leaf_f, all_of(host_predictor))
  )
  write.csv(vp$part$indfract, str_glue("../result/vp_tro_c",comp,".csv"))
}
