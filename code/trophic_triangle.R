library(tidyverse)

load("../data/full_data_set.RData")

calc_centroids_trophic_triangle = function(df_trophic_mode) {
  
  # according to biplot(pct)
  
  all_patho = df_trophic_mode$TM_patho + 
    df_trophic_mode$TM_patho_sapro / 2 +
    df_trophic_mode$TM_patho_symbio / 2 +
    df_trophic_mode$TM_patho_sapro_symbio / 3
  
  all_symbio = df_trophic_mode$TM_symbio + 
    df_trophic_mode$TM_sapro_symbio / 2 +
    df_trophic_mode$TM_patho_symbio / 2 +
    df_trophic_mode$TM_patho_sapro_symbio / 3

  all_sapro = df_trophic_mode$TM_sapro + 
    df_trophic_mode$TM_patho_sapro / 2 +
    df_trophic_mode$TM_sapro_symbio / 2 +
    df_trophic_mode$TM_patho_sapro_symbio / 3

  pid = df_trophic_mode$SMD_pid
  
  x = data.frame(
    patho = all_patho,
    symbio = all_symbio,
    sapro = all_sapro
  ) %>% scale()
  
  center = sweep(x, 2, apply(x, 2, min),'-')
  R = apply(x, 2, max) - apply(x, 2, min)
  x_star = sweep(center, 2, R, "/") %>% as.data.frame()
  
  theta = (c(1,2,3)-1) * 2*pi/3
  xm = sweep(x_star, 2, cos(theta), "*")
  ym = sweep(x_star, 2, sin(theta), "*")
  
  centroids = t(apply(cbind(xm, ym), 1, function(x) pracma::poly_center(x[1:3], x[4:6])))
  area = apply(cbind(xm, ym), 1, function(x) pracma::polyarea(x[1:3], x[4:6]))
  
  df = tibble(
    area = area,
    cent.x = centroids[,1],
    cent.y = centroids[,2],
    SMD_pid = pid
  )
  x_star$SMD_pid = pid
  ret = list(df_cent=df, xs = x_star)
}

f = function(l) {
  geom_polygon(aes(
    x = c(l,l,l)*cos((c(1,2,3)-1) * 2*pi/3),
    y = c(l,l,l)*sin((c(1,2,3)-1) * 2*pi/3),
    group=c(1,1,1)
  ), color="gray", fill="white",alpha=0)
}

f2 = function(df.cl, term) {
  df = select(df.cl, all_of(term), compartment)
  df.test = pairwise.wilcox.test(df[[term]], df$compartment)
  df.test = as.data.frame(df.test$p.value)
  df.test$row = as.factor(2:7)
  df.test = pivot_longer(df.test, cols=1:6, names_to="col") %>% na.omit()
  
  mat = matrix(NA, 6, 6)
  for (c in 1:6) {
    for (r in (c+1):7) {
      col = filter(df, compartment==c)
      row = filter(df, compartment==r)
      delta = mean(col[[term]]) - mean(row[[term]])
      mat[r-1,c] = delta
    }
  }
  mat = as.data.frame(mat)
  colnames(mat) = 1:6
  mat$row = 2:7
  mat = pivot_longer(mat, cols=1:6, names_to="col") %>% na.omit()
  mat$row = as.factor(mat$row)
  mat$col = as.factor(mat$col)
  
  df = left_join(mat, df.test, by=c("row","col")) %>%
    filter(row %in% c(5,6,7) & col %in% c(1,2,3,4))
  df$row = ordered(df$row, c(2,3,1,4,6,7,5))
  df$col = ordered(df$col, c(2,3,1,4,6,7,5))
  df$sig = cut(df$value.y, c(1, 0.05, 0), labels=c("*",""), right=F, include.lowest=T)
  ggplot(df, aes(col, row)) +
    geom_tile(aes(fill=value.x),color="black") + 
    geom_text(aes(label=sig),vjust=0.8,size=12) +
    ggtitle(term) +
    scale_fill_gsea()+
    theme_bw()
  ggsave(str_glue("../result/",term,"_test.jpg"))
}

ret = calc_centroids_trophic_triangle(t_leaf)

xs = as_tibble(ret$xs)
cent = ret$df_cent
df = left_join(xs, cent, by="SMD_pid")

df.cl = left_join(df, select(smd, SMD_pid, compartment), by="SMD_pid")
df.cl.sum = df.cl %>% select(-SMD_pid) %>% group_by(compartment) %>% summarise_all(mean)

for (i in c(1,2,3,4,5,6,7)) {
  theta = (c(1,2,3)-1) * 2*pi/3
  ggplot() + 
    f(1) + f(0.75) + f(0.5) + f(0.25) +
    geom_segment(aes(x=0,y=0,xend=1*cos(theta[1]),yend=1*sin(theta[1]))) +
    geom_segment(aes(x=0,y=0,xend=1*cos(theta[2]),yend=1*sin(theta[2]))) +
    geom_segment(aes(x=0,y=0,xend=1*cos(theta[3]),yend=1*sin(theta[3]))) +
    geom_text(aes(x=1*cos(theta[1])-0.05,y=1*sin(theta[1])-0.1,label="Patho"))+
    geom_text(aes(x=1*cos(theta[2])+0.15,y=1*sin(theta[2]),label="Symbio"))+
    geom_text(aes(x=1*cos(theta[3])+0.15,y=1*sin(theta[3]),label="Sapro"))+
    geom_polygon(aes(
      x = as.double(df.cl.sum[i,2:4])*cos(theta),
      y = as.double(df.cl.sum[i,2:4])*sin(theta),
      group=c(1,1,1)
    ), fill=rainbow(7)[i], alpha=0.5) +
    geom_point(aes(x=as.double(df.cl.sum[i,6]),y=as.double(df.cl.sum[i,7]))) +
    theme_bw() +
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()
    ) +
    ggtitle(str_glue("Compartment ",i))
  ggsave(str_glue("../result/tro_tri_",i,".jpg"))
}

for (term in c("area", "cent.x", "cent.y", "patho", "sapro", "symbio")) {
  f2(df.cl, term)
}

library(lmerTest)

host_predictor = c(
  "TRA_leafherb", "TRA_leafbio", 
  "TRA_SSL", "TRA_SLA", "TRA_BI",
  "TRA_trichome", "TRA_saponins", "TRA_lignin")#, "TRA_CN")
clim_predictor = c("cpc1", "cpc2")
soil_predictor = c(
  "SOL_water", "SOL_PH", 
  "SOL_SOM", "SOL_SCN", 
  "SOL_NH", "SOL_NO", "SOL_AP")
all_predictor = c(clim_predictor, host_predictor, soil_predictor)

pack_c_leaf = df.cl %>% 
  left_join(select(smd, all_of(all_predictor), starts_with("SMD")), by="SMD_pid") %>% 
  na.omit()

sformula = str_glue(
  "area ~ ",
  paste0(all_predictor, collapse=" + "),
  " + ",
  paste0(all_predictor, ":SMD_origin", collapse = " + "),
  "+(1|SMD_site)")
fit=  lmer(as.formula(sformula), data=pack_c_leaf)

sformula = str_glue(
  "cbind(cent.x, cent.y) ~ ",
  paste0(all_predictor, collapse=" + "),
  " + ",
  paste0(all_predictor, ":SMD_origin", collapse = " + "))



lm(area ~ SMD_origin + as.factor(compartment) + SMD_site, data=pack_c_leaf) %>% anova()
lm(cent.x ~ SMD_origin + as.factor(compartment) + SMD_site, data=pack_c_leaf) %>% anova()
lm(cent.y ~ SMD_origin + as.factor(compartment) + SMD_site, data=pack_c_leaf) %>% anova()