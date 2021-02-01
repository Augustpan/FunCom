library(tidyverse)

# summarize vp_tro_c
summ = matrix(0, 4, 7)
for (comp in c(1,2,3,4,5,6,7)) {
  df = read.csv(str_glue("../result/vp_tro_c", comp, ".csv"))
  df$Adj.R.squared[df$Adj.R.squared < 0] = 0
  summ[,comp] = df$Adj.R.squared
}
summ = as.data.frame(summ)
colnames(summ) = paste0("comp", 1:7)
rownames(summ) = df$X
write.csv(summ, "../result/vp_tro_c_summ.csv")

# summarize vp_otu_c
summ = matrix(0, 4, 7)
for (comp in c(1,2,3,4,5,6,7)) {
  df = read.csv(str_glue("../result/vp_otu_c", comp, ".csv"))
  df$Adj.R.squared[df$Adj.R.squared < 0] = 0
  summ[,comp] = df$Adj.R.squared
}
summ = as.data.frame(summ)
colnames(summ) = paste0("comp", 1:7)
rownames(summ) = df$X
write.csv(summ, "../result/vp_otu_c_summ.csv")

# summarize rda_otu_c
summ = matrix(0, 19, 14)
for (comp in c(1,2,3,4,5,6,7)) {
  df = read.csv(str_glue("../result/rda_otu_c", comp, ".csv"))
  summ[,comp*2-1] = df$Variance
  summ[,comp*2] = df$Pr..F.
}
summ = as.data.frame(summ)
colnames(summ) = paste0("comp", 1:7)
rownames(summ) = df$X
write.csv(summ, "../result/rda_otu_c_summ.csv")

# summarize rda_tro_c
summ = matrix(0, 19, 14)
for (comp in c(1,2,3,4,5,6,7)) {
  df = read.csv(str_glue("../result/rda_tro_c", comp, ".csv"))
  summ[,comp*2-1] = df$Variance
  summ[,comp*2] = df$Pr..F.
}
summ = as.data.frame(summ)
colnames(summ) = paste0("comp", 1:7)
rownames(summ) = df$X
write.csv(summ, "../result/rda_tro_c_summ.csv")

df1 = read.csv("../result/ss1_otu.csv")
df2 = read.csv("../result/ss2_otu.csv")
summ = data.frame(term=df1$X, estimate_1=df1$Estimate, p_value_1=df1$Pr...t..,
                  estimate_2=df2$Estimate, p_value_2=df2$Pr...t..)
write_csv(summ, "../result/ss_otu.csv")

df1 = read.csv("../result/ss1_tro.csv")
df2 = read.csv("../result/ss2_tro.csv")
summ = data.frame(term=df1$X, estimate_1=df1$Estimate, p_value_1=df1$Pr...t..,
                  estimate_2=df2$Estimate, p_value_2=df2$Pr...t..)
write_csv(summ, "../result/ss_tro.csv")