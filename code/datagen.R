# This R script transform raw dataset (stored in .csv files) into easy-to-load from (.RData)
# In the downstream analysis, we won't handle rawdata, but work with .RData generated here

library(tidyverse)
library(adespatial)
library(SoDA)
library(vegan)

create.dbmem = function (coord = NULL, D.mat = NULL, nsites, randtest=T) {
  if (is.null(coord) & is.null(D.mat)) 
    stop("Geographic information must be provided in 'coord' or in 'D.mat'")
  if (is.null(D.mat)) 
    D.mat <- dist(coord)
  D.mat <- as.matrix(D.mat)
  if (!is.null(coord)) {
    n <- nrow(coord)
  }
  else {
    n <- nrow(D.mat)
  }
  if (sum(nsites) != n) 
    stop("Vector nsites does not sum to nrow(coord) or nrow(D.mat)")
  if (min(nsites) == 1) 
    stop("At least one group contains a single site")
  out <- matrix(0, n, n)
  end <- 0
  end.mem <- 0
  for (k in 1:length(nsites)) {
    start <- end + 1
    end <- end + nsites[k]
    tmp <- as.dist(D.mat[start:end, start:end])
    res <- dbmem(tmp, MEM.autocor = "positive")
    if (randtest) {
      test.mem = moran.randtest(res, attributes(res)$listw)
      tb_res = as_tibble(res)
      sel_res = select(tb_res, test.mem$names[test.mem$pvalue<0.05])
    }
    dbMEM <- as.matrix(sel_res)
    n.mem <- ncol(dbMEM)
    out[start:end, (end.mem + 1):(end.mem + n.mem)] <- dbMEM
    end.mem <- end.mem + n.mem
  }
  out <- out[, 1:end.mem]
  if (is.null(rownames(coord))) 
    rownames(out) <- rownames(out, do.NULL = FALSE, prefix = "Site.")
  else rownames(out) <- rownames(coord)
  colnames(out) <- colnames(out, do.NULL = FALSE, prefix = "dbMEM.")
  out
}

alpha_div = function(x, prefix) {
  otu.table = x %>% select(starts_with("OTU"))
  ret = tibble(
    richness = specnumber(otu.table),
    shannon = diversity(otu.table),
    evenness = shannon / log(richness),
    depth = rowSums(otu.table),
    richness.rarefied = rarefy(otu.table, min(depth)))
  colnames(ret) = paste0("DIV_", prefix, ".",colnames(ret))
  ret
}

# load in sample metadata
smd = read_csv(
  file = "../data/sample_metadata.csv", 
  col_types = cols(
    .default = col_double(),
    SMD_pid = col_character(),
    SMD_origin = col_factor(),
    SMD_site = col_character(),
    SMD_population = col_factor()), 
  progress=F)

# load in OTU tables
leaf = read_csv(
  file = "../data/otu_table_leaf.csv", 
  col_types=cols(
    .default = col_double(),
    SMD_sid = col_character(),
    SMD_pid = col_character()),
  progress=F)

root = read_csv(
  file = "../data/otu_table_root.csv", 
  col_types=cols(
    .default = col_double(),
    SMD_sid = col_character(),
    SMD_pid = col_character()),
  progress=F)

soil = read_csv(
  file = "../data/otu_table_soil.csv", 
  col_types=cols(
    .default = col_double(),
    SMD_sid = col_character(),
    SMD_pid = col_character()),
  progress=F)

# calculate dbMEMs and bind them into smd
cartesian_coord = geoXY(smd$SMD_lat, smd$SMD_lon, unit=1000)
sep = T
if (sep) {
  dbMEMs = create.dbmem(coord=cartesian_coord, nsites=c(420, 420), randtest=T)
  nmem.nat = sum(dbMEMs[1,]!=0)
  nmem.int = sum(dbMEMs[421,]!=0)
  dbMEMs.nat = dbMEMs[,1:nmem.nat]
  dbMEMs.int = dbMEMs[,(nmem.nat+1):ncol(dbMEMs)]
  colnames(dbMEMs.nat) = paste0("dbMEMs.nat.",1:nmem.nat)
  colnames(dbMEMs.int) = paste0("dbMEMs.int.",1:nmem.int)
  smd = cbind(smd, dbMEMs.nat, dbMEMs.int)
} else {
  dbMEMs = dbmem(cartesian_coord, thresh=270, MEM.autocor = "positive")
  test.mem = moran.randtest(dbMEMs, attributes(dbMEMs)$listw)
  selected_dbMEMs = select(dbMEMs, test.mem$names[test.mem$pvalue<0.05])
  smd = cbind(smd, selected_dbMEMs)
}

clim_predictor = c("V1", "V4", "V5", "V6", "V12", "V14", "V15", "V18")
cpc = smd %>%
  select(all_of(clim_predictor)) %>%
  scale() %>%
  prcomp()
# climate PC1 and PC2
smd$cpc1 = cpc$x[,1]
smd$cpc2 = cpc$x[,2]

# calculate alpha diversity for each of three community
leaf = cbind(leaf, alpha_div(leaf, "leaf"))
root = cbind(root, alpha_div(root, "root"))
soil = cbind(soil, alpha_div(soil, "soil"))

smd = smd %>% 
  left_join(select(leaf, SMD_pid, starts_with("DIV")), by="SMD_pid") %>%
  left_join(select(root, SMD_pid, starts_with("DIV")), by="SMD_pid") %>%
  left_join(select(soil, SMD_pid, starts_with("DIV")), by="SMD_pid")

leaf = select(leaf, SMD_pid, starts_with("OTU"))
root = select(root, SMD_pid, starts_with("OTU"))
soil = select(soil, SMD_pid, starts_with("OTU"))

# drop S5 S6 S7
to_drop = (smd %>% filter(SMD_site %in% c("S5", "S6", "S7")))$SMD_pid
smd = smd %>% filter(!SMD_pid %in% to_drop)
leaf = leaf %>% filter(!SMD_pid %in% to_drop)
root = root %>% filter(!SMD_pid %in% to_drop)
soil = soil %>% filter(!SMD_pid %in% to_drop)

# drop empty OTU
leaf = leaf[,colSums(leaf==0)!=nrow(leaf)]
root = root[,colSums(root==0)!=nrow(root)]
soil = soil[,colSums(soil==0)!=nrow(soil)]

# assign guild annotations
guild = read_tsv("../data/taxonomy.guilds.tsv")
guild = select(
  guild,
  otu_name = "OTU ID",
  trophic_mode="Trophic Mode", 
  guild="Guild", 
  growth_morph="Growth Morphology"
)

otu_leaf = leaf %>%
  select(starts_with("OTU")) %>%
  t() %>%
  as.data.frame()

# drop low abundance OTUs
otu_leaf = otu_leaf[rowSums(otu_leaf) > 5,]
colnames(otu_leaf) = leaf$SMD_pid
otu_leaf = as.data.frame(otu_leaf > 0)
otu_leaf$otu_name = rownames(otu_leaf)

# merge OTUs
trophic_mode_leaf = otu_leaf %>%
  left_join(guild, by="otu_name") %>%
  group_by(trophic_mode) %>%
  select(-otu_name, -guild, -growth_morph) %>%
  summarise_all(sum)

# add prefix
trophic_mode = paste0("TM_", c("unassigned", "patho", "patho_sapro", "patho_sapro_symbio", "patho_symbio", "sapro", "sapro_symbio", "symbio"))

# transpose
trophic_mode_leaf = trophic_mode_leaf %>% 
  select(-trophic_mode) %>%
  t() %>%
  as.data.frame()
colnames(trophic_mode_leaf) = trophic_mode

trophic_mode_leaf$SMD_pid = rownames(trophic_mode_leaf)
t_leaf = trophic_mode_leaf
##

pack_leaf = left_join(leaf, smd, by="SMD_pid")
otu_leaf = select(pack_leaf, starts_with("OTU")) 
otu_leaf = otu_leaf[rowSums(otu_leaf) > 5,] > 0
ca_leaf = cca(otu_leaf)
ss = ca_leaf$CA$u
hc = cbind(ss[,1], ss[,2]) %>%
  dist() %>%
  hclust(method="ward.D")
cl = cutree(hc, 7)
pack_leaf$compartment = cl
smd = smd %>%
  left_join(select(pack_leaf, SMD_pid, compartment), by="SMD_pid")

save(smd, leaf, t_leaf, root, soil, file="../data/full_data_set.RData")