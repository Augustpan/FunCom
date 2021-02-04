library(tidyverse)

load("../data/full_data_set.RData")

ncol(leaf %>% select(starts_with("OTU")))

nrow(leaf)

leaf %>% select(starts_with("OTU")) %>% sum()

guild = read_tsv("../data/taxonomy.guilds.tsv")
guild = select(
  guild,
  otu_name = "OTU ID",
  trophic_mode="Trophic Mode", 
  guild="Guild", 
  growth_morph="Growth Morphology"
)

cs = select(leaf, starts_with("OTU")) %>% colSums()
otu = names(cs[cs>5])
f = filter(guild, otu_name %in% otu)
filter(, trophic_mode=="-") %>% nrow()


s_leaf = left_join(leaf, smd, by="SMD_pid")
nat = filter(s_leaf, SMD_origin=="nat") %>% select(starts_with("OTU"))
nat= nat[,colSums(nat) > 0] %>% colnames()

int = filter(s_leaf, SMD_origin=="int") %>% select(starts_with("OTU"))
int = int[,colSums(int) > 0] %>% colnames()

intersect(nat, int) %>% length()
setdiff(nat, int) %>% length()
setdiff(int, nat) %>% length()

