
library(dplyr)
library(viridis)
library(patchwork)
library(networkD3)
source("/public2/Bio/pipeline/SingleCell_Collections/SCellWare/v3.1/R/Seurat_lib.R", chdir = TRUE)
tf = readLines("tf.list")
target = readLines("target.list")
edge = read.table("tf_target.egde.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

edge = edge[edge$tf %in% tf, ]
edge = edge[edge$target %in% target, ]

edge$target = paste0(edge$target, " ")
target = paste0(target, " ")

df = data.frame(
  TF = factor(edge$tf, tf),
  target = factor(edge$target, target),
  value = 1,
  stringsAsFactors = FALSE
)

color_genes = SetColor(factor(c(tf, target), levels = c(tf, target)), "onlinereport")

library(ggalluvial)
df_lodes <- to_lodes_form(df, key = "x", value = "stratum", id = "alluvium", axes = 1:2)

library(ggrepel)
p = ggplot(df_lodes, aes(x = x, stratum = stratum, alluvium = alluvium, fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_flow(width = 0.2, knot.pos = 0.1) +
  geom_stratum(alpha = .9, color = NA, width = 1/7) +
  geom_text(stat = "stratum", size = 3, color = "black") +
  scale_fill_manual(values = color_genes) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank())+
  guides(fill = FALSE)
ggsave(p, filename = "Sankey.pdf", width = 6, height = 10)

