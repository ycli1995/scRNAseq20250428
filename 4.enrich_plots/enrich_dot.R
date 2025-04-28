#!/public2/Bio/pipeline/Toolkit_ycli/bin/Rscript4_ycli

args = commandArgs(TRUE)

infile = args[1]
outpfx = args[2]

library(ggplot2)
library(ggpubr)
library(RColorBrewer)

color_name = "Pvalue"

dt = data.table::fread(infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

#dt = dt[order(dt$ratio), ]
#dt$Descrption = factor(dt$Descrption, unique(dt$Descrption))

dt$class[dt$class %in% "Biological Process"] = "BP"
dt$class[dt$class %in% "Molecular Function"] = "MF"
dt$class[dt$class %in% "Cellular Component"] = "CC"

all_go_class = c("BP", "MF", "CC")
if (any(all_go_class %in% dt$class)) {
  dt$class = factor(dt$class, intersect(all_go_class, dt$class))  
}
dt$Count = dt$num
dt$color = dt[, color_name]

dt = dt[order(dt$ratio), ]
dt$Descrption = factor(as.character(dt$Descrption), unique(dt$Descrption))

p1 = ggplot(dt, aes(x = ratio, y = Descrption)) +
  geom_point(aes(size = Count, color = color)) +
  scale_colour_gradient(low = "red", high = "blue", name = color_name) +
  scale_x_continuous(limits = c(0, max(dt$ratio + 0.01))) +
  scale_size_continuous(range = c(1, 5), limits = c(0, max(dt$Count)))
if (is.factor(dt$class)) {
  p1 = p1 + facet_wrap(~class, ncol = 1, strip.position = "right", scales = "free_y")
}
p1 = p1 + ylab("") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    strip.text = element_text(size = 14),
    panel.grid = element_blank()
  ) + 
  Seurat::RotatedAxis()
h = max(5, 1.2 + length(levels(dt$Descrption)) * 0.12)
w = max(6, 6 + max(nchar(levels(dt$Descrption))) * 0.01)
ggsave(p1, filename = paste0(outpfx, ".DotPlot.pdf"), height = h, width = w)

mytheme <- theme(
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 12),
  plot.title = element_text(size = 12,
  hjust = 0.5,
  face = "bold"),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5,
    r = 10,
    l = 5.5,
    b = 5.5)
)

dt = dt[order(dt[, color_name]), ]
dt$Descrption = factor(as.character(dt$Descrption), unique(dt$Descrption))
max_y1 = ceiling(max(-log10(dt$color)))
scale_factor = ceiling(max(dt$Count) / max_y1)
p2 = ggplot() +
  geom_bar(data = dt, aes(x = -log10(color), y = Descrption), stat = "identity", width = 0.8, fill = 'pink') + 
  scale_x_continuous(
    limits = c(0, max_y1),
    sec.axis = sec_axis(~. * scale_factor + 50, name = "Count")
  ) +
  geom_line(data = dt, aes(x = Count / scale_factor, y = Descrption), color = "grey", group = 1) +
  geom_point(data = dt, aes(y = Descrption, x = Count / scale_factor)) + 
  labs(y = "", x = paste0("-log10(", color_name, ")")) +
  theme_classic() + mytheme
if (is.factor(dt$class)) {
  p2 = p2 + facet_wrap(~class, ncol = 1, strip.position = "right", scales = "free_y")
}
h = max(5, 1.2 + length(levels(dt$Descrption)) * 0.12)
w = max(5, 6 + max(nchar(levels(dt$Descrption))) * 0.01)
ggsave(p2, filename = paste0(outpfx, ".BarPlot.pdf"), height = h, width = w)

