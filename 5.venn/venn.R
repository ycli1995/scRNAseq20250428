#!/Bio/bin/Rscript
# =======================
# linpeng
# plin@genedenovo.com
# Wed Feb 19 15:44:38 CST 2020
# =======================

args <- commandArgs()
bin <- dirname(sub('--file=', '',  args[grep('--file=', args)]))

# import packages
library(VennDiagram)
library(UpSetR)
library(reshape2)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Input a table [default %default]"),
  make_option(c("-c", "--compare"), type = "character", default = NULL,
              help = "Compare str (e.g. 'A&B&C') [default %default]"),
  make_option(c("-l", "--list"), type = "character", default = NULL,
              help = "List file (e.g. a.glist,b.glist), [default %default]"),
  make_option(c("-o", "--outdir"), type = "character", default = 'Venn',
              help = "Output dir [default %default]"),
  make_option(c("-d", "--outdir2"), type = "logical", default = T,
              help = "Output dir [default %default]"),
  make_option(c("-m", "--map"),  type = "character", default = NULL,
              help = "Group file [default %default]"),
  make_option(c("-a", "--annot"),  type = "character", default = NULL,
              help = "Annot file [default %default]"),
  make_option(c("-p", "--plot"),  type = "character", default = NULL,
              help = "Script to plot Venn [default %default]")
)
opts <- parse_args(OptionParser(usage = "%prog [options]", option_list = option_list), positional_arguments = F)

# int the fill colors 
my_colors = c("#FF8080","#807FFF","#BF4080","purple","yellow")

# read the group shared info
data <- data.frame()
pairwises <- c()
if (!is.null(opts$input) && !is.null(opts$compare)) {
  # read data table
  data <- read.table(opts$input,header=T,check=F,comment="",sep="\t",quote="")
  
  if (!is.null(opts$map)) {
    group <- read.table(opts$map,header=T,check=F,row=1,comment="",sep="\t",quote="")
    data.df <- melt(data)
    colnames(data.df)[1] <- "id"
    data.df$variable <- group[data.df$variable,1]
    data <- dcast(data.df, id~variable, sum )
  }
  row.names(data) <- data[, 1]
  data <- data[,-1]
  # read the compare
  pairwises <- unlist(strsplit(opts$compare, ","))

} else if(!is.null(opts$list)) {
  opts$outdir2 <- F
  files <- unlist(strsplit(opts$list, ","))
  arr <- sub(".g?list", "", basename(files))
  pairwises <- paste(arr, collapse = "&")

  id <- c()
  alist <- list()
  for (i in 1:length(arr)) {
    glist <- read.table(files[i], header = F, check = F, comment = "", sep = "\t", quote = "", stringsAsFactors = F)
    glist <- glist[, 1]
    alist[[i]] <- glist
    id <- c(id, glist)
  }
  id <- unique(id)
  data <- matrix(ncol = length(arr), nrow = length(id))
  for (i in 1:length(arr)) {
    data[, i] <- as.numeric(id %in% alist[[i]])
  }
  row.names(data) <- id
  colnames(data) <- arr

} else {
  usage = "
Version: 2.0
Descs:   An R script to batch draw venn diagram
Usage:   Rscript venn.R -i <table> -c <compare str> -m <group_file> -o <Venn>
         Rscript venn.R -l list,list,list 
Note:    compare str was defined the in the meta16 conf file, (two_group_diff or multi_group_diff)\n"
  cat(usage)
  quit("no")
}

print(str(pairwises))
print(str(data))

if (!is.null(opts$annot)) {
  opts$annot <- read.table(opts$annot,header=T,check=F,comment="",sep="\t",quote="")
  rownames(opts$annot) <- opts$annot[, 1]
  if (colnames(opts$annot)[ncol(opts$annot)] == 'Taxonomy') {
    opts$annot <- opts$annot[, c(1, ncol(opts$annot)), drop = F]
  }
}

# all otus and  all groups 
otus   = rownames(data)
groups = colnames(data)

# a function to draw venn diagram, create png (600 dpi) and pdf two format file
venn_plot <- function(x,outname,cex=2) {
  fill_colors = my_colors[1:length(x)]
  
  # png format 
  filename = paste0(outname,"/venn.png")
  #venn.diagram(x,filename=filename,col="white",fill=fill_colors,lwd=.5,cex=cex / 5,cat.cex=cex / 5, ext.text = TRUE, margin = 0,
  #                         alpha=.5,resolution=600, width=1200,height=1200,imagetype="png")
  # pdf format
  venn = venn.diagram(x,filename=NULL,col="white",fill=fill_colors,margin = 0, ext.text = T,
                      cex=cex,alpha=.5, cat.cex=cex)
  pdf_name = paste0(outname,"/venn.pdf")
  pdf(file=pdf_name,height=7,width=14)
  grid.newpage()
  pushViewport(viewport(width=unit(0.5, "npc")))
  grid.draw(venn)
  dev.off()
  system(paste("convert -density 300", pdf_name, filename))
}

# a function to draw venn diagram with UpSetR and output the elements of each sets
upset_plot <- function (x, outdir) {
  x.all = sort(unique(unlist(x)))
  x.mat = matrix(0, nrow = length(x.all), ncol = length(x), dimnames = list(genes=x.all, groups=names(x)))

  for(i in 1:length(x)) x.mat[x[[i]], i] <- 1
  
  combination <- nrow(unique(data.frame(x.mat)))
  width <- 8
  height <- 5
  if (combination > 5) width <- width + (combination - 5) * .2
  if (length(x) > 4) height <- height + (length(x) - 4) * .5
  if (length(x) > 5) width <- width + length(x) - 5
  png_name <- paste0(outdir, "/upset.png")
  pdf_name <- paste0(outdir, "/upset.pdf")

  text.scale = c(1.75, 1.5, 1.75, 1.5, 1.75, 2)
  if (ncol(x.mat) <= 8) {
    pdf(pdf_name, width = width, height = height, onefile=TRUE)
    print(upset(as.data.frame(x.mat), nsets = length(x), nintersects = 40, set_size.show = TRUE, set_size.numbers_size = 8, set_size.scale_max = max(lengths(x)) * 1.75, text.scale = text.scale, point.size = 5, line.size = 2), newpage = FALSE)
    dev.off()
    png(png_name, width = width, height = height, units = 'in', res = 250)
    print(upset(as.data.frame(x.mat), nsets = length(x), nintersects = 40, set_size.show = TRUE, set_size.numbers_size = 8, set_size.scale_max = max(lengths(x)) * 1.75, text.scale = text.scale, point.size = 5, line.size = 2))
    dev.off()
    #cmd <- paste('/usr/bin/convert -density 300', pdf_name, png_name)
    #system(cmd)
  } else {
    print('[WARN] compare group > 8, ignore upset')
  }

  # output the elements
  x.venn <- split(rownames(x.mat), apply(x.mat, 1, paste, collapse=""))
  names_lst <- 
    lapply(names(x.venn), function(.ele) 
    colnames(x.mat)[as.logical(as.numeric(strsplit(.ele, "")[[1]]))])
  
  names    = sapply(names_lst, function(x) paste(x, collapse='_'))
  numbers  = rep(0,length(names))
  elements = rep("",length(names))

  print(str(x.venn))
  for (i in 1:length(names)) {
    num <- length(names_lst[[i]])
    if (ncol(x.mat) >= 5 & num > 1 & num < ncol(x.mat)) next
    f   <- paste0(outdir, '/', names[i], '.annot.xls')
    lst <- x.venn[[i]]
    arr <- names_lst[[i]]
    if (!is.null(opts$annot)) {
      if (!is.null(opts$list)) {
        dat <- opts$annot[lst, , drop=F]
      } else {
        dat <- data.frame(opts$annot[lst,  1, drop=F], data[lst, arr, drop=F], 
                          opts$annot[lst, -1, drop=F], check.names = F)
      }
      write.table(dat, file = f, row = F, col = T, quote = F, sep = "\t")
    } else if (!is.null(opts$input) & is.null(opts$plot)) {
      dat <- data.frame(id = lst, data[lst, arr, drop=F], check.names = F)
      write.table(dat, file = f, row = F, col = T, quote = F, sep = "\t")
    }
  }
}

outdir <- opts$outdir
if(!dir.exists(outdir)) dir.create(outdir, recursive = T)
outdir <- normalizePath(outdir)

for (i in pairwises) {
  compares <- unlist(strsplit(i,"&"))
  x = list()
  for (gname in compares){
    x[[gname]] <- otus[as.numeric(data[,gname]) > 0]
  }
  if (opts$outdir2) {
    name <- gsub(pattern="&",replacement="_vs_",i)
    name <- paste(outdir,name,sep="/")
    if(!dir.exists(name)) dir.create(name)
  } else {
    name <- outdir
  }
  
  # prepare data for venn plot
  dat <- data.frame(Groups = compares)
  dat$Tags <- apply(data[, compares], 2, function(x) paste(otus[x!=0], collapse = "\t") )
  fil <- paste0(name, "/input.list")
  write.table(dat, file = fil, sep = "\t", quote = F, col = F, row = F)

  if (!is.null(opts$plot) && opts$plot != 'sbv') {
    outprefix <- paste0(name, '/venn')
    cmd <- paste('perl', opts$plot, 'venn', '-file', fil, '-outprefix', outprefix, '-width 900 -height 900 -margin 50')
    print(cmd)
    system(cmd)
    unlink(fil)
  } else if (!is.null(opts$plot) && opts$plot == 'sbv') {
    svg <- paste0(name, "/venn.svg")
    png <- paste0(name, "/venn.png")
    sbv <- paste0(bin, "/../SBV/bin/sbv.pl")
    cof <- paste0(bin, "/venn.conf")
    cmd <- paste('perl', sbv, 'venn', '-conf', cof, '-out', svg, fil)
    cat('[CMD]', cmd, '\n')
    system(cmd)
    system(paste('convert', svg, png))
    #warning("VennDiagram can't draw >5 sets venn diagram, skip ...")
  } else {
    venn_plot(x,name)
  }
  upset_plot(x, name)
  unlink(paste0(name, "/*.log"))
  unlink("VennDiagram*.log")
}

