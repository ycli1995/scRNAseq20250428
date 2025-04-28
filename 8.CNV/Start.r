options(future.globals.maxSize = 60 * 1024^3)
suppressWarnings(library(rlang, warn.conflicts = F))
suppressWarnings(library(yaml,  warn.conflicts = F))

# yml handler #
handlers = list('bool#no'  = function(x) if ( x %in% c('false', 'FALSE', 'no')  ) FALSE else x, 
                'bool#yes' = function(x) if ( x %in% c('true',  'TRUE',  'yes') ) TRUE  else x)

# word color #
suppressWarnings(library(crayon, warn.conflicts = F))
Pa    = crayon::cyan
Er    = crayon::red$bold
Sa    = crayon::blue
No    = crayon::magenta$bold
Wa    = crayon::yellow
timer = function() crayon::yellow(as.character(Sys.time()))

# make percentage #
makePct = function(x , digits = 3) {
  if(!anyNA(as.numeric(x))) x = paste0(round(as.numeric(x), digits = digits) *100, '%')
  x
}

# make minMax #
minMax = function(x) (x - min(x)) / (max(x) - min(x))

# cbinds #
cbinds = function(F1, F2, fill = 0) {
  if(any(dim(F1) == 0)) return(F2)
  if(any(dim(F2) == 0)) return(F1)
  rowall = union(rownames(F1), rownames(F2))
  if(length(setdiff(rowall, rownames(F1)))){
    SF1r           = matrix(fill, nrow = length(setdiff(rowall, rownames(F1))), ncol = ncol(F1))
    rownames(SF1r) = setdiff(rowall, rownames(F1))
    colnames(SF1r) = colnames(F1)
    F1             = rbind(F1, SF1r)
    rm(SF1r)
  }
  if(length(setdiff(rowall, rownames(F2)))){
    SF2r           = matrix(fill, nrow = length(setdiff(rowall, rownames(F2))), ncol = ncol(F2))
    rownames(SF2r) = setdiff(rowall, rownames(F2))
    colnames(SF2r) = colnames(F2)
    F2             = rbind(F2, SF2r)
    rm(SF2r)
  }
  F2 = F2[rownames(F1), , drop = F]
  cbind(F1, F2)
}

# change name #
change = function(i) gsub(' ', '-', gsub('\\.', '-', i))

## Pipeline latest version ##
linkdataLatestVersion  = 1.1
trustVDJLatestVersion  = 1
cellstateLatestVersion = 1.2
cnvLatestVersion       = 1.3
fastqIOLatestVersion   = 1.2

## Color ##
col20 = c('#00468B', '#5377A7', '#6C6DA4', '#925E9F', '#759EDD', 
          '#0099B4', '#42C1BB', '#76D1B1', '#0A7C2E', '#B8D24D', 
          '#EDE447', '#FAB158', '#FDAF91', '#FF7777', '#FD0000', 
          '#AD002A', '#AE8691', '#DEB8A1', '#4C4E4E', '#5B4232')

