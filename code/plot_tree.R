#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Usage plot_tree.R input output.n", call.=FALSE)
} 

library(ggtree)
library(ape)
library(tidyverse)
# source https://rdrr.io/bioc/treeio/src/R/nhx.R
get_nhx_feature <- function(nhx_features) {
  nameSET <- strsplit(nhx_features, split=":") %>% unlist %>%
    gsub("=.*", "", .) %>% unique
  lapply(nhx_features, get_nhx_feature_internal, nameSET=nameSET) %>%
    do.call(rbind, .) %>% as.data.frame(., stringsAsFactors = FALSE)
}

get_nhx_feature_internal <- function(feature, nameSET) {
  x <- strsplit(feature, ":") %>% unlist
  name <- gsub("=.*", "", x)
  val <- gsub(".*=", "", x)
  
  names(val) <- name
  y <- character(length(nameSET))
  for (i in seq_along(nameSET)) {
    if (nameSET[i] %in% name) {
      y[i] <- val[nameSET[i]]
    } else {
      y[i] <- NA
    }
  }
  names(y) <- nameSET
  return(y)
}

read.nhx <- function(file) {
treetext <- readLines(file, warn=FALSE)
treetext <- treetext[treetext != ""]
treetext <- treetext[treetext != " "]

if (length(treetext) > 1) {
  treetext <- paste0(treetext, collapse = '')
}
treetext <- gsub(" ", "", treetext)

phylo <- read.tree(text=treetext)
nnode <- Nnode(phylo, internal.only=FALSE)
nlab <- paste("X", 1:nnode, sep="")
tree2 <- treetext

pattern <- "(\\w+)?(:?\\d*\\.?\\d*[Ee]?[\\+\\-]?\\d*)?\\[&&NHX.*?\\]"
for (i in 1:nnode) {
  tree2 <- sub(pattern, paste0("\\", nlab[i], "\\2"), tree2)
}

phylo2 <- read.tree(text = tree2)
node <- match(nlab, sub(".+(X\\d+)$","\\1",
                        c(phylo2$tip.label, phylo2$node.label)))
node <- node[!is.na(node)] 
nhx.matches <- gregexpr(pattern, treetext)

matches <- nhx.matches[[1]]
match.pos <- as.numeric(matches)
if (length(match.pos) == 1 && (match.pos == -1)) {
  nhx_tags <- data.frame(node = 1:nnode)
} else {
  match.len <- attr(matches, 'match.length')
  
  nhx_str <- substring(treetext, match.pos, match.pos+match.len-1)
  
  nhx_features <- gsub("^[^\\[]*", "", nhx_str) %>%
    gsub("\\[&&NHX:", "", .) %>%
    gsub("\\]", "", .)
  
  nhx_tags <- get_nhx_feature(nhx_features)
  fields <- names(nhx_tags)
  for (i in ncol(nhx_tags)) {
    if(any(grepl("\\D+", nhx_tags[,i])) == FALSE) {
      ## should be numerical varialbe
      nhx_tags[,i] <- as.character(nhx_tags[,i])
    }
  }
  nhx_tags$node <- as.integer(node)
}

# Order rows by row number to facilitate downstream manipulations
nhx_tags <- nhx_tags[order(nhx_tags$node),]

new("treedata",
    #file = filename(file),
    phylo = phylo,
    data = as_tibble(nhx_tags)
)
}


lines = readLines(args[1])
treetext = lines[1]
tree <- read.nhx(textConnection(treetext))

P <- ggtree(tree) + geom_tiplab(geom = "label", size = 6, hjust = -0.2, fill = 'pink') + 
  geom_label(aes(x=branch, label=S), hjust = 0.55, vjust = 1, fill='lightgreen', size = 5) + 
  geom_label(aes(label=N), fill='lightblue', size = 6) 
outfile = args[2]
outfile = paste(outfile, "pdf", sep = ".")
ggsave(outfile, width = 49, height = 28, dpi = 320)