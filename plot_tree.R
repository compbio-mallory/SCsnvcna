library(ggtree)
library(ape)
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

file <- textConnection(treetext)
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

treetext = "((|1|3|:2.2[&&NHX:S=:L=:M=|+20|+21|+22|+47|+48|+54|+55|+66|+70|+71|+86|:N=CO8PA12],
(|10|11|12|28|:2.6[&&NHX:S=+4|:L=:M=|+6|+53|+59|+60|+61|+62|+63|+64|+65|+80|+81|-83|+85|:N=CO8PA24],
(|2|4|7|:5.2[&&NHX:S=:L=:M=|-0|-1|-3|-4|-5|-14|-15|-16|-26|-27|-28|-29|-30|-39|-40|-41|-42|-43|-44|-45|-46|-68|-69|-78|+89|+91|:N=CO8PA13],
(|:1.2[&&NHX:S=:L=:M=|+7|+20|-52|+74|+75|+79|:N=CO8PA19],
(|:1.8[&&NHX:S=:L=:M=|-3|-4|-5|+72|+80|-87|-88|-90|+91|:N=CO8PA2],
(|:0.6[&&NHX:S=:L=:M=|+79|+80|+81|:N=CO8PA10],
((|0|17|21|22|:0.6[&&NHX:S=:L=-2|:M=|-6|-53|-73|:N=CO8PA23],
((|:0.2[&&NHX:S=:L=:M=|+89|:N=CO8PA6],
(|:0[&&NHX:S=:L=:M=|:N=CO8PA3],
|6|8|20|:3.2[&&NHX:S=:L=-3|:M=|-5|+6|-21|-22|-23|-24|-25|-53|+68|+72|+73|-78|-79|-80|-81|+89|:N=CO8PA16]):
0.8[&&NHX:S=:L=:M=|-1|-6|-73|+85|:N=25]):
0.8[&&NHX:S=:L=:M=|+54|+55|+78|+79|:N=26],
(|:0.6[&&NHX:S=:L=:M=|+7|+68|-73|:N=CO8PA8],
(|13|:1.4[&&NHX:S=:L=:M=|+54|+55|-79|-87|-88|+89|-90|:N=CO8PA17],
(|:1.0[&&NHX:S=:L=:M=|+5|-6|+75|-76|-77|:N=CO8PA20],
|5|9|14|24|25|26|27|29|30|31|32|33|:1.2[&&NHX:S=:L=:M=|+45|-50|-51|-52|+85|+91|:N=CO8PA14]):
0.8[&&NHX:S=:L=:M=|-73|+74|+82|+83|:N=22]):
0.8[&&NHX:S=:L=:M=|-1|-45|-46|-53|:N=23]):
0.4[&&NHX:S=+0|:L=:M=|-80|-81|:N=24]):
0.4[&&NHX:S=:L=:M=|-5|+79|:N=27],
(|:0.6[&&NHX:S=:L=:M=|-5|+72|+74|:N=CO8PA7],
|:0.6[&&NHX:S=:L=:M=|-87|-88|-90|:N=CO8PA11]):
0.4[&&NHX:S=:L=:M=|-79|+85|:N=21]):
0.4[&&NHX:S=:L=:M=|-82|-83|:N=28],
(|:0.8[&&NHX:S=:L=:M=|-1|+47|+48|+49|:N=CO8PA5],
(|:2.6[&&NHX:S=:L=:M=|-6|+7|+16|-53|+68|-73|+75|-78|-79|-80|-81|+89|+91|:N=CO8PA18],
|15|16|18|19|23|:1.2[&&NHX:S=:L=:M=|-17|-18|-19|+43|-46|+74|:N=CO8PA1]):
1.0[&&NHX:S=:L=:M=|-31|-32|-33|+74|+85|:N=19]):
0.2[&&NHX:S=:L=:M=|+78|:N=20]):
0.6[&&NHX:S=:L=:M=|+65|-68|-89|:N=29]):
0.8[&&NHX:S=:L=:M=|+6|+8|+79|-81|:N=30]):
2.8[&&NHX:S=+3|:L=:M=|+21|+22|+53|+59|+60|+61|+62|+63|+64|+78|+80|+84|+86|+89|:N=31]):
4.6[&&NHX:S=+1|:L=:M=|-8|+9|+10|+11|+12|+13|+20|+23|+24|+25|+31|+32|+33|+34|+36|+37|+38|-72|+73|+80|+81|+82|+83|:N=32]):
2.6[&&NHX:S=+2|7|:L=:M=|-2|-7|+35|+36|+37|+38|+52|-56|-57|-58|-67|+76|+77|:N=33]):
4.6[&&NHX:S=:L=:M=|+3|+4|+5|+14|+15|+16|+26|+27|+28|+29|+30|+52|+53|+72|+73|+74|+75|+78|+80|+81|+87|+88|+90|:N=34]):
17.4[&&NHX:S=+8|9|14|6|5|10|13|11|12|15|:L=:M=|+0|+1|+2|+3|+4|+5|+6|+7|+8|+9|+10|+11|+12|+13|+14|+15|+16|+17|+18|+19|+20|+21|+22|+23|+24|+25|+26|+27|+28|+29|+30|+31|+32|+33|+34|+35|+36|+37|+38|+39|+40|+41|+42|+43|+44|+45|+46|+47|+48|+49|+50|+51|+52|+53|+54|+55|+56|+57|+58|+59|+60|+61|+62|+63|+64|+65|+66|+67|+68|+69|+70|+71|+72|+73|+74|+75|+76|+77|+78|+79|+80|+81|+82|+83|+84|+85|+86|:N=35])[&&NHX:N=diploid];"
tree <- read.nhx(textConnection(treetext))
ggtree(tree) + geom_tiplab(geom = "label", fill = 'pink') + 
  geom_label(aes(x=branch, label=S), hjust = 0, vjust = 1.7, fill='lightgreen') + 
  geom_label(aes(x=branch, label=L), hjust = 0, vjust = -0.7, fill='gray') +
  #geom_label(aes(x=branch, label=M), hjust= 0, vjust = 0) +
  geom_label(aes(label=N), fill='steelblue') 

