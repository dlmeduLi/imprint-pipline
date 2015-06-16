library("data.table")
library("reshape2")
library("ggplot2")
library("BSgenome.Mmusculus.UCSC.mm9")

PlotTagMethRange <- function(genome, tagmeth.data, chr.name, range.start, range.end) {
  
  # get cg positions
  
  refseq <- genome[[chr.name]][range.start:range.end]
  cg.pos <- range.start + start(matchPattern("CG", refseq)) - 2
  
  # get tags in the range
  
  tagmeth.range <- tagmeth.data[tagmeth.data$chr == chr.name & 
                                  (tagmeth.data$pos + tagmeth.data$len / 2) > range.start &
                                  (tagmeth.data$pos + tagmeth.data$len / 2) < range.end]
  
  # convert tag cg positions
  
  tagmeth <- tagmeth.range
  tagmeth$methylated <- lapply(strsplit(tagmeth.range$methylated, ","), as.integer)
  tagmeth$unmethylated <- lapply(strsplit(tagmeth.range$unmethylated, ","), as.integer) 
  tagmeth$undetermined <- lapply(strsplit(tagmeth.range$undetermined, ","), as.integer)
  
  tagmeth$meth.len <- unlist(lapply(tagmeth$methylated, length))
  tagmeth$unmeth.len <- unlist(lapply(tagmeth$unmethylated, length))
  tagmeth$undt.len <- unlist(lapply(tagmeth$undetermined, length))
  
  tagmeth$meth.count <- tagmeth$meth.len
  tagmeth$meth.count[is.na(tagmeth$methylated)] <- 0
  tagmeth$unmeth.count <- tagmeth$unmeth.len
  tagmeth$unmeth.count[is.na(tagmeth$unmethylated)] <- 0
  tagmeth$undt.count <- tagmeth$undt.len
  tagmeth$undt.count[is.na(tagmeth$undetermined)] <- 0
  
  tagmeth$meth.score <- tagmeth$meth.count / (tagmeth$meth.count + tagmeth$unmeth.count)
  tagmeth <- tagmeth[order(-tagmeth$meth.score), ]
  tagmeth$tagname <- NULL
  tagmeth$tagname <- c(1:nrow(tagmeth))
  
  # prepare plot data
  
  tagmeth.plot.data <- tagmeth[, list(index = match(c(unlist(methylated), unlist(unmethylated), unlist(undetermined)) + pos, cg.pos), 
                                      type = c(rep("M", meth.len), rep("U", unmeth.len), rep("X", undt.len))),
                               by = tagname]
  tagmeth.plot.data <- tagmeth.plot.data[!is.na(tagmeth.plot.data$index)]
  
  # plot
  
  ggplot(tagmeth.plot.data, aes(index, tagname, shape = type)) + geom_point() + scale_shape_manual(values=c(19, 1, 0))
  
}

tagmeth.filename <- "~/work/imprinting/41-WT_CTTGTA_L006_001.csv"
genome <- Mmusculus

# Read csv file

tagmeth.data <- fread(tagmeth.filename, sep="\t",  
                         colClasses=c("factor","character","character", "character",
                                      "integer","integer", "character", "character", "character"))
chr.name <- "chr19"
range.start <- 5406369
range.end <- 5407700

PlotTagMethRange(genome, tagmeth.data, chr.name, range.start, range.end)


