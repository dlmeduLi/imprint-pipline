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

##sliding windows
swsCalc<-function(x, win=list(L=75, R=75)){
  library(zoo)
  xLen <- length(x)
  winLen<-win$L + win$R + 1
  sws <- rollsum(x,winLen)
  sws_head<-tail(cumsum(x[1:(winLen-1)]),win$L)
  sws_tail<-rev(tail(cumsum(rev(x[(xLen - winLen + 2):xLen])),win$R))
  sws <- append(sws_head, sws)
  sws <- append(sws, sws_tail)
  return(sws)
}

tagmeth.filename <- "~/work/imprinting/41-WT_CTTGTA_L006_001.csv"
genome <- Mmusculus

# Read csv file

tagmeth.data <- fread(tagmeth.filename, sep="\t",  
                         colClasses=c("factor","character","character", "character",
                                      "integer","integer", "character", "character", "character"))
chr.name <- "chr1"
range.start <- 95721123
range.end <- 95723462

# PlotTagMethRange(genome, tagmeth.data, chr.name, range.start, range.end)

# calculate the midpoint of tags

tagmeth.data$mid <- round(tagmeth.data$pos + tagmeth.data$len / 2)

# calculate methy score

tagmeth <- tagmeth.data
tagmeth$methylated <- lapply(strsplit(tagmeth.data$methylated, ","), as.integer)
tagmeth$unmethylated <- lapply(strsplit(tagmeth.data$unmethylated, ","), as.integer) 
tagmeth$undetermined <- lapply(strsplit(tagmeth.data$undetermined, ","), as.integer)

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

# build tagmeth vector

vector.chr <- c(1:length(genome[[chr.name]]))
tagmeth.counter <- data.frame(pos = vector.chr)

low.score <- 0.2
high.score <- 0.8

tagmeth.low <- tagmeth[tagmeth$chr == chr.name & tagmeth$meth.score < low.score]
tagmeth.high <- tagmeth[tagmeth$chr == chr.name & tagmeth$meth.score > high.score]

tagmeth.low.stat <- as.data.frame(table(tagmeth.low$mid))
tagmeth.low.stat[,1] <- as.integer(as.character(tagmeth.low.stat[,1]))
tagmeth.high.stat <- as.data.frame(table(tagmeth.high$mid))
tagmeth.high.stat[,1] <- as.integer(as.character(tagmeth.high.stat[,1]))

#tagmeth.low$count <- rep(1, nrow(tagmeth.low))
#tagmeth.low.stat <- dcast(tagmeth.low, mid ~ count, sum, value.var = "")
#tagmeth.high$count <- rep(1, nrow(tagmeth.high))
#tagmeth.high.stat <- dcast(tagmeth.high, mid ~ count, sum)

tagmeth.counter$low.count <- rep(0, nrow(tagmeth.counter))
tagmeth.counter$low.count[tagmeth.low.stat[,1]] <- tagmeth.low.stat[,2]
tagmeth.counter$high.count <- rep(0, nrow(tagmeth.counter))
tagmeth.counter$high.count[tagmeth.high.stat[,1]] <- tagmeth.high.stat[,2]

low.sum <- swsCalc(tagmeth.counter$low.count, win=list(L=25, R=25))
high.sum <- swsCalc(tagmeth.counter$high.count, win=list(L=25, R=25))

# write wig file

wig.header <- paste("fixedStep chrom=", chr.name, " start=1 step=1 span=1", sep = "")
low.file = file("~/work//imprinting/low.wig", "w")
writeLines(c(wig.header, low.sum), low.file)
close(low.file)

high.file = file("~/work//imprinting/high.wig", "w")
writeLines(c(wig.header, high.sum), high.file)
close(high.file)