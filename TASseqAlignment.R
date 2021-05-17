# Date: 13 March 2019
# Thomas van Emden
# Sequence alignments of pNSU21 and pNSU70
#------------------------------------------------------------------------
# clear workspace
rm(list = ls())

# Install needed libraries
#install.packages("devtools")
#library(devtools)
#install_github("evolvedmicrobe/dotplot", build_vignettes = FALSE)
#install.packages("ggplot2")
#install.packages("seqRFLP")

# Load libraries
library("ggplot2")
library("dotplot")
#library("Biostrings")
library("seqRFLP")

# Define date for naming of plots
#datum <-format(Sys.Date(),format="%y%m%d")

# Define function from https://github.com/evolvedmicrobe/dotplot
# This gives control over the ggplot settings (probably also possible through package but I don't know how)
dotPlotg <- function (seq1, seq2, wsize = 10, wstep = 1, nmatch = -1)
{
  if (length(seq1[1]) > 1)
    stop("seq1 should be provided as a single string")
  if (length(seq2[1]) > 1)
    stop("seq2 should be provided as a single string")
  if (wsize < 1)
    stop("non allowed value for wsize")
  if (wstep < 1)
    stop("non allowed value for wstep")
  if (nmatch < 1)
    nmatch = wsize
  if (nmatch > wsize)
    stop("nmatch > wsize is not allowed")
  xy <- mkDotPlotDataFrame(seq1, seq2, wsize, wstep, nmatch)
  ggplot2::ggplot(xy, ggplot2::aes(x=x, y=y)) + 
    ggplot2::geom_point(shape=15, size = 0.01) + #+ # size=... is what I like to change
    ggplot2::theme(text=element_text(size=12,  family="sans")) +
    ggplot2::xlim(550, 700) +
    ggplot2::ylim(550, 700)
}

#------------------------------------------------------------------------
# Load in sequences of the different pNSUs
pnsu21 <- readChar("pNSU21.txt", file.info("pNSU21.txt")$size)
pnsu70 <- readChar('pNSU70.txt', file.info('pNSU70.txt')$size)

# Make sequence a bit shorter for quicker script
#pnsu21 <- substr(pnsu21, 1,200)
#pnsu70 <- substr(pnsu70, 2851, 2957)

# Set arguments for the plots
wsize = 87
nmatch = 58
dotsize = 0.01

# Make the plots
pdf("21x21.pdf")
dotPlotg(pnsu21, pnsu21, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="pNSU21", y = "pNSU21", title="pNSU21 x pNSU21", size= dotsize)
dev.off()

pdf("70x21.pdf")
dotPlotg(pnsu70, pnsu21, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="pNSU70", y = "pNSU21", title="pNSU70 x pNSU21", size= dotsize)
dev.off()

pdf("70x70.pdf", width = 8.6/2.54, height = 8.6/2.54)
dotPlotg(pnsu70, pnsu70, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="pNSU70", y = "pNSU70", title="pNSU70 x pNSU70", size= dotsize)
dev.off()

#print(theme_bw())
#text = c("This is text", "Text with\nmultiple lines", "Some more text")

#------------------------------------------------------------------------
# Redefine function for this part, probably not a good idea, but quick fix and just make sure you run whole script
dotPlotg <- function (seq1, seq2, wsize = 10, wstep = 1, nmatch = -1)
{
  if (length(seq1[1]) > 1)
    stop("seq1 should be provided as a single string")
  if (length(seq2[1]) > 1)
    stop("seq2 should be provided as a single string")
  if (wsize < 1)
    stop("non allowed value for wsize")
  if (wstep < 1)
    stop("non allowed value for wstep")
  if (nmatch < 1)
    nmatch = wsize
  if (nmatch > wsize)
    stop("nmatch > wsize is not allowed")
  xy <- mkDotPlotDataFrame(seq1, seq2, wsize, wstep, nmatch)
  ggplot2::ggplot(xy, ggplot2::aes(x=x, y=y)) + 
    ggplot2::geom_point(shape=15, size = 0.01) + #+ # size=... is what I like to change
    ggplot2::theme(text=element_text(size=12,  family="sans")) +
    ggplot2::xlim(0, 50000)+
    ggplot2::ylim(0, 50000)
}

chrI <- gsub("[\n]", "", readChar("chromosome1_ASM294v2.txt", file.info("chromosome1_ASM294v2.txt")$size))
chrII <- gsub("[\n]", "", readChar("chromosome2_ASM294v2.txt", file.info("chromosome2_ASM294v2.txt")$size))

Tel1L <- paste0(paste(replicate(9000, "B"), collapse = ""), substr(chrI,1,41000))
Tel1Lb <- gsub("B", "D", Tel1L)
Tel1R <- paste0(paste(replicate(17000, "D"), collapse = ""), revComp(substr(chrI,5546134,5579133)))
Tel1Rb <- gsub("D", "E", Tel1R)
Tel2L <- paste0(paste(replicate(28571, "E"), collapse = ""), substr(chrII, 1, 50000-28571))
Tel2Lb <- gsub("E", "B", Tel2L)
Tel2R <- revComp(substr(chrII,4489805,4539804))

write((paste0(">Tel1L\n", substr(chrI,1,41000))),"Tel1L.txt")
write((paste0(">Tel1R\n", revComp(substr(chrI,5546134,5579133)))),"Tel1R.txt")
write((paste0(">Tel2L\n", substr(chrII, 1, 50000-28571))),"Tel2L.txt")
write((paste0(">Tel2R\n", revComp(substr(chrII,4489805,4539804)))),"Tel2R.txt")

write(paste0(">Tel1L\n", 
             substr(chrI,1,41000),
             "\n\n",
             ">Tel1R\n",
             revComp(substr(chrI,5546134,5579133)),
             "\n\n",
             ">Tel2L\n",
             substr(chrII, 1, 50000-28571),
             "\n\n",
             ">Tel2R\n",
             revComp(substr(chrII,4489805,4539804))),
             "All_Tels.txt"
             )

print(nchar(Tel2R))
#revComp(substr(chrI,5529134,5579133))) original values I'd like to keep is things go south
# Simple format
#pdf("1Rx2R.pdf", width = 4/2.54, height = 4/2.54)
#dotPlotg(Tel1R, Tel2R, wsize=wsize, nmatch = nmatch) + 
#  theme_bw() + 
#  labs(x="Tel1R", y = "Tel2R", title="Tel1R x Tel2R", size= dotsize)
#dev.off()

# 1L column
pdf("1Lx1R.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel1L, Tel1R, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  #labs(x="Tel1L", y = "Tel1R", title="Tel1L x Tel1R", size= dotsize)+
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
        )
dev.off()

pdf("1Lx2L.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel1L, Tel2L, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

pdf("1Lx2R.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel1L, Tel2R, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

# 1R column
pdf("1Rx1L.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel1R, Tel1L, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

pdf("1Rx2L.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel1R, Tel2L, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

pdf("1Rx2R.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel1R, Tel2R, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

# 2L column
pdf("2Lx1R.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel2L, Tel1R, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  #labs(x="Tel2L", y = "Tel1R", title="Tel2L x Tel1R", size= dotsize)+
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

pdf("2Lx1L.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel2L, Tel1L, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

pdf("2Lx2R.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel2L, Tel2R, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

# 2R column
pdf("2Rx1L.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel2R, Tel1L, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

pdf("2Rx2L.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel2R, Tel2L, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

pdf("2Rx1R.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel2R, Tel1R, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

pdf("2Rx2R.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel2R, Tel2R, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

pdf("2Lx2L.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel2L, Tel2Lb, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

pdf("1Rx1R.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel1R, Tel1Rb, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

pdf("1Lx1L.pdf", width = 17.8/4/2.54, height = 17.8/4/2.54)
dotPlotg(Tel1L, Tel1Lb, wsize=wsize, nmatch = nmatch) + 
  theme_bw() + 
  labs(x="", y = "", title="", size= dotsize)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        title = element_blank()
  )
dev.off()

