#!/usr/bin/env Rscript
# version 0.0.1
#########################################################################################
# intersection.R                                                                        #
#########################################################################################
#                                                                                       #
# This program is free software: you can redistribute it and/or modify it under the     #
# terms of the GNU General Public License as published by the Free Software Foundation, #
# either version 3 of the License, or (at your option) any later version.               #
#                                                                                       #
#########################################################################################
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, ##
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A       ##
## PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT  ##
## HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION   ##
## OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE      ##
## SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                              ##
#########################################################################################

# Gobals informations
author = "Corentin Hochart"
email = "corentin.hochart.pro@gmail.com"
version = '0.0.1'

# Initialisation
if (!require("optparse")) install.packages("optparse")
args <- commandArgs(trailingOnly = F)
script.path <- dirname(sub("--file=","",args[grep("--file",args)]))[1]

outputDir<-paste(Sys.Date(),'intersection',sep = '_')

# Get options
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Input file matrix ", metavar="FILE"),
  make_option(c("-p", "--pairwise"), type="character", default = "n",
              help="Input file is a pairwise matrix [y/N]", metavar="[y/n]"),
  make_option(c("-c", "--correlation_coeficient"), type="character", default="kendal",
              help="correlation coefficient [kendal]", metavar="METHOD"),
  make_option(c("-a","--agglomeration"), type="character", default="complete",
              help="agglomeration for hclust [complete]", metavar="METHOD"),
  make_option(c("-d","--distance"), type="character", default="euclidian",
              help="distance matrix computation [euclidian]", metavar="METHOD"),
  make_option(c("--cluster_number"), type="numeric", default = 4,
              help="Cluster number define to color the dendogram [4]", metavar="INTEGER"),
  make_option(c("-o","--output"), type="character", default="pairwise_heatmap.png",
              help="Ouptut file in png format [pairwise_heatmap.png]", metavar="FILE")
);

opt_parser = OptionParser(usage = "Usage: %prog -f [FILE]",option_list=option_list,
                          description= "
Description: Produce a pairwise intersection heatmap with dendrogram from features abundances by sample data.
                          
Written by Corentin Hochart (corentin.hochart@gmail.com).
Released under the terms of the GNU General Public License v3. 
intersection.R version 0.0.1.")

opt = parse_args(opt_parser)

# Control arguments
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("No input file.")
}

# Librairies
if (!require("gplots")) install.packages("gplots")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("dendextend")) install.packages("dendextend")

# Main code

## Create pairwise matrix 
if (opt$pairwise == "n" || opt$pairwise == "N"){
  matrix<-read.csv(file=opt$file,sep = ",",row.names=1,check.names=FALSE)
  
  matrix[matrix > 0] <- 1
  write.table(matrix,
              file = paste(outputDir,'PresenceAbscenceTable.csv',sep=''),
              sep=",")
  
  pairwise_matrix<-matrix(0,nrow = ncol(matrix), ncol = ncol(matrix))
  colnames(pairwise_matrix)<-colnames(matrix)
  rownames(pairwise_matrix)<-colnames(matrix)
  for (i in 1:nrow(matrix)){
    cluster<-rownames(matrix[i,])
    for (j in 1:ncol(matrix)){
      S1<-colnames(matrix[i,])[j]
      for (k in 1:ncol(matrix)){
        S2<-colnames(matrix[i,])[k]
        if (matrix[i,j] == 1 && matrix[i,k] == 1){
          pairwise_matrix[j,k]<-pairwise_matrix[j,k] + 1
        }
      }
    }
  }
  
  write.table(pairwise_matrix,
              file = paste(outputDir,'pairwiseTable.csv',sep=''),
              sep=",")
} else {
  pairwise_matrix<-read.csv(file=opt$file,sep = ",",row.names=1,check.names=FALSE)
}

## Compute correlation 
if (opt$correlation_coeficient != 'None'){
  pairwise_matrix<-cor(pairwise_matrix,method = opt$correlation_coeficient)
} else {
  pairwise_matrix<-as.matrix(pairwise_matrix)  
}

## Compute distance matrix and cluster 
distance<-dist(pairwise_matrix, method = opt$distance)    
hcluster<-hclust(distance, method = opt$agglomeration)

## Create dendogram object 
dend1<-as.dendrogram(hcluster)
colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))
dend1 <- color_branches(dend1, k = opt$cluster_number, col = colfunc(opt$cluster_number))
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

## Plot heatmap and dendogram 
printf <- function(...) print(sprintf(...))
title <- printf("Pairwise analysis\n(%s/%s/%s)",
                opt$correlation_coeficient,opt$distance, opt$agglomeration)

png(opt$output,width = 8, height = 8,  units = 'in', res = 300)
heatmap.2(pairwise_matrix,
          # dendrogram control
          Rowv = dend1,
          Colv = dend1,
          breaks = seq(-1, 1,0.015),
          # title
          main = title,
          
          # image plot
          col = colorRampPalette(brewer.pal(9,"RdYlBu"))(133),
          
          # block sepration
          colsep=1:ncol(pairwise_matrix),
          rowsep=1:nrow(pairwise_matrix),
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          
          # Row/Column Labeling
          RowSideColors = col_labels,  
          ColSideColors = col_labels,
          colRow = col_labels,
          colCol = col_labels,
          
          # level trace
          trace = "none" ,
          
          # color key + density info
          key = TRUE,
          keysize = 1.2,
          key.title = "Color Key",
          key.par = list(cex=0.5)
          # density.info="density"
)
dev.off()
