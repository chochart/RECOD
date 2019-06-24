#!/usr/bin/env Rscript
# version 0.0.1
#########################################################################################
# CCA_ggplot.R                                                                          #
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
start_time <- Sys.time()
if (!require("optparse")){install.packages("optparse",dep = TRUE ,repos = "http://cran.us.r-project.org")}
args <- commandArgs(trailingOnly = F)
script.path <- dirname(sub("--file=","",args[grep("--file",args)]))[1]

# Get options
option_list = list(
  make_option(c("-a", "--abundTable"), type="character", default=NULL,
              help="Features abundance table ", metavar="FILE"),
  make_option(c("-e", "--envTable"), type="character", default=NULL,
              help="Environmental data table ", metavar="FILE"),
  make_option(c("-m", "--metaTable"), type="character", default=NULL,
              help="Meta data table ", metavar="FILE"),
  make_option(c("-p","--prefix"), type="character", default='CCA_plot',
              help="Output prefix for svg file [Default: 'CCA_plot' ]", metavar="FILE"),
  make_option(c("--features_color_by"), type="character", default = NULL,
              help="Parameter to color features"),
  make_option(c("--sites_shape_by"), type="character", default = NULL,
              help="Parameter to shape features.
                Must be in this format 'FieldNumber,FieldName'"),
  make_option(c("--site_seperator_field"), type="character", default = '_',
              help="Field seperator for sites [Default: '_' ]"),
  make_option(c("--significant_environmental_variable"), type="logical", action="store_true", default=FALSE,
              help="Use adonis to find significant environmental variables and only use these in cca. Use
            all environmmental variables by default"),
  make_option(c("--significant_co"), type="numeric", default=0.01,
              help="Cutoff for Pr(>F) significance probability value associated with the F Value 
              [Default: 0.01]"),
  make_option(c("--scaling"), type="integer", default = 0,
              help="Scaling for species and site scores. Either species (2) or site (1) scores are scaled 
            by eigenvalues, and the other set of scores is left unscaled, or with 3 both are scaled 
            symmetrically by square root of eigenvalues. This scaling is know as Hill scaling. The 
            type of scores are 'none', 'sites', 'species', or 'symmetric', which correspond to the 
            values 0, 1, 2, and 3 respectively. For more details help(scores) [Default: '0' ]"),
  make_option(c("--scaling_arrows"), type="logical", action="store_true", default=FALSE,
                help="Scaling arrows using vegan::ordiArrowMul() function. Arrows are scaled by their
              correlation (square root of the column r2) so in that way '1weak' predictors have shorter
              arrows than 'strong' predictors"),
  make_option(c("--arrowmul"), type="character", default=NULL,
              help="Expand arrows in the graph. Arrows will be scaled automatically to fit the graph 
            {NULL,min,max,optimal} [Default :NULL]"
  ),
  make_option(c("-s", "--separator"), type="character", default='comma',
              help="Input files fields separator {comma|tab} [comma]"),
  make_option(c("-r", "--reproductible"), type="logical", action="store_true", default=FALSE,
              help="Enable reproductibility of the analysis")
);

opt_parser = OptionParser(usage = "usage: %prog -a [ABUNDTABLE] -m [envTable] -a [metaTable]",option_list=option_list,
                          description= "
Description: Produce CCA plot from a features abundance, meta-data and features annotation files. 

Written by Corentin Hochart (corentin.hochart.pro@gmail.com).
Released under the terms of the GNU General Public License v3. 
CCA_plot.R version 0.0.1.
Inspired by R code for ecological data analysis by Umer Zeeshan Ijaz.
")

opt = parse_args(opt_parser)

# Controle positionnal argument
if (is.null(opt$abundTable)){
  print_help(opt_parser)
  stop("No abundance file.")
}
file <- opt$abundTable
if( file.access(file) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", file))
}

if (is.null(opt$envTable)){
  print_help(opt_parser)
  stop("No meta-data file.")
}
file <- opt$envTable
if( file.access(file) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", file))
}

if (is.null(opt$metaTable)){
  print_help(opt_parser)
  stop("No annotated features file.")
}
file <- opt$metaTable
if( file.access(file) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", file))
}

if(opt$separator == 'comma'){
  opt$separator = ','
} else if(opt$separator == 'tab'){
  opt$separator = '\t'
} else{
  stop(sprintf("Uncorrect specified separator ( %s ). Must be {comma|tab}", opt$separator))
}

# Libraries
if (!require("ggplot2")){install.packages("ggplot2",dep = TRUE ,repos = "http://cran.us.r-project.org")}
if (!require("gridExtra")){install.packages("gridExtra",dep = TRUE ,repos = "http://cran.us.r-project.org")}
if (!require("svglite")){install.packages("svglite",dep = TRUE ,repos = "http://cran.us.r-project.org")}
if (!require("vegan")){install.packages("vegan",dep = TRUE ,repos = "http://cran.us.r-project.org")}

# Functions
draw<- function(df1,df2,df3,fcb,ssh,abscisse,ordinate,arrowmul){
  
  p<-ggplot()
  ### Plot sites 
  if(is.null(ssh)){
    shape_by = c("Sites")
    shape_legend = "Sites"
  } else{
    shape<-strsplit(ssh,",")
    indice<-as.numeric(shape[[1]][1])+2
    shape_by<-factor(df1[,indice])
    shape_legend = shape[[1]][2]
  }
  
  p<-p+geom_point(data=df1,aes(x,y,shape=shape_by), color = 'darkgrey')
  p<-p+scale_shape_manual(values=1:nlevels(shape_by))

  ### Define Theme
  p<-p+theme(
    panel.border = element_rect(colour = "grey", fill=NA),
    panel.background = element_blank()
  )
  
  xlabel=paste('CCA Axis',abscisse,sep=' ')
  ylabel=paste('CCA Axis',ordinate,sep=' ')
  
  ### Draw features
  if(is.null(fcb)){
    color_by = c("Features")
    color_legend = "Features"
  } else{
    color_by<-df3[,fcb]
    color_legend = fcb
  }
  
  p1<-p+geom_point(data=df3,aes(x,y,color=color_by))
  p1<-p1+labs(x=xlabel, y=ylabel,colour=color_legend, shape=shape_legend)
  
  ### Draw biplots
  if(!is.null(arrowmul)){
    xmindf<-abs(min(df1$x,df3$x))
    xmaxdf<-max(df1$x,df3$x)
    xmindf2<-abs(min(df2$x))
    xmaxdf2<-max(df2$x)
    xminMulti<-xmindf/xmindf2
    xmaxMulti<-xmaxdf/xmaxdf2
    
    ymindf<-abs(min(df1$y,df3$y))
    ymaxdf<-max(df1$y,df3$y)
    ymindf2<-abs(min(df2$y))
    ymaxdf2<-max(df2$y)
    yminMulti<-ymindf/ymindf2
    ymaxMulti<-xmaxdf/ymaxdf2

    multi<-1
    if(arrowmul=='min'){
      multi<-min(c(xminMulti,xmaxMulti,yminMulti,ymaxMulti))
    } else if (arrowmul=='max'){
      multi<-max(c(xminMulti,xmaxMulti,yminMulti,ymaxMulti))
    } else if (arrowmul=='max') {
      multi<-min(c(max(xminMulti,xmaxMulti),max(yminMulti,ymaxMulti)))
    }
    df2<-df2*multi
  }
  
  p1<-p1+geom_segment(data=df2, aes(x = 0, y = 0, xend = x, yend = y),
                    arrow = arrow(length = unit(0.2, "cm")))
  
  p1<-p1+geom_text(data=as.data.frame(df2*1.1),aes(x, y, label = rownames(df2)))

  if(i == 1 || j == 2){
    # Create features text plot only for x = CCA1 and y = CCA2 
    
    ### Draw features
    p2<-p+geom_text(data=df3,aes(x,y,label=color_by,color=color_by),
                    hjust=0, vjust=0,size=3)
    p2<-p2+labs(x=xlabel, y=ylabel,colour=color_legend, shape=shape_legend)

    ### Draw biplots
    p2<-p2+geom_segment(data=df2, aes(x = 0, y = 0, xend = x, yend = y),
                      arrow = arrow(length = unit(0.2, "cm")))
    
    p2<-p2+geom_text(data=as.data.frame(df2*1.1),aes(x, y, label = rownames(df2)))
      
  }else{
    p2<-NULL
  }

  return(list(p1,p2))
  
}

logger<-function(message){
  today <- Sys.time()
  message(sprintf("[%s] INFO: %s",today,message))
}

# Main code 
user<-Sys.getenv("USERNAME")
logger(sprintf("Hi %s! Let's do some good jobs together.",user))

if (opt$reproductible){
  set.seed(2)
}

## Load data
logger("Load data")
abund_table<-read.csv(opt$abundTable,row.names=1,check.names=FALSE,sep = opt$separator)
# Filter out any features that have zero entries
abund_table<-abund_table[apply(abund_table[,-1], 1, function(x) !all(x==0)),]
abund_table[is.na(abund_table)]<-0
abund_table<-t(abund_table)
# Filter out any samples that have zero entries
abund_table<-abund_table[apply(abund_table[,-1], 1, function(x) !all(x==0)),]

env_table<-read.csv(opt$envTable,row.names=1,check.names=FALSE,sep = opt$separator)
env_table<-env_table[ (rownames(env_table)  %in%  rownames(abund_table)),]
env_table<-env_table[rownames(abund_table),]

meta_table<-read.csv(opt$metaTable,row.names=1,check.names=FALSE,sep = opt$separator)
meta_table<-meta_table[ (rownames(meta_table)  %in%  colnames(abund_table)),]

## Compute cca
logger("Compute CCA")

if(opt$significant_environmental_variable){
  logger("Find significant environmental variable")
  abund_table.adonis <- adonis(abund_table ~ ., data=env_table)
  
  bestEnvVariables<-rownames(abund_table.adonis$aov.tab)[abund_table.adonis$aov.tab$"Pr(>F)"<=opt$significant_co]
  bestEnvVariables<-bestEnvVariables[!is.na(bestEnvVariables)]

  logger(sprintf("Significant environmental varialbe: %s",do.call(paste,c(as.list(bestEnvVariables),sep=" + "))))
  eval(parse(text=paste("sol <- cca(abund_table ~ ",
                        do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),
                        ",data=env_table)",
                        sep="")))

  filename<-paste(opt$prefix,'adonisRes.txt',sep = "_")
  sink(filename)
  print(abund_table.adonis$aov.tab)
  sink()
  
  
} else {
  sol<-cca(abund_table ~ ., data=env_table)
  
}

if(nrow(as.data.frame(sol$CCA$eig)) > 2 ){
  choices<-c(1,2,3)
  nbAxes<-3
} else {
  choices<-c(1,2)
  nbAxes<-2
}
scrs <- scores(sol, scaling = opt$scaling , hill=TRUE, display=c("sp","wa","lc","bp","cn"),choices = choices)

## Extract data
df_sites<-data.frame(scrs$sites,t(as.data.frame(strsplit(rownames(scrs$sites),opt$site_seperator_field))))
meta_table<-meta_table[rownames(scrs$species),]
df_features<-merge(scrs$species, meta_table, by=0, all=T)
rownames(df_features)<-df_features$Row.names
df_features<-df_features[,-1]

## Launch drawing
logger("Draw plots")
p<-list()
k<-1 
for(i in seq(1,nbAxes)){
  for(j in seq(1,nbAxes)){
    df_sites_tmp<-df_sites[,c(i,j,(nbAxes+1):ncol(df_sites))]
    colnames(df_sites_tmp)<-c('x','y',colnames(df_sites[,(nbAxes+1):ncol(df_sites)]))
    
    biplot_tmp<-scrs$biplot[,c(i,j)]
    if (opt$scaling_arrows){
      multiplier <- vegan:::ordiArrowMul(biplot_tmp)
      df_arrows_tmp<- biplot_tmp*multiplier
    } else{
      df_arrows_tmp<- biplot_tmp
    }
    colnames(df_arrows_tmp)<-c("x","y")
    df_arrows_tmp=as.data.frame(df_arrows_tmp)
    
    df_features_tmp<-df_features[,c(i,j,(nbAxes+1):ncol(df_features))]
    colnames(df_features_tmp)<-c('x','y',colnames(df_features[,(nbAxes+1):ncol(df_features)]))

    p[[k]]<-draw(df1=df_sites_tmp,df2=df_arrows_tmp,df3=df_features_tmp,fcb=opt$features_color_by,ssh=opt$sites_shape_by,
                 abscisse=i,ordinate=j,arrowmul=opt$arrowmul)
    k<-k+1
    
  }
}


grid_plot<-list()
k<-1
for(i in seq(1,(nbAxes*nbAxes))){
  grid_plot[[k]]<-p[[k]][[1]]
  k<-k+1
}

## Save files
logger("Save plots")
filename1<-paste(opt$prefix,'featuresPoint.svg',sep = "_")
filename2<-paste(opt$prefix,'featuresText.svg',sep = "_")
filenameTot<-paste(opt$prefix,'featuresPointAllGrid.svg',sep = "_")
dir.create(dirname(filename1), showWarnings = FALSE)

logger(sprintf("Create %s",filename1))
svglite(filename1,width=12,height=9)
print(p[[2]][[1]])
dev.off()

logger(sprintf("Create %s",filename2))
svglite(filename2)
print(p[[2]][[2]])
dev.off()

logger(sprintf("Create %s",filenameTot))
svglite(filenameTot,width=12,height=9)
do.call(grid.arrange,grid_plot)
dev.off()

logger("Analysis done")
end_time <- Sys.time()
logger(sprintf("Runtime: %s sec",end_time - start_time))
