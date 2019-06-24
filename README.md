# RECOD
R code for ecological data analysis


## CCA_ggplot.R

This script uses the analysis of variance using distance matrices to find the best set of environmental variables that describe the community structure.

### How to use ?

```
Usage: scripts/CCA_ggplot.R -a [ABUNDTABLE] -m [envTable] -a [metaTable]

Options:
	-a FILE, --abundTable=FILE
		Features abundance table 

	-e FILE, --envTable=FILE
		Environmental data table 

	-m FILE, --metaTable=FILE
		Meta data table 

	-p FILE, --prefix=FILE
		Output prefix for svg file [Default: 'CCA_plot' ]

	--features_color_by=FEATURES_COLOR_BY
		Parameter to color features

	--sites_shape_by=SITES_SHAPE_BY
		Parameter to shape features.
                Must be in this format 'FieldNumber,FieldName'

	--site_seperator_field=SITE_SEPERATOR_FIELD
		Field seperator for sites [Default: '_' ]

	--significant_environmental_variable
		Use adonis to find significant environmental variables and only use these in cca. Use
            all environmmental variables by default

	--significant_co=SIGNIFICANT_CO
		Cutoff for Pr(>F) significance probability value associated with the F Value 
              [Default: 0.01]

	--scaling=SCALING
		Scaling for species and site scores. Either species (2) or site (1) scores are scaled 
            by eigenvalues, and the other set of scores is left unscaled, or with 3 both are scaled 
            symmetrically by square root of eigenvalues. This scaling is know as Hill scaling. The 
            type of scores are 'none', 'sites', 'species', or 'symmetric', which correspond to the 
            values 0, 1, 2, and 3 respectively. For more details help(scores) [Default: '0' ]

	--scaling_arrows
		Scaling arrows using vegan::ordiArrowMul() function. Arrows are scaled by their
              correlation (square root of the column r2) so in that way '1weak' predictors have shorter
              arrows than 'strong' predictors

	--arrowmul=ARROWMUL
		Expand arrows in the graph. Arrows will be scaled automatically to fit the graph 
            {NULL,min,max,optimal} [Default :NULL]

	-s SEPARATOR, --separator=SEPARATOR
		Input files fields separator {comma|tab} [comma]

	-r, --reproductible
		Enable reproductibility of the analysis

	-h, --help
		Show this help message and exit
```


### Citations 

* Inspired by R code for ecological data analysis by Umer Zeeshan Ijaz. 
  * B Torondel, JHJ Ensink, O Gundogdu, UZ Ijaz, J Parkhill, F Abdelahi, V-A Nguyen, S Sudgen, W Gibson, AW Walker, and C Quince. 
Assessment of the influence of intrinsic environmental and geographical factors on the bacterial ecology of pit latrines
Microbial Biotechnology, 9(2):209-223, 2016. DOI:10.1111/1751-7915.12334

## interraction.R 

This scipt allows to compute and visualize intersections of genomic sets as clustered heatmap.

### How to use ?

```
Usage: scripts/intersection.R -f [FILE]

Options:
	-f FILE, --file=FILE
		Input file matrix 

	-p [Y/N], --pairwise=[Y/N]
		Input file is a pairwise matrix [N]

	-c METHOD, --correlation_coeficient=METHOD
		correlation coefficient [kendal]

	-a METHOD, --agglomeration=METHOD
		agglomeration for hclust [complete]

	-d METHOD, --distance=METHOD
		distance matrix computation [euclidian]

	--cluster_number=INTEGER
		Cluster number define to color the dendogram [4]

	-o FILE, --output=FILE
		Ouptut file in png format [pairwise_heatmap.png]

	-h, --help
		Show this help message and exit
```
