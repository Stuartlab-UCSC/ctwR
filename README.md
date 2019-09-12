# ctwR
Vignette creating a Cell Type Worksheet from a Seurat Object in R.

### What is a Cell Type Worksheet?
A Cell Type Worksheet is an application designed to ease the burden of manual cell type annotation from single cell
mRNA sequencing experiments. It lets you explore the specificity of markers across clusters and label the clusters
with a cell type annotation.

The web application provides three interactive components for this goal:

1. An editable dot plot visualizing marker specificity and cell type annotation across all clusters.
2. A scatter plot visualizing gene expression across all cells.
3. A table of gene metric rankings per cluster.

Here's a rough visual of the layout of the application, the gene metrics are explored via the table at the bottom.
![Alt text](cell_atlas_layout.png)

This readme servers as instructions for how to make Cell Type Worksheet with the Seurat Package.
### CTW format
The Cell Type Worksheet format is simply a compressed directory with a minimum of 4 tab delimited files:

1. Expression Matrix
2. Cell to Cluster Assignment
3. XY Coordinates
4. Gene Metrics Per Cluster

### Example
For this example we will assume you have already followed the Seurat pipeline down past the
dimensionality reduction stage and the clustering stage.

The first step of the process is to make a worksheet directory on your machine. Your name for the
worksheet is encoded in the directory name.

```
mkdir /path/to/worksheet-name
```
You'll be writing necessary data to the worksheet directory.

The next step is to determine what genes you would like to explore in the
Cell Type Worksheet. We suggest using Seurat's FindAllMarkers function. In an R session:
```R
library("Seurat")

# If working with a lot of cells consider a multiprocess strategy for
# calculating the differentials.
#library("future")
#options(future.globals.maxSize = 10 * 1024 ^ 3)
#plan(strategy = "multiprocess", workers = availableCores())

# Read in your seurat object
sboj <- readRDS("/path/to/your/seuratv3.rds")

# This returns makers based on significance of the differential
markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = -Inf, logfc.threshold = -Inf, return.thresh = Inf)
```

The next step is to gather the values that will fill your dotplot. The easiest way
to do this is to use Seurat's DotPlot function.

```R
dotplot = DotPlot(sobj, 
     features = markers, cols.use = c("blue","red", "green"),
     x.lab.rot = T, plot.legend = F,
     dot.scale = 8, do.return = T
     )

markers = dotplot$data

# Rename the columns of the table to the correct format.
colnames(markers) <- c("avg.exp", "pct.exp", "gene", "cluster", "avg.exp.scaled")
# Reorder the columns of the table to the correct format.
markers <- markers[c("gene", "avg.exp.scaled", "pct.exp", "avg.exp", "cluster")]
```

Now we have all the data we need inside the seurat object and in the markers dataframe.

Our next step is to write them all to our worksheet directory.

```R
# Wrapper to access UMAP xy coordinates from Seurat
getUMAP <- function(object){
  object@reductions$umap@cell.embeddings[,1:2]
}
# Write xy coordinates to file.
write.table(getUMAP(sobj),"/path/to/worksheet-name/xys.tsv", sep="\t")

# Write cell to cluster assignment to file.
write.table(Idents(sobj),"/path/to/worksheet-name/clustering.tsv", sep="\t")

# Accessing expression matrix.
exp = GetAssayData(object = fetalCombined, slot = "data")
# Filter to markers found.
exp = exp[rownames(exp) %in% as.character(markers),]
# Write the expression matrix to file.
exp <- data.frame(exp)
write.table(exp, "/path/to/worksheet-name/exp.tsv", sep="\t", row.names=T)
```

Now you have written out four key types of information for the Cell Type Worksheet:
1. Expression Matrix
2. Cell to Cluster Assignment
3. XY Coordinates of dimensionality reduction algorithm.
4. Gene Metrics Per Cluster

Now it's time to tar and gzip your worksheet directory.
```
tar -cvzf worksheet-name.ctw.tgz /path/to/worksheet-name
```
Your worksheet-name.ctw.tgz file is ready to upload to the Cell Type Workbench.
