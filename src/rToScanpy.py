
# Seurat 3 to AnnData conversion in python - Henry
# Motivation - run output in python, able to convert to AnnData or Cell Type Worksheet .tsv files
# Seurat 3 object is an S4 R object
# it has "slots" which you can access with the @ accessor
# In those slots are objects which are either accessible with either @ or $
# or just e.g. matrices or DataFrames or lists or vectors themselves.
# Relevant Seurat3 structure for our purposes:
# @assays
#      $RNA
#          @counts - dgCMatrix, which corresponds to python csc_matrix
#             $x - values of non-empty cells
#             $i - indices of non-empty cells
#             $p - pointers of non-empty cells
#             $Dimnames - row/colnames (not a part of the csc_matrix object);
#                     this is an R list of two vectors
# @meta.data - metadata, rows are samples; R DataFrame
# @reductions - dimensionality reductions; R matrices
#     $umap - uniform manifold and projection
#         @cell.embeddings - x,y
#     $pca - principal components analysis
#         @cell.embeddings - all principal components
# AnnData structure for our purposes, access slots with .
# X = data in full or sparse format, obs x var dimensions
# obs = metadata, with rows as samples
# var = gene names from @assays$RNA@counts$Dimnames[[1]]
#     could also be a dataframe with gene info
# obsm = layouts (x, y)
#     ["pca"] = PCA
#     ["umap"] = UMAP
# Cell Type Worksheet files
# clustering.tsv - per cell cluster assignment
# columns: cell, cluster
# clusters.tsv - per cluster info
# columns: cluster, cell_count, cell_type (optional)
# xys.tsv - x, y scatter plot coordinates
# columns: cell, x, y
# exp.tsv
# columns: gene, cell1, cell2, ..., celln
​
​
​
import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
base = importr("base")
dollar = base.__dict__["$"]
at = base.__dict__["@"]
​
def rpy2DataFrameToPandasDataFrame(df):
    # convert an "R/rpy2 DataFrame" to a pandas dataframe.
    # first extract column and row IDs from dataframe.
    cols = list(base.colnames(df))
    rows = list(base.rownames(df))
    data = pd.DataFrame(columns=cols, index=rows)
    # now put the data into the pandas dataframe
    for j in range(len(cols)):
        data.loc[:, cols[j]] = list(dollar(df, cols[j]))
    return data
​
from scipy.sparse import csc_matrix
bracket = base.__dict__["["]
doubleBracket = base.__dict__["[["]
​
def rpy2dgCMatrixToScipySparseCSCMatrix(dgCM):
    # R sparse matrix dgCMatrix to pythonic sparse matrix csc_matrix.
    # csc_matrix((data, indices, indptr), [shape=(M, N)])
    # note - shape is input as a list not a tuple or you will get an error
    # dgCMatrix matches CSC matrix in python, with indptr = p and indices = i
    # first get the column and row IDs from dgCMatrix.
    rows = list(doubleBracket(at(dgCM, "Dimnames"), 1))
    cols = list(doubleBracket(at(dgCM, "Dimnames"), 2))
    # then get the components of the csc_matrix from the dgCMatrix.
    dgCMData = list(at(dgCM, "x"))
    dgCMIndices = list(at(dgCM, "i"))
    dgCMIndptr = list(at(dgCM, "p"))
    # note that the csc_matrix has no slots for row/column names.
    cscM = csc_matrix((dgCMData, dgCMIndices, dgCMIndptr), [len(rows), len(cols)])
    return cscM, rows, cols
​
from anndata import AnnData
​
def rpy2Seurat3ToAnnData(s3):
    # Convert Seurat3 to AnnData
    # TODO: detect UMAP, tSNE, PCA
    # TODO: scaled data?
    # first, convert the data
    expression, exprRows, exprCols = rpy2dgCMatrixToScipySparseCSCMatrix(at(dollar(at(s3, "assays"), "RNA"), "counts"))
    # get the metadata
    meta = rpy2DataFrameToPandasDataFrame(at(s3, "meta.data"))
    # create the AnnData object
    ad = AnnData(X=expression.transpose(), obs=meta, var=exprRows)
    # then get the dimensionality reductions
    ad.obsm["umap"] = np.array(at(dollar(at(s3, "reductions"), "umap"), "cell.embeddings"))
    ad.obsm["pca"] = np.array(at(dollar(at(s3, "reductions"), "pca"), "cell.embeddings"))
    return ad
​
import math
from collections import Counter
import time
​
def rpy2Seurat3ToCTWInput(s3, outDir, transform, metaClusterField,
                          chunkSize=5000,
                          metaCellTypeField = None):
    # Convert Seurat3 to Cell Type Worksheet Input .tsv files
    # TODO: graceful failure mode if transform not available
    # TODO: alternatively, user supplies transform name and we are agnostic
    # First, make the directory if it doesn't already exist.
    if not os.path.isdir(outDir):
        os.mkdir(outDir)
        print("made directory " + outDir)
    startTime = time.time()
    # Convert the data to a pythonic format.
    expression, exprRows, exprCols = rpy2dgCMatrixToScipySparseCSCMatrix(at(dollar(at(s3, "assays"), "RNA"), "counts"))
    print("Done converting dgCMatrix to csr_matrix, elapsed time " + str(time.time() - startTime))
    meta = rpy2DataFrameToPandasDataFrame(at(s3, "meta.data"))
    # Output xys.tsv
    # Get the input data.
    if transform == "umap":
        xys = np.array(at(dollar(at(s3, "reductions"), "umap"), "cell.embeddings"))
    elif transform == "pca":
        # PCA likely has many PCs, we only want the first two.
        xys = np.array(at(dollar(at(s3, "reductions"), "pca"), "cell.embeddings"))[:, 0:2]
    elif transform == "tsne":
        xys = np.array(at(dollar(at(s3, "reductions"), "tsne"), "cell.embeddings"))
    # Print to file.
    xys = pd.DataFrame(xys, columns = ["x", "y"], index = meta.index)
    xys.to_csv(path_or_buf = outDir + "xys.tsv", sep = "\t", header = True,
               index = True, index_label = "cell")
    print("Done outputting xys, elapsed time " + str(time.time() - startTime))
    # Output clustering.tsv
    # Which is the cell to cluster correspondence.
    clustering = pd.DataFrame({"cluster": meta.loc[:, metaClusterField]}, columns = ["cluster"], index = meta.index)
    clustering.to_csv(path_or_buf = outDir + "clustering.tsv", header = True, index = True, index_label = "cell", sep = "\t")
    # Output clusters.tsv
    # Which is the cluster metadata file. Just two fields - cell count and cell type (which is optional).
    # Get the unique clusters and count them.
    clusterCounter = Counter(meta.loc[:, metaClusterField])
    clusterNames = list(set(meta.loc[:, metaClusterField]))
    clusterCounts = [clusterCounter[clusterNames[i]] for i in range(len(clusterNames))]
    
    # If the user specifies a cell type field, get that.
    if metaCellTypeField == None:
        clusters = pd.DataFrame({"cell_count": clusterCounts}, index = clusterNames, columns = ["cell_count"])
    else:
        clusterCellTypes = clusterNames
        for clusterIndex, cluster in clusterNames:
            clusterCellTypes[clusterIndex] = list(set(meta.loc[meta.loc[:, metaClusterField] == cluster, metaClusterField]))[0]
        clusters = pd.DataFrame({"cell_count": clusterCounts, "cell_type": clusterCellTypes}, columns = ["cell_count", "cell_type"])
    # Output to file.
    clusters.to_csv(path_or_buf = outDir + "clusters.tsv", header=True, index=True, index_label = "cluster", sep = "\t")
    print("Done outputting cluster info, elapsed time " + str(time.time() - startTime))
    # Output expression as chunks
    # Calculate how many chunks will be needed. A smaller chunk size uses less RAM.
    nChunk = math.ceil(len(exprRows)/chunkSize)
    for thisChunk in range(nChunk):
        # Which genes are in this chunk
        thisRange = range(thisChunk * chunkSize, (thisChunk + 1) * chunkSize)
        thisIndex = exprRows[(thisChunk * chunkSize):(((thisChunk + 1) * chunkSize) - 0)] # subset of genes, in chunk_size chunks
        if thisChunk == (nChunk - 1): # last chunk - will be smaller than or equal to chunkSize
            thisExpression = pd.DataFrame(expression[(thisChunk * chunkSize):, ].todense(), index = thisIndex, columns = meta.index)
        else:
            thisExpression = pd.DataFrame(expression[thisRange, ].todense(), index = thisIndex, columns = meta.index)
        if thisChunk == 0: # for the first chunk, output the header.
            thisExpression.to_csv(path_or_buf = outDir + "exp" + str(thisChunk) + ".tsv",
                                  sep = "\t", header = True, index = True, index_label = "gene")
        else: # otherwise don't so we don't have to trim it later.
            thisExpression.to_csv(path_or_buf = outDir + "exp" + str(thisChunk) + ".tsv",
                                  sep = "\t", header = False, index = True)
        print("done exporting chunk " + str(thisChunk))
    print("Done outputting exp.tsv, elapsed time " + str(time.time() - startTime))
    # Combine the chunks into one file, then delete the chunks.
    with open((outDir + 'exp.tsv'),'wb') as wfd:
        myCatFiles = []
        for thisChunk in range(nChunk):
            myCatFiles = myCatFiles + [(outDir + "exp" + str(thisChunk) + ".tsv")]
        for f in myCatFiles:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)
        for thisFile in myCatFiles:
            os.remove(thisFile)
    print("All done, elapsed time " + str(time.time() - startTime))
