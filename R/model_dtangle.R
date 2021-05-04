#-----------------------------------------------
#                 LICENSE
#-----------------------------------------------
# Copyright 2019 Novartis Institutes for BioMedical Research Inc.
# Licensed under the GNU General Public License, Version 3 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# https://www.r-project.org/Licenses/GPL-3
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


#'@title Dtangle
#'
#'@description \code{model_dtangle} computes the cell type proportions from
#'heterogeneous bulk RNA-seq data. \code{model_dtangle} is a wrapper for 
#'the function \code{dtangle} from the package \pkg{limma}.
#'
#'@param sigMatrix Signature matrix: a genes (rows) by cell types (columns)
#'data frame containing TPM-normalized signature genes.
#'
#'@param m Bulk RNA-seq: a genes (rows) by samples (columns) data frame
#'containing TPM-normalized gene expression values.
#'
#'@param ncores Number of cores to use for parallel processing
#'
#'@return Returns a list comprising a data frame containing the estimated
#'coefficients (rows) for each sample (columns).
#'
#'@importFrom dtangle dtangle
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@noRd
model_dtangle <- function(m, sigMatrix, ncores){

    # preprocess
    df <- order_genes(m, sigMatrix)

    # deconvolute
    res <- dtangle(
        Y = log2(as.matrix(t(df$m))+1), 
        references = log2(as.matrix(t(df$sigMatrix))+1),
        marker_method = "ratio")

    # check cell-types without markers
    nas <- names(which(is.na(unlist(lapply(res$markers,function(x) x[[1]])))))

    # remove cell-types without markers
    if(length(nas)>0)
    {
        df <- order_genes(m, df$sigMatrix[,!(colnames(df$sigMatrix) %in% nas)])
        res2 <- dtangle(
            Y = log2(as.matrix(t(df$m))+1), 
            references = log2(as.matrix(t(df$sigMatrix))+1),
            marker_method = "ratio")
        res$estimates[,colnames(res2$estimates)] <- res2$estimates
    }

    # return coefficients
    return(list(coeff = round(100*res$estimates,2)))
}
