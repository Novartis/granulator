#-----------------------------------------------
#                 LICENSE
#-----------------------------------------------
# Copyright 2018 Novartis Institutes for BioMedical Research Inc.
# Licensed under the GNU General Public License, Version 3 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# https://www.r-project.org/Licenses/GPL-3
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


#'@title Robust weighted least squares
#'
#'@description \code{model_rls} infers the cell type proportions from 
#'heterogenous bulk RNA-seq data using a robust linear model.
#'
#'@details \code{model_rls} applies the \code{\link[MASS]{rlm}} function 
#'for multivariate robust linear regression to estimate cell sub-type proportions 
#'from bulk RNA-seq data. In latent linear modeling, coefficients correspond to 
#'the proportions of a cell sub-types.
#'
#'@param m Bulk RNAseq: a genes (rows) by samples (columns) data frame
#'containing transcript-per-million (TPM)-normalized gene expression
#'values
#'
#'@param sigMatrix Signature matrix: A data frame or a list of data frames.
#'Each signature matrix should be a genes (rows) by cell types (columns)
#'data frame containing TPM-normalized gene expression values of signature
#'genes.
#'
#'@param ncores Number of cores to use for parallel processing
#'
#'@return Coefficients: a data frame containing estimates for cell sub-type
#'proportions (columns) for individuals (rows).
#'
#'@importFrom MASS rlm psi.huber
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@noRd
model_rls <- function(m, sigMatrix, ncores){

    # preprocess
    df <- order_genes(m, sigMatrix)

    # deconvolute
    rlm <- apply(df$m, 2, function(x){
        suppressMessages(
                suppressWarnings(
                    rlm(as.matrix(df$sigMatrix), as.vector(x),
                    maxit = 1000, psi = psi.huber)))
    })

    # return coefficients
    coeff <- t(data.frame(lapply(rlm, function(x) x$coefficients*100)))
    rownames(coeff) <- colnames(df$m)

    return(list(coeff = round(coeff,2)))
}
