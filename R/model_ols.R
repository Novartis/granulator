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


#'@title Ordinary least squares model
#'
#'@description \code{model_ols} infers the cell type proportions from
#'heterogenous bulk RNA-seq data using a simple linear regression model
#'without constraints.
#'
#'@details \code{model_ols} applies the \code{\link[stats]{lsfit}}
#'function for ordinary least squares model regression to estimate cell sub-type
#'proportions from bulk RNA-seq data.
#'
#'@param sigMatrix Signature matrix: a genes (rows) by cell types (columns)
#'data frame containing TPM-normalized signature genes.
#'
#'@param m Bulk RNA-seq: a genes (rows) by samples (columns) data frame
#'containing TPM-normalized gene expression values.
#'
#'@param ncores Number of cores to use for parallel processing.
#'
#'@return Returns a list comprising a data frame containing the estimated
#'coefficients (rows) for each sample (columns).
#'
#'@importFrom stats lsfit
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@noRd
model_ols <- function(m, sigMatrix, ncores)
{
    # preprocess
    df <- order_genes(m, sigMatrix)

    # deconvolute: lm
    fit <- apply(df$m, 2, function(x) {lsfit(
        x = as.matrix(df$sigMatrix), 
        y = as.vector(x), 
        intercept = FALSE)})

    # return coefficients
    coeff <- t(data.frame(lapply(fit, function(x) x$coefficients * 100)))
    rownames(coeff) <- colnames(df$m)

    return(list(coeff = round(coeff,2)))
}
