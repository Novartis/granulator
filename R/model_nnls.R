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


#'@title Non-negative least squares model
#'
#'@description \code{model_nnls} infers cell type proportions from 
#'bulk RNA-seq data using a linear algorithm with a non-negative constraint.
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
#'@importFrom nnls nnls
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@noRd
model_nnls = function(m, sigMatrix, ncores){

    # preprocess
    df <- order_genes(m, sigMatrix)

    # deconvolute
    fit <- apply(df$m, 2, function(x) {
        nnls(A = as.matrix(df$sigMatrix), b = as.vector(x))})

    # return coefficients
    coeff <- data.frame(lapply(fit, function(x) x[['x']]))*100
    colnames(coeff) <- colnames(df$m)
    rownames(coeff) <- colnames(df$sigMatrix)
    coeff <- t(coeff)

    return(list(coeff = round(coeff,2)))
}
