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


#'@title Quadratic programming without constraints
#'
#'@description \code{model_qprog} infers the cell type proportions from 
#'heterogenous bulk RNA-seq data using quadratic programming without constraints.
#'
#'@details \code{model_qprog} applies \code{\link[limSolve]{Solve}} function 
#'from the \pkg{limSolve} package to deconvolve bulk RNA-seq data by quadratic 
#'programming. The function \code{\link[limSolve]{Solve}} requires a signature 
#'matrix and bulk RNA-seq data as input.
#'
#'@param m Bulk RNAseq: a genes (rows) by samples (columns) data frame
#'containing transcript-per-million (TPM)-normalized gene expression values
#'
#'@param sigMatrix Signature matrix: A data frame or a list of data frames.
#'Each signature matrix should be a genes (rows) by cell types (columns)
#'data frame containing TPM-normalized gene expression values of signature 
#'genes.
#'
#'@param ncores Number of cores to use for parallel processing.
#'
#'@return Returns a list comprising a data frame containing the estimated
#'coefficients (rows) for each sample (columns).
#'
#'@importFrom limSolve lsei
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@importFrom limSolve Solve
#'
#'@noRd
model_qprog <- function(m, sigMatrix, ncores){

    # preprocess
    df <- order_genes(m, sigMatrix)

    # deconvolute
    lim_res <- lapply(seq_len(ncol(df$m)), function(i){
        # parameters to the objective function
        A <- as.matrix(df$sigMatrix)
        B <- as.vector(df$m[,i])
        return( Solve(A, B) )
    })

    # return coefficients
    names(lim_res) <- colnames(df$m)
    lim_wo_constraints <- list()
    lim_wo_constraints$coeff <- t(data.frame(lapply(lim_res, 
        function(x) x*100)))
    lim_wo_constraints <- lapply(lim_wo_constraints, 
        function(x) `rownames<-`(x,colnames(df$m)))

    return(list(coeff = round(lim_wo_constraints$coeff,2)))
}
