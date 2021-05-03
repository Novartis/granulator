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


#'@title Quadratic program with constraints
#'
#'@description \code{model_qprogwc} infers the cell type proportions from 
#'heterogenous bulk RNA-seq data using quadratic programming under non-negativity 
#'and sum-to-one constraint.
#'
#'@details \code{model_qprogwc} applies the function \code{\link[limSolve]{lsei}}
#'from the \pkg{limSolve} package to deconvolve bulk RNA-seq data by quadratic
#'programming. The function \code{\link[limSolve]{lsei}} requires a signature
#'matrix and bulk RNA-seq data as input and imposes non-negativity and 
#'sum-to-one constraints.
#'
#'@param m Bulk RNAseq: a genes (rows) by samples (columns) data frame
#'containing transcript-per-million (TPM)-normalized gene expression values.
#'
#'@param sigMatrix Signature matrix: A data frame or a list of data frames.
#'Each signature matrix should be a genes (rows) by cell types (columns) data
#'frame containing TPM-normalized gene expression values of signature genes.

#'@param ncores Number of cores to use for parallel processing.
#'
#'@return Returns a list comprising a data frame containing the estimated
#'coefficients (rows) for each sample (columns).
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@importFrom limSolve lsei
#'
#'@noRd
model_qprogwc <- function(m, sigMatrix, ncores) {

    # preprocess
    df <- order_genes(m, sigMatrix)

    # deconvolute
    lim_res <- lapply(seq_len(ncol(df$m)), function(i) {
        # parameters to the objective function
        A <- as.matrix(df$sigMatrix)
        B <- as.vector(df$m[,i])

        # equality constraint (weights sum to 1)
        E <- matrix(rep(1, ncol(df$sigMatrix)), nrow = 1)
        F <- 1

        # inequality constraints (all weights nonnegative)
        G <- diag(1, ncol(df$sigMatrix))
        H <- rep(0, ncol(df$sigMatrix))
        return( lsei(A = A, B = B, E = E, F = F, G = G, H = H) )
    })

    # return coefficients
    names(lim_res) <- colnames(df$m)
    lim_with_constraints <- list()
    lim_with_constraints$coeff <- t(data.frame(lapply(lim_res, 
        function(x) x$X*100)))
    lim_with_constraints <- lapply(lim_with_constraints, 
        function(x) `rownames<-`(x,colnames(df$m)))

    return(list(coeff = round(lim_with_constraints$coeff,2)))
}
