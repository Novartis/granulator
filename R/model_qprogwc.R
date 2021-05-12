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
#'@description \code{model_qprogwc} solves the linear model x*coeff=y using
#'quadratic programming under non-negativity and sum-to-one constraint.
#'
#'@details \code{model_qprogwc} applies the function \code{\link[limSolve]{lsei}}
#'from the \pkg{limSolve} package.
#'
#'@param y Matrix with gene name son rows or columns.
#'
#'@param x Matrix with gene name son rows or columns.

#'@param ncores Number of cores to use for parallel processing.
#'
#'@return Coefficients: a data frame containing estimates for 
#'coefficients for the linear model x*coeff=y.
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@importFrom limSolve lsei
#'
#'@noRd
model_qprogwc <- function(y, x, ncores) {

    # preprocess
    df <- order_names(y, x)

    # deconvolute
    lim_res <- lapply(seq_len(ncol(df$y)), function(z) {
        # parameters to the objective function
        A <- df$x
        B <- as.vector(df$y[,z])

        # equality constraint (weights sum to 1)
        E <- matrix(rep(1, ncol(df$x)), nrow = 1)
        F <- 1

        # inequality constraints (all weights nonnegative)
        G <- diag(1, ncol(df$x))
        H <- rep(0, ncol(df$x))
        return( lsei(A = A, B = B, E = E, F = F, G = G, H = H) )
    })

    # return coefficients
    names(lim_res) <- colnames(df$y)
    lim_with_constraints <- list()
    lim_with_constraints$coeff <- t(data.frame(lapply(lim_res, 
        function(z) z$X)))
    lim_with_constraints <- lapply(lim_with_constraints, 
        function(z) `rownames<-`(z,colnames(df$y)))

    return(list(coeff = lim_with_constraints$coeff))
}
