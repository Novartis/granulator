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
#'@description \code{model_qprog} solves the linear model x*coeff=y using
#'quadratic programming without constraints.
#'
#'@details \code{model_qprog} applies \code{\link[limSolve]{Solve}} function 
#'from the \pkg{limSolve} package.
#'
#'@param y Matrix with gene name son rows or columns.
#'
#'@param x Matrix with gene name son rows or columns.
#'
#'@param ncores Number of cores to use for parallel processing.
#'
#'@return Coefficients: a data frame containing estimates for 
#'coefficients for the linear model x*coeff=y.
#'
#'@importFrom limSolve lsei
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@importFrom limSolve Solve
#'
#'@noRd
model_qprog <- function(y, x, ncores){

    # preprocess
    df <- order_names(y, x)

    # deconvolute
    lim_res <- lapply(seq_len(ncol(df$y)), function(z){
        # parameters to the objective function
        A <- df$x
        B <- as.vector(df$y[,z])
        return( Solve(A, B) )
    })

    # return coefficients
    names(lim_res) <- colnames(df$y)
    lim_wo_constraints <- list()
    lim_wo_constraints$coeff <- t(data.frame(lapply(lim_res, 
        function(z) z)))
    lim_wo_constraints <- lapply(lim_wo_constraints, 
        function(z) `rownames<-`(z,colnames(df$y)))

    return(list(coeff = lim_wo_constraints$coeff))
}
