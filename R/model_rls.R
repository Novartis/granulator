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
#'@description \code{model_rls}  solves the linear model A*coeff=y using
#'a robust linear model.
#'
#'@details \code{model_rls} applies the \code{\link[MASS]{rlm}} function 
#'for multivariate robust linear regression.
#'
#'@param y Matrix with gene name son rows or columns.
#'
#'@param x Matrix with gene name son rows or columns.
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
model_rls <- function(y, x, ncores){

    # preprocess
    df <- order_names(y, x)

    # deconvolute
    rlm <- apply(df$y, 2, function(z){
        suppressMessages(
                suppressWarnings(
                    rlm(df$x, as.vector(z),
                    maxit = 1000, psi = psi.huber)))
    })

    # return coefficients
    coeff <- t(data.frame(lapply(rlm, function(z) z$coefficients)))
    rownames(coeff) <- colnames(df$y)

    return(list(coeff = coeff))
}
