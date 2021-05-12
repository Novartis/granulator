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
#'@description \code{model_nnls}  solves the linear model A*coeff=y using
#'a linear algorithm with a non-negative constraint.
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
#'@importFrom nnls nnls
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@noRd
model_nnls = function(y, x, ncores){

    # preprocess
    df <- order_names(y, x)

    # deconvolute
    fit <- apply(df$y, 2, function(z) {
        nnls(A = df$x, b = as.vector(z))})

    # return coefficients
    coeff <- data.frame(lapply(fit, function(z) z$x))
    colnames(coeff) <- colnames(df$y)
    rownames(coeff) <- colnames(df$x)
    coeff <- t(coeff)

    return(list(coeff = coeff))
}
