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
#'@description \code{model_ols}  solves the linear model A*coeff=y using
#'a simple linear regression model without constraints.
#'
#'@details \code{model_ols} applies the \code{\link[stats]{lsfit}}
#'function for ordinary least squares model regression.
#'
#'@param y Matrix with gene name son rows or columns.
#'
#'@param x Matrix with gene name son rows or columns.
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
model_ols <- function(y, x, ncores)
{
    # preprocess
    df <- order_names(y, x)

    # deconvolute: lm
    fit <- apply(df$y, 2, function(z) {lsfit(
        x = df$x, 
        y = as.vector(z), 
        intercept = FALSE)})

    # return coefficients
    coeff <- t(data.frame(lapply(fit, function(z) z$coefficients)))
    rownames(coeff) <- colnames(df$y)

    return(list(coeff = coeff))
}
