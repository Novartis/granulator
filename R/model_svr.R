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


#'@title Support vector regression
#'
#'@description \code{model_svr} solves the linear model x*coeff=y using
#'support vector regression.
#'
#'@details \code{model_svr} applies the \code{\link[e1071]{svm}} function.
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
#'@import parallel
#'@importFrom e1071 svm
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@noRd
model_svr <- function(y, x, ncores){

    # preprocess
    df <- order_names(y, x)

    # hyperparameter
    nu <- c(0.25, 0.5, 0.75)

    # deconvolute
    res <- mclapply(seq_len(ncol(df$y)), function(z) {
        lapply(nu, function(i) {
            suppressMessages(
                suppressWarnings(
                    svm(x = df$x,
                        y = as.vector(df$y[,z]), 
                        nu = i, 
                        type = 'nu-regression',
                        scale = TRUE, 
                        kernel = 'linear')
            ))
        })
    }, mc.cores = ncores)

    # rmse for each model
    rmse <-  lapply(res, function(z) { 
        lapply(z, function(i){
            sqrt(sum(i$residuals ^ 2)/length(i$residuals))
        })
    })

    # find minimum rmse for each sample
    rmse <- unlist(lapply(rmse, function(z) which.min(unlist(z))))

    # extract models with min rmse
    mod <- lapply(seq_along(res), function(z) {res[[z]][[rmse[z]]] })

    # return coefficients
    coefs <- lapply(mod, function(z) {as.numeric(t(z$coefs) %*% z$SV)})
    coefs <- as.data.frame(t(data.frame(coefs)))
    rownames(coefs) <- colnames(df$y)
    colnames(coefs) <- colnames(df$x)

    return(list(coeff = coefs))
}
