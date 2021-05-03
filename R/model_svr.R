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
#'@description \code{model_svr} deconvolves the estimated cell sub-type
#'proportions from bulk RNA-seq data.
#'
#'@details \code{model_svr} applies the \code{\link[e1071]{svm}} function.
#'
#'@param m Bulk RNAseq: a genes (rows) by samples (columns) data frame
#'containing transcript-per-million (TPM)-normalized gene expression
#'values.
#'
#'@param sigMatrix Signature matrix: A data frame or a list of data frames.
#'Each signature matrix should be a genes (rows) by cell types (columns) data 
#'frame containing TPM-normalized gene expression values of signature genes.
#'
#'@param ncores Number of cores to use for parallel processing.
#'
#'@return Coefficients: a data frame containing estimates for cell sub-type
#'proportions (columns) for individuals (rows).
#'
#'@import parallel
#'@importFrom e1071 svm
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@noRd
model_svr <- function(m, sigMatrix, ncores){

    # preprocess
    df <- order_genes(m, sigMatrix)

    # hyperparameter
    nu <- c(0.25, 0.5, 0.75)

    # deconvolute
    res <- mclapply(seq_len(ncol(df$m)), function(x) {
        lapply(nu, function(i) {
            suppressMessages(
                suppressWarnings(
                    svm(x = as.matrix(df$sigMatrix),
                        y = as.vector(df$m[,x]), 
                        nu = i, 
                        type = 'nu-regression',
                        scale = TRUE, 
                        kernel = 'linear')
            ))
        })
    }, mc.cores = ncores)

    # rmse for each model
    rmse <-  lapply(res, function(x) { 
        lapply(x, function(i){
            sqrt(sum(i$residuals ^ 2)/length(i$residuals))
        })
    })

    # find minimum rmse for each sample
    rmse <- unlist(lapply(rmse, function(x) which.min(unlist(x))))

    # extract models with min rmse
    mod <- lapply(seq_along(res), function(i) {res[[i]][[rmse[i]]] })

    # return coefficients
    coefs <- lapply(mod, function(x) {as.numeric(t(x$coefs) %*% x$SV)})
    coefs <- as.data.frame(t(data.frame(coefs)))*100
    rownames(coefs) <- colnames(df$m)
    colnames(coefs) <- colnames(df$sigMatrix)
    coefs <- apply(coefs, 2, function(x) {ifelse(x < 0, 0, x)})

    return(list(coeff = round(coefs,2)))
}
