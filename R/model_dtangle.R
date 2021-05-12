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


#'@title Dtangle
#'
#'@description \code{model_dtangle}  solves the linear model x*coeff=y using
#'the function \code{dtangle} from the package \pkg{limma}.
#'
#'@param sigMatrix Signature matrix: a genes (rows) by cell types (columns)
#'data frame containing TPM-normalized signature genes.
#'
#'@param y Matrix with gene name son rows or columns.
#'
#'@param x Matrix with gene name son rows or columns.
#'
#'@return Coefficients: a data frame containing estimates for 
#'coefficients for the linear model x*coeff=y.
#'
#'@importFrom dtangle dtangle
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@noRd
model_dtangle <- function(y, x, ncores){

    # preprocess
    df <- order_names(y, x)

    # deconvolute
    res <- dtangle(
        Y = log2(t(df$y)+1), 
        references = log2(t(df$x)+1),
        marker_method = "ratio")

    # check cell-types without markers
    nas <- names(which(is.na(unlist(lapply(res$markers,function(x) x[[1]])))))

    # remove cell-types without markers
    if(length(nas)>0)
    {
        df <- order_names(df$y, df$x[,!(colnames(df$x) %in% nas)])
        res2 <- dtangle(
            Y = log2(t(df$y)+1), 
            references = log2(t(df$x)+1),
            marker_method = "ratio")
        res$estimates[,colnames(res2$estimates)] <- res2$estimates
    }

    # return coefficients
    return(list(coeff = res$estimates))
}
