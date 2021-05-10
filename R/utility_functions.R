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

# Collections of functions required to preprocess the input data for various
# deconvolution methods.


#'@title Fix column names
#'@description \code{\link{fix_names}} replace underscore from column names.
#'@param matrix Matrix with column names to be fixed.
#'@return matrix: a matrix or data frame with column names without underscore.
#'@author Vincent Kuettel, Sabina Pfister
#'@noRd
fix_col_names <- function(matrix){
    colnames(matrix) <- gsub("_",".",colnames(matrix))
    return(matrix)
}

#'@title Extaction of overlapping  between two matrices
#'@description \code{\link{order_names}} finds intersecting samples between
#'two matrices.names
#'@param y Matrix or data.frame.
#'@param x Matrix or data.frame.
#'@return Returns a list of matrices filtered by overlapping names
#'@author Vincent Kuettel, Sabina Pfister
#'@noRd
order_names <- function(y, x){
    if(length(intersect(rownames(y), rownames(x)))>0){
        names <- intersect(rownames(y), rownames(x))
        y <- as.matrix(y[names,,drop=FALSE])
        y <- as.matrix(y[order(rownames(y)), order(colnames(y)),drop=FALSE])
        x <- as.matrix(x[names,])
        x <- as.matrix(x[order(rownames(x)), order(colnames(x))])
    } else if(length(intersect(colnames(y), colnames(x)))>0){
        names <- intersect(colnames(y), colnames(x))
        y <- as.matrix(y[,names,drop=FALSE])
        y <- as.matrix(y[order(rownames(y)), order(colnames(y)),drop=FALSE])
        x <- as.matrix(x[,names])
        x <- as.matrix(x[order(rownames(x)), order(colnames(x))])
    } else
        stop("No overlapping names found.") 
    return(list(y=as.matrix(y), x=as.matrix(x)))
}

#'@title get supported deconvolution methods acronyms
#'@description \code{get_decon_methods} returns vector of deconvolution methods
#'supported in this package.
#'@return vector containing the acronyms of deconvolution methods
#'@noRd
get_decon_methods <- function(){
    return(c('ols','nnls','qprog','qprogwc','rls','svr','dtangle'))
}

#'@title extract sample (donor) names from row names of ABIS data
#'@description \code{\link{subject_names}} processes cell type labels in the 
#'default input data.
#'@param df data frame
#'@return data frame with unique sample names
#'@noRd
subject_names <- function(df){
    sub_names <- gsub('\\_.*', '', colnames(df))
    sub_names <- unlist(vapply(sub_names, 
                    function(x) if (substr(x, 1, 1) == 'X')
                        {gsub('^.', '', x)}
                        else {x = x}, character(1)))
    sub_names <- sub_names[order(sub_names)]
    colnames(df) <- sub_names
    return(df)
}

#'@title Capital first letter
#'@description \code{\link{firstup}} returns string with capital first letter
#'@param string string of characters
#'@return string of characters with capital first letter
#'@noRd
firstup <- function(x)
{
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    return( x )
}

#'@title emulate ggplot2 color palette
#'@description \code{\link{gg_color_hue}} returns ggplot2 colors mapped to values
#'@param vec vector of instances to be mapped to colors
#'@return string vector of ggplot2 colors mapped to values
#'@importFrom grDevices hcl
#'@noRd
gg_color_hue <- function(vec, l=65)
{
    n <- length(unique(vec))
    hues = seq(15, 375, length = n + 1)
    res <- hcl(h = hues, l = l, c = 100)[seq_len(n)]
    names(res) <- unique(vec)[order(unique(vec))]
    return ( res )
}

#'@title Mean and standard deviation
#'@description \code{\link{smean_sdl}} returns the mean plus or minus the standard deviation.
#'Reimplementation of the function \code{smean.sdl} from the package \pkg{Hmisc}.
#'@param x numeric vector
#'@param na.rm Boolean: indicate if NAs should be removed (default: TRUE)
#'@return vector containg mean, lower sd, and upper sd
#'@noRd
smean_sdl <- function(x, na.rm=TRUE)
{
    if(na.rm)
        x <- x[!is.na(x)]
        
    n <- length(x)
    if(n == 0)
        return(c(Mean=NA, Lower=NA, Upper=NA))

    xbar <- sum(x)/n
    sd <- sqrt(sum((x - xbar)^2)/(n-1))
    out <- c(y=xbar, ymin=xbar - sd, ymax=xbar + sd)
    return( out )
}

#'@title Counts to TPMs
#'@description Convert counts to TPMs.
#'@param counts vector of counts
#'@param effLen gene length
#'@return vector of TPMs
#'@noRd
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  ifelse(counts == 0, 0, exp(rate - denom + log(1e6)))
}
