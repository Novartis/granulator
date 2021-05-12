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


#'@title Deconvolution from bulk RNAseq
#'
#'@description \code{deconvolute} predicts cell type proportions from bulk
#'RNAseq data by applying multiple deconvolution methods.
#'
#'@param m Bulk RNAseq: a genes (rows) by samples (columns) matrix
#'containing transcript-per-million (TPM)-normalized gene expression values.
#'
#'@param sigMatrix Reference profile: a matrix or a named list of matrices.
#'Each signature matrix should be a genes (rows) by cell types (columns) data
#'frame containing TPM-normalized gene expression values of signature genes.
#'
#'@param methods Deconvolution methods: a character vector containing the names
#'of the deconvolution methods to be applied. By default, all methods are run.
#'Functions are either reimplementations of published methods or wrapper
#'functions for published packages:
#'\itemize{
#'\item{ols: ordinary least squares}
#'\item{nnls: non negative least squares
#'regression model. Adapted from Abas et al. (2009)}
#'\item{qprog: quadratic programming without constraints}
#'\item{qprogwc: quadratic programming non-negative and sum-to-one constraints.
#'Adapted from Gong et al. (2015)}
#'\item{dtangle: wrapper for the cell deconvolution function \code{\link[dtangle]{dtangle}}
#'form the package \pkg{dtangle}}
#'\item{rls: robust linear regression. Adapted from Monaco
#'et al. (2019)}
#'\item{svr: support vector regression. Adapted from Newman et al.
#'(2015)}
#'}
#'
#'@param use_cores Number of cores to use for parallel processing
#'
#'@return Returns a list containing two elements: \itemize{
#'\item{coefficients: estimated cell type coefficients}
#'\item{proportions: estimated cell type proportions in percentage}
#'\item{combinations: combination of methods and signatures tested}
#'}
#'
#'@importFrom tidyr crossing everything unite
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@examples
#'# load demo PBMCS data
#'load_ABIS()
#'
#'# generate list of reference profiles to be tested
#'sigMatrix <- list(
#'sig1 = sigMatrix_ABIS_S0, 
#'sig2 = sigMatrix_ABIS_S1)
#'
#'# deconvolute
#'decon <- deconvolute(m = bulkRNAseq_ABIS, 
#'sigMatrix = sigMatrix)
#'
#'@export
deconvolute <- function(m, sigMatrix, 
    methods = get_decon_methods(),
    use_cores = 1){

    # input check
    if (!is.vector(methods))
        stop('methods should be a vector.')

    if (!all(methods %in% get_decon_methods()))
        stop('unknown methods. Available methods: ', paste0(get_decon_methods(),collapse=", ")) 

    if (!is.matrix(m))
        stop('m should be a matrix.')

    if (any(is.na(m)))
        stop("m should not contain missing values.")

    if (any(m<0))
        stop("m should not contain negative values.")

    if (!is.matrix(sigMatrix) & !is.list(sigMatrix))
        stop('sigMatrix should be a matrix or list of multiple matrices.')

    if (is.list(sigMatrix) & length(sigMatrix)==0)
        stop('no valid sigMatrix provided.')

    if (is.list(sigMatrix))
        if (any(unlist(lapply(sigMatrix,function(x) !is.matrix(x)))))
            stop('sigMatrix should be a matrix or list of multiple matrices.')

    if (is.matrix(sigMatrix))
        sigMatrix = list(sigMatrix)

    for(i in seq_along(sigMatrix)){
        if (any(is.na(sigMatrix[[i]])))
            stop("sigMatrix should not contain missing values.")
        if (any(sigMatrix[[i]]<0))
            stop("sigMatrix should not contain missing values.")
        if (nrow(sigMatrix[[i]]) <= ncol(sigMatrix[[i]]))
            stop("in sigMatrix the number of genes should be equal or greater then number of cell types.")
    }

    # fix signature names
    if (is.null(names(sigMatrix))){
        signatures <- vapply(seq_along(sigMatrix), function(i){
        paste0('sig', i)}, FUN.VALUE = character(1))
        names(sigMatrix) <- signatures
    }

    # fix cell type names
    sigMatrix <- lapply(sigMatrix,fix_col_names)

    # create data frame containing methods and signatures combinations
    df <- crossing(method = methods, signature = names(sigMatrix))
    df <- unite(df,'model',remove=FALSE)

    # run methods
    res <- lapply(seq_len(nrow(df)), function(i) {
        message('Running deconvolution method "', df$method[[i]], '" on signature matrix "', df$signature[[i]], '"')
        get(paste0('model_',df$method[[i]]))(y = m, x = sigMatrix[[df$signature[[i]]]], ncores=use_cores)
    })

    # name the results
    res_names <- paste0(df$method, '_', df$signature)
    names(res) <- res_names

    # extract coefficients
    coefficients <- lapply(res_names, function(x) as.data.frame(round(res[[x]]$coeff,2)))
    names(coefficients) <- res_names

    # compute proportions
    proportions <- lapply(coefficients, function(x) as.data.frame(round((100/max(rowSums(x), na.rm=TRUE))*x,2)))

    # output
    res <- list(
        coefficients =  coefficients,
        proportions = proportions, 
        combinations = as.data.frame(df))

    return(res)
}
