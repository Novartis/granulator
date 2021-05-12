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


#'@title Convert raw counts to TPM
#'
#'@description \code{get_TPM} is used to convert raw counts to TPMs,
#'which is the most suitable normalization for deconvolution.
#'
#'@param counts Bulk RNAseq: a genes (rows) by samples (columns) matrix
#'containing gene raw counts.
#'
#'@param effLen Vector of gene lengths.
#'
#'@return Returns a transcript-per-million (TPM)-normalized matrix.
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@examples
#'# get TPMs from raw counts and gene lengths.
#'mat <- round(matrix(rexp(200, rate=.01), ncol=20))
#'len <- round(matrix(rexp(10, rate=.001), ncol=1))+10
#'tpm <- get_TPM(mat,as.vector(len))
#'
#'@export
get_TPM <- function(counts, effLen){

    # input check
    if (!is.vector(effLen))
        stop('effLen should be a vector.')

    if (sum(is.na(effLen))>0)
        stop("effLen should not contain missing values.")

    if (any(effLen<1))
        stop("effLen should not contain values lower than 1.")

    if (!is.matrix(counts))
        stop('counts should be a matrix.')

    if (sum(is.na(counts))>0)
        stop("counts should not contain missing values.")

    if (any(counts<0))
        stop("counts should not contain negative values.")

    # return tpms
    return( as.matrix(apply(counts,2,countToTpm,effLen)) )
}
