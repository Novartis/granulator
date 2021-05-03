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


#'@title Load demo PBMCs deconvolution data
#'
#'@description \code{load_ABIS} is used to load a demo dataset for the
#'deconvolution of PBMCs samples from published data under the accession
#'number GSE107011. The dataset consists of the following datasets: 
#'\itemize{
#'\item{bulkRNAseq_ABIS: PBMCs expression profiles}
#'\item{sigMatrix_ABIS_S0: Signature matrix for deconvolution of PBMCs 
#'in 17 cell types}
#'\item{sigMatrix_ABIS_S1: Signature matrix for deconvolution of PBMCs
#'in 13 cell types}
#'\item{sigMatrix_ABIS_S2: Signature matrix for deconvolution of PBMCs
#'in 11 cell types}
#'\item{sigMatrix_ABIS_S3: Signature matrix for deconvolution of PBMCs
#'in 9 cell types}
#'\item{groundTruth_ABIS: PBMCS true cell type proprotions}
#'}
#'
#'@return Returns string confirming successfull loading of the data.
#'
#'@importFrom utils data
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@examples
#'# load data
#'load_ABIS()
#'
#'@export
load_ABIS <- function(){

    # load data
    data("bulkRNAseq_ABIS")
    data("sigMatrix_ABIS_S0")
    data("sigMatrix_ABIS_S1")
    data("sigMatrix_ABIS_S2")
    data("sigMatrix_ABIS_S3")
    data("groundTruth_ABIS")

    return("ABIS dataset was loaded successfully.")
}
