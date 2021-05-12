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


#'@title Deconvolution methods acronyms
#'
#'@description \code{get_decon_methods} returns supported deconvolution 
#'methods acronyms.
#'
#'@return vector containing the acronyms of deconvolution methods.
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@examples
#'# get available deconvolution methods
#'get_decon_methods()
#'
#'@export
get_decon_methods <- function(){

    # return available deconvolution methods
    return(c('ols','nnls','qprog','qprogwc','rls','svr','dtangle'))
}
