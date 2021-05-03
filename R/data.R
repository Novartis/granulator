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


#'@title Signature matrix for deconvolution of PBMCs in 17 cell types
#'@description A dataset containing the TPM-normalized RNA-seq gene expression
#'values for signature genes of 17 PBMCs.
#'@format A matrix with 1296 rows (genes) and 17 variables (cell types)
#'@usage data(sigMatrix_ABIS_S0)
#'@references Monaco et al. (2019) Cell Reports 26, 1627–1640
#'(\href{https://www.cell.com/cell-reports/pdf/S2211-1247(19)30059-2.pdf}{Cell 
#'Reports})
#'@source \href{https://github.com/giannimonaco/ABIS/tree/master/data}{Github}
"sigMatrix_ABIS_S0"

#'@title Signature matrix for deconvolution of PBMCs in 13 cell types
#'@description A dataset containing the TPM-normalized RNA-seq gene expression
#'values for signature genes of 17 PBMCs.
#'@format A matrix with 1296 rows (genes) and 13 variables (cell types)
#'@usage data(sigMatrix_ABIS_S1)
#'@references Monaco et al. (2019) Cell Reports 26, 1627–1640
#'(\href{https://www.cell.com/cell-reports/pdf/S2211-1247(19)30059-2.pdf}{Cell 
#'Reports})
"sigMatrix_ABIS_S1"

#'@title Signature matrix for deconvolution of PBMCs in 11 cell types
#'@description A dataset containing the TPM-normalized RNA-seq gene expression
#'values for signature genes of 17 PBMCs.
#'@format A matrix with 1296 rows (genes) and 11 variables (cell types)
#'@usage data(sigMatrix_ABIS_S2)
#'@references Monaco et al. (2019) Cell Reports 26, 1627–1640
#'(\href{https://www.cell.com/cell-reports/pdf/S2211-1247(19)30059-2.pdf}{Cell 
#'Reports})
"sigMatrix_ABIS_S2"

#'@title Signature matrix for deconvolution of PBMCs in 9 cell types
#'@description A dataset containing the TPM-normalized RNA-seq gene expression
#'values for signature genes of 17 PBMCs.
#'@format A matrix with 1296 rows (genes) and 9 variables (cell types)
#'@usage data(sigMatrix_ABIS_S3)
#'@references Monaco et al. (2019) Cell Reports 26, 1627–1640
#'(\href{https://www.cell.com/cell-reports/pdf/S2211-1247(19)30059-2.pdf}{Cell 
#'Reports})
"sigMatrix_ABIS_S3"

#'@title PBMCs expression profiles (ABIS dataset)
#'@description Public dataset (GSE107011) containing the TPM-normalized gene
#'expression values from bulk RNAseq of PBMCs of 12 healthy individuals. We 
#'include here only genes selected in the signature matrices.
#'@format A matrix with 1296 rows (genes) and 12 variables (samples)
#'@usage data(bulkRNAseq_ABIS)
#'@references Monaco et al. (2019) Cell Reports 26, 1627–1640
#'(\href{https://www.cell.com/cell-reports/pdf/S2211-1247(19)30059-2.pdf}{Cell
#'Reports})
#'@source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse107011}{
#'GEO}
"bulkRNAseq_ABIS"

#'@title PBMCS true cell type proprotions (ABIS dataset)
#'@description Public dataset (GSE107011)containing the true proportions for 
#'all combinations of cell types (PBMCs) for 12 individuals.
#'@format A matrix with 12 rows (samples) and 24 variables (cell types)
#'@usage data(groundTruth_ABIS)
#'@references Monaco et al. (2019) Cell Reports 26, 1627–1640
#'(\href{https://www.cell.com/cell-reports/pdf/S2211-1247(19)30059-2.pdf}{Cell 
#'Reports})
#'@source \href{https://github.com/giannimonaco/ABIS/tree/master/data}{Github}
"groundTruth_ABIS"
