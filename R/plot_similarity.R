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


#'@title Plot reference profile similariy matrix
#'
#'@description \code{plot_similarity} plots cell type similarity matrix by
#'computing the Kendall rank correlations between cell type expression profiles.
#'Kendall rank correlation is used to test the similarities in the ordering of data 
#'when it is ranked by quantities, and  provides a less inflated measure of accuracy 
#'than Pearson correlation by accounting for ties in the data.
#'
#'@param sigMatrix Signature matrix: a data frame or a named list of data frames.
#'Each signature matrix should be a genes (rows) by cell types (columns) data
#'frame containing TPM-normalized gene expression values of signature genes.
#'
#'@return Plot showing the Kendall rank correlations similariy matrix.
#'
#'@import ggplot2
#'@importFrom magrittr %>%
#'@importFrom stats cor 
#'@importFrom dplyr select filter arrange mutate pull
#'@importFrom tidyr gather
#'@importFrom tibble as_tibble
#'@importFrom pheatmap pheatmap
#'@importFrom ggplotify as.ggplot
#'@importFrom grDevices colorRampPalette
#'@importFrom cowplot plot_grid
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@examples
#'# load demo PBMCS data
#'load_ABIS()
#'
#'# generate list of reference profiles to be tested
#'sigMatrix <- list(sig1 = sigMatrix_ABIS_S0, 
#'sig2 = sigMatrix_ABIS_S2)
#'
#'# plot similarity
#'plot_similarity(sigMatrix = sigMatrix)
#'
#'@export
plot_similarity <- function(sigMatrix){

    # input check
    if (!is.matrix(sigMatrix) & !is.list(sigMatrix))
        stop('Signature matrix should be a matrix or list of multiple matrices; got ', class(sigMatrix))

    if (is.matrix(sigMatrix))
        sigMatrix = list(sig = sigMatrix)

    # fix cell type names
    sigMatrix <- lapply(sigMatrix,fix_col_names)

    # plot
    gList <- list()
    for (i in seq_along(sigMatrix)){

        name <- names(sigMatrix)[i]
        reference <- log2(as.matrix(t(sigMatrix[[i]]))+1)

        breaksList = seq(0, 1, by = 0.05)

        g <- as.ggplot(pheatmap(mat = cor(t(reference), method = "kendall"),
            clustering_method = "ward.D2",
            color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)),
            breaks = breaksList,
            fontsize_row = 3+7/length(sigMatrix), 
            fontsize_col = 3+7/length(sigMatrix),
            legend = TRUE,
            treeheight_row = 0,
            treeheight_col = 0,
            silent=TRUE)) +
            ggtitle(paste0(name,' (k=',round(kappa(t(reference))),')')) +
            coord_equal() +
            theme(plot.title = element_text(size = 12, face = 'bold', hjust = 0.5, margin=margin(0,0,5,0)), aspect.ratio = 0.8)

        gList <- c(gList,list(g))

    }

    p <- plot_grid(plotlist = gList)

    return(p)
}
