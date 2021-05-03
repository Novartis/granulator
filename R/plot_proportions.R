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


#'@title Plot estimated cell type proportions
#'
#'@description \code{plot_proportions} plots the estimated cell type 
#'proportions as computed by a given method and signature matrix.
#'
#'@param deconvoluted Output object from function \code{deconvolute}.
#'
#'@param method Character string with name of method to be regressed.
#'
#'@param signature Character string with name of signature to be regressed.
#'
#'@return Plot showing regression of estimated versus measured cell type 
#'coefficients.
#'
#'@import ggplot2
#'@importFrom magrittr %>%
#'@importFrom dplyr select mutate
#'@importFrom tibble as_tibble
#'@importFrom rlang .data
#'
#'@author Vincent Kuettel, Sabina Pfister
#'
#'@examples
#'# load demo PBMCS data
#'load_ABIS()
#'
#'# deconvolute
#'decon <- deconvolute(m = bulkRNAseq_ABIS, 
#'sigMatrix = sigMatrix_ABIS_S0)
#'
#'# plot cell type proportions
#'plot_proportions(deconvoluted = decon, 
#'method = 'svr', signature = 'sig1')
#'
#'@export
plot_proportions <- function(deconvoluted, method = 'svr', signature = 'sig1') {

    # check
    if(!(paste0(method,'_',signature) %in% deconvoluted$combinations$model))
        stop('Method or signature name are not valid names.')

    # data
    dat <- deconvoluted$coefficients[[paste0(method,'_',signature)]]

    # format
    dat <- as_tibble(dat) %>%
        mutate(Sample = rownames(dat)) %>%
        gather("CellType", "Proportions", -.data$Sample)

    # plot
    g <- ggplot(data = dat, aes(y = .data$Proportions, x = .data$Sample, fill = .data$CellType)) +
        geom_bar(position="stack", stat="identity") +
        scale_fill_manual(values = gg_color_hue(unique(dat$CellType),l=50)) +
        theme_bw() +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            plot.title = element_text(size = 12, face = 'bold', hjust = 0.5)
        ) +
        ylim(c(0,100)) +
        labs(
            x = 'samples', 
            y = 'estimated cell type proportions (%)')

    return(g)
}
