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


#'@title Plot estimated cell type coefficients against the ground truth
#'
#'@description \code{plot_regress} depicts the measured cell type 
#'proportions (x-axis) vs. the estimated proportions (y-axis).
#'
#'@param benchmarked List: output object from function \code{benchmarked}.
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
#'@importFrom tidyr nest
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
#'# bechmark
#'bench <- benchmark(deconvoluted = decon, 
#'ground_truth = groundTruth_ABIS)
#'
#'# plot regress
#'plot_regress(benchmarked = bench, 
#'method = 'svr', signature = 'sig1')
#'
#'@export
plot_regress <- function(benchmarked, method = 'svr', signature = 'sig1') {

    # check
    if(!(paste0(method,'_',signature) %in% benchmarked$combinations$model))
        stop('method or signature name are not valid names.')

    # data
    data <- benchmarked$data[[paste0(method,'_',signature)]]

    # plot
    g <- ggplot(data = data, aes(y = .data$prediction, x = .data$proportion)) +
        geom_point() +
        geom_smooth(method = "lm", formula = y ~ x, se = TRUE, col = 'blue') +
        facet_wrap(~celltype, nrow = 3, scales = 'free') +
        theme_bw() +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10),
            plot.title = element_text(size = 12, face = 'bold', hjust = 0.5)
        ) +
        labs(
            x = 'measured cell type proportions (%)', 
            y = 'estimated cell type proportions (%)')

    return(g)
}
