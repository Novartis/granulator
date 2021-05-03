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


#'@title Plot benchmarking analysis scores
#'
#'@description \code{plot_benchmark} plots the median correlation scores 
#'between estimated and measured cell types across methods and cell types.
#'
#'@param benchmarked List: output object from function \code{benchmarked}.
#'
#'@param metric Character: the metric of evaluation. Options include Pearson
#'Correlation Coefficient ('pcc'), Concordance Correlation Coefficient ('ccc'), 
#'Coefficient of Determination ('adj.r2') and Root Mean Square Error ('rmse') 
#'of the linear regression model.
#'
#'@return Plot showing correlations across algorithms and cell types.
#'
#'@import ggplot2
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
#'bench<- benchmark(deconvoluted = decon, 
#'ground_truth = groundTruth_ABIS)
#'
#'# plot bechmark
#'plot_benchmark(benchmarked = bench, 
#'metric = 'pcc')
#'
#'@export
plot_benchmark <- function(benchmarked, metric = 'pcc'){

    # check
    if(length(metric)>1)
        stop('Only one single metric should be provided; got ', length(metric))

    if (!(metric %in% c('pcc','ccc','adj.r2','rmse')))
        stop('Unkown metric ', metric)

    # plot
    stats <- benchmarked$summary
    data <- stats[,c('method','signature','celltype',metric)]
    colnames(data) <- c('method','signature','celltype','correlation')

    g <- ggplot(data = data, mapping = aes(
        x = .data$celltype, 
        y = .data$correlation, 
        colour = .data$method,
        group = .data$celltype)) +
        stat_summary(
            fun.data=smean_sdl, 
            geom="errorbar", 
            color="black", 
            width=0.4) +
        stat_summary(
            fun=mean, 
            geom='crossbar', 
            color="black", 
            mapping=aes(ymin=.data$..y.., ymax=.data$..y..), 
            width = 0.6, 
            size = 0.4) +
        geom_jitter(
                position=position_jitter(0.2), 
                show.legend = TRUE, size = 2)

    if(metric=="pcc" | metric=="ccc" | metric=="adj.r2")
        g <- g + coord_cartesian(ylim=c(0,1))

    g <- g + facet_wrap(signature~., ncol=1) +
        theme_bw() + 
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.text.x = element_text(size = 10, angle = 90, 
                hjust=0.95,vjust=0.2),
            axis.text.y = element_text(size = 10),
            plot.title = element_text(size = 12, face = 'bold', hjust = 0.5),
            axis.title = element_text(size = 12),
            axis.title.y = element_text(
                margin = margin(t = 0, r = 20, b = 0, l = 0)),
            axis.title.x = element_text(
                margin = margin(t = 20, r = 0, b = 0, l = 0)),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10)
        ) +
        xlab("") +
        ylab(paste0(firstup(
            gsub("pcc","Pearson Correlation Coefficient",
                gsub("ccc","Concordance Correlation Coefficient",
                    gsub("adj.r2","Adjusted R Squared",
                        gsub("rmse","Root Mean Square Error",metric))))))) +
        scale_color_discrete(name = 'Methods')

    return(g)
}
