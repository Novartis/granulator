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


#'@title Plot of correlations between deconvolution methods
#'
#'@description \code{plot_correlate} is used to visualize the results 
#'obtained by \code{correlation_analysis}.
#'
#'@details \code{plot_correlate} plots the correlation of cell type 
#'proportions across methods in form of a heatmap or a violin plot. 
#'If methods agree, cell type proportions of the same cell type should 
#'by strongly correlated. For cell types with weak correlation across 
#'methods, corresploding estimated cell type proportions should be 
#'interpreted with caution.
#'
#'@param correlated output object from \code{correlate}
#'
#'@param method plot type ("heatmap" or "boxplot")
#'
#'@param legend boolean to display color legend
#'
#'@return Returns a heatmap or violin plot showing the correlation 
#'distribution of by different methods/signature matrices for each 
#'cell type
#'
#'@import ggplot2
#'@importFrom stats na.omit
#'@importFrom tibble as_tibble
#'@importFrom dplyr mutate filter
#'@importFrom tidyr gather
#'@importFrom pheatmap pheatmap
#'@importFrom ggplotify as.ggplot
#'@importFrom grDevices colorRampPalette
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
#'# correlate
#'correl <- correlate(deconvoluted = decon)
#'
#'# plot correlate
#'plot_correlate(correlated = correl, 
#'method="heatmap")
#'
#'@export
plot_correlate <- function(correlated, method="heatmap", legend=TRUE)
{
    cor_big <- as.matrix(correlated$correlation)

    # remove NAs
    cor_big <- cor_big[
        !apply(cor_big, 1, function(x) sum(na.omit(x))==1),
        !apply(cor_big, 2, function(x) sum(na.omit(x))==1)]

    # plot
    if(method=="boxplot")
    {
        df <- cor_big %>% 
            as_tibble() %>%
            mutate(CellTypeA = rownames(cor_big)) %>%
            gather("CellTypeB","Correlation", -.data$CellTypeA) %>%
            mutate(CellTypeA = gsub("_.*","",.data$CellTypeA)) %>%
            mutate(CellTypeB = gsub("_.*","",.data$CellTypeB)) %>%
            filter(.data$CellTypeA==.data$CellTypeB) %>%
            filter(.data$Correlation!=1)

        if(legend)
            g <- ggplot(df, aes(.data$CellTypeA, .data$Correlation, fill = .data$CellTypeA)) + 
                geom_violin(scale = "width", width=0.8, trim=FALSE) + 
                geom_boxplot(fill="white",alpha=0.3) +
                ylim(c(0,1)) +
                xlab("") +
                ylab("cell type proportions correlation") + 
                guides(fill=guide_legend(title="CellType")) +
                theme_bw() + 
                theme(
                    panel.grid.major.x = element_blank(), 
                    panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    axis.text.x = element_blank(),
                    plot.title = element_text(size = 12, face = 'bold', hjust = 0.5),
                    axis.ticks.x = element_blank())
        else
            g <- ggplot(df, aes(.data$CellTypeA, .data$Correlation, fill = .data$CellTypeA)) + 
                geom_violin(scale = "width", width=0.8, trim=FALSE) + 
                geom_boxplot(fill="white",alpha=0.5) +
                ylim(c(0,1)) + 
                xlab("") + 
                ylab("cell type proportions correlation") +  
                theme_bw() + 
                theme(
                    panel.grid.major.x = element_blank(), 
                    panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    axis.text.x = element_text(angle=90),
                    plot.title = element_text(size = 12, face = 'bold', hjust = 0.5),
                    legend.position = "none")
    } 
    else if(method=="heatmap")
    {
        df <- data.frame(CellType = gsub("_.*","",colnames(cor_big)), row.names = colnames(cor_big))

        breaksList = seq(-1, 1, by = 0.1)

        if(legend)
            g <- as.ggplot(pheatmap(mat = cor_big,
                    color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)),
                    breaks = breaksList,
                    annotation_col = df,
                    annotation_legend = TRUE,
                    annotation_names_row = FALSE,
                    annotation_names_col= FALSE,
                    annotation_colors = list(CellType = gg_color_hue(df$CellType,l=50)),
                    show_rownames = FALSE,
                    show_colnames = FALSE,
                    fontsize_col = 5,
                    border_color = NA,
                    clustering_distance_rows = 'euclidean',
                    treeheight_row = 0,
                    treeheight_col = 3,
                    silent=TRUE)) +
                    theme(plot.title = element_text(size = 12, face = 'bold', hjust = 0.5, margin=margin(0,0,5,0))) +
                    ggtitle('Cell type proprotions correlation')
        else
            g <- as.ggplot(pheatmap(mat = cor_big,
                    color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)),
                    breaks = breaksList,
                    annotation_col = df,
                    annotation_legend = FALSE,
                    annotation_names_row = FALSE,
                    annotation_names_col= FALSE,
                    show_rownames = FALSE,
                    show_colnames = FALSE,
                    fontsize_col = 5,
                    border_color = NA,
                    clustering_distance_rows = 'euclidean',
                    treeheight_row = 0,
                    treeheight_col = 3,
                    silent=TRUE)) +
                    theme(plot.title = element_text(size = 12, face = 'bold', hjust = 0.5, margin=margin(0,0,5,0))) +
                    ggtitle('Cell type proprotions correlation')
    } else {
        stop('Choose either "boxplot" or "heatmap" as plotting option.')
    }
    return(g)
}
