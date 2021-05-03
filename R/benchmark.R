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


#'@title Regress estimated cell type coefficients against the ground truth
#'
#'@description \code{regress} computes regression between estimated cell 
#'type proportions and the measured cell type  proportions (ground truth).
#'
#'@param deconvoluted Output object of the function \code{deconvolute}.
#'
#'@param ground_truth A matrix containing measured cell type 
#'proportions. Samples names are inlcuded in rownames. In the case multiple 
#'signatures were tested, the corresponding data frames containing the 
#'measured proportions can be given as a list.
#'
#'@return Returns a list containing thres elements: \itemize{
#'\item{data: a list of data frames with celltype matched estimated and 
#'predicted proportions}
#'\item{stats: a list of data frames with regression statistics comprising 
#'Pearson Correlation Coefficient ('pcc'), Concordance Correlation Coefficient
#'('ccc'), Coefficient of Determination ('adj.r2') and Root Mean Square Error 
#'('rmse')}
#'\item{summary: a data frame with summary statistics by cell type}
#'\item{rank: ranking of deconvolution alghoritms by highest all-to-all 
#'correlation of coefficients}
#'\item{summay: summary statistics of regression coefficients by method, 
#'signature and cell type}
#'\item{rank: ranking of methods and signatures by highest average regression 
#'coefficient}
#'\item{combinations: combination of methods and signatures tested}
#'}
#'
#'@importFrom magrittr %>%
#'@importFrom dplyr select mutate bind_rows left_join group_by summarize ungroup arrange
#'@importFrom tidyr gather nest
#'@importFrom tibble as_tibble rownames_to_column
#'@importFrom purrr map
#'@importFrom purrr map_dbl
#'@importFrom epiR epi.ccc
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
#'@export
benchmark <- function(deconvoluted, ground_truth){

    # input check
    if (!is.matrix(ground_truth))
        stop('Bulk RNA-seq should be a matrix; got ', class(ground_truth))

    # fix cell type names
    ground_truth <- fix_col_names(ground_truth)

    # results 
    res <- list(
            data = lapply(deconvoluted$coefficients, function(x) NULL),
            stats = lapply(deconvoluted$coefficients, function(x) NULL),
            summary = NULL)

    # stats
    for (i in seq_along(deconvoluted$coefficients))
    {  
        x <- names(deconvoluted$coefficients)[i]
        m <- deconvoluted$combinations$method[i]
        s <- deconvoluted$combinations$signature[i]

        # coefficients
        pred <- deconvoluted$coefficients[[x]]
        prop <- ground_truth

        # remove nas
        pred[is.na(pred)] <- 0
        prop[is.na(prop)] <- 0

        # match columns
        pred <- pred[,intersect(colnames(pred),colnames(prop))]
        prop <- prop[,colnames(pred)]

        # match samples
        pred <- pred[intersect(rownames(pred),rownames(prop)),]
        prop <- prop[rownames(pred),]

        # check input
        if (ncol(pred)==0)
            stop('No matching cell types were found!')

        if (nrow(pred)==0)
            stop('No matching samples were found!')

        # data
        prop <- prop %>%
            as.data.frame() %>%
            rownames_to_column(var = 'sample') %>%
            gather('celltype', 'proportion', -.data$sample)

        pred <- pred %>%
            as.data.frame() %>%
            rownames_to_column(var = 'sample') %>%
            gather('celltype', 'prediction', -.data$sample)

        res$data[[x]] <- left_join(prop, pred, by = c("sample", "celltype"))

        # stats
        stats <- res$data[[x]] %>% 
            as_tibble() %>%
            select(-.data$sample) %>%
            mutate(proportion = .data$proportion/100) %>%
            mutate(prediction = .data$prediction/100) %>%
            nest(thisdata = c(.data$prediction, .data$proportion)) %>%
            mutate(
                method = m,
                signature = s,
                pcc = map_dbl(.data$thisdata, ~ suppressMessages(suppressWarnings(cor(.x$prediction, .x$proportion)))),
                ccc = map_dbl(.data$thisdata, ~ suppressMessages(suppressWarnings(epi.ccc(.x$prediction, .x$proportion, ci = "z-transform",conf.level = 0.95,rep.measure = FALSE)[['rho.c']][['est']]))),
                lmmodel = map(.data$thisdata, ~lm(prediction ~ proportion, data = .x)),
                adj.r2 = map_dbl(.data$lmmodel, ~signif(summary(.x)$adj.r.squared,2)),
                rmse = map_dbl(.data$lmmodel, ~signif(sqrt(mean(.x$residuals^2)),2))) %>%
            select(-.data$thisdata,-.data$lmmodel) %>%
            as.data.frame()

        res$stats[[x]] <- stats
    }

    # summary
    res$summary <- bind_rows(res$stats) %>%
        select(.data$signature,.data$method,.data$celltype,.data$pcc,.data$ccc,.data$adj.r2,.data$rmse) %>%
        mutate(pcc = round(.data$pcc,4)) %>%
        mutate(ccc = round(.data$ccc,4)) %>%
        mutate(adj.r2 = round(.data$adj.r2,4)) %>%
        mutate(rmse = round(.data$rmse,4)) %>%
        as.data.frame()

    # rank
    res$rank <- res$summary %>%
        group_by(.data$method, .data$signature) %>%
        summarize(
                mean_pcc=round(mean(.data$pcc),4),
                mean_ccc=round(mean(.data$ccc),4),
                mean_adj.r2=round(mean(.data$adj.r2),4),
                mean_rmse=round(mean(.data$rmse),4)) %>%
        ungroup() %>%
        select(.data$signature,.data$method,.data$mean_pcc,.data$mean_ccc,.data$mean_pcc,.data$mean_adj.r2,.data$mean_rmse) %>%
        arrange(desc(.data$mean_pcc)) %>%
        as.data.frame()

    # combinations
    res$combinations <- deconvoluted$combinations

    return(res)
}
