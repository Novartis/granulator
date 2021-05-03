context("Deconvolution methods")

test_that("Deconvolution methods", {

    # load data
    load_ABIS()

    for(m in c('ols','nnls','qprog','qprogwc','rls','svr','dtangle')){

        # deconvolution
        decon <- deconvolute(bulkRNAseq_ABIS,sigMatrix_ABIS_S0,m)
        expect_equal(class(decon[[1]][[1]]),"data.frame")
        expect_equal(class(plot_proportions(decon, method = m, signature = 'sig1'))[1],"gg")
        expect_equal(class(plot_deconvolute(decon))[1],"gg")

        # benchmark
        bench <- benchmark(decon, groundTruth_ABIS)
        expect_equal(class(bench[[1]][[1]]),"data.frame")
        expect_equal(class(plot_regress(bench, method = m, signature = 'sig1'))[1],"gg")
        expect_equal(class(plot_benchmark(bench))[1],"gg")

        # correlation
        correl <- correlate(decon)
        expect_equal(class(correl[[1]]),"data.frame")
        expect_equal(class(plot_correlate(correl))[1],"gg")

    }
})

