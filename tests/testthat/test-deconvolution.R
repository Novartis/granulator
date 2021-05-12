context("Deconvolution methods")

test_that("Deconvolution methods", {

    # load data
    load_ABIS()

    for(m in get_decon_methods()){

        # tpms
        mat <- round(matrix(rexp(200, rate=.01), ncol=20))
        len <- round(matrix(rexp(10, rate=.001), ncol=1))+10
        tpm <- get_TPM(mat,as.vector(len))
        expect_equal(class(tpm)[1],"matrix")

        # deconvolution
        decon <- deconvolute(bulkRNAseq_ABIS,sigMatrix_ABIS_S0,m)
        expect_equal(class(decon[[1]][[1]]),"data.frame")
        expect_equal(class(decon[[2]][[1]]),"data.frame")
        expect_equal(class(decon[[3]]),"data.frame")
        expect_equal(class(plot_proportions(decon, method = m, signature = 'sig1'))[1],"gg")
        expect_equal(class(plot_deconvolute(decon))[1],"gg")

        # benchmark
        bench <- benchmark(decon, groundTruth_ABIS)
        expect_equal(class(bench[[1]][[1]]),"data.frame")
        expect_equal(class(bench[[2]][[1]]),"data.frame")
        expect_equal(class(bench[[3]]),"data.frame")
        expect_equal(class(bench[[4]]),"data.frame")
        expect_equal(class(bench[[5]]),"data.frame")
        expect_equal(class(plot_regress(bench, method = m, signature = 'sig1'))[1],"gg")
        expect_equal(class(plot_benchmark(bench))[1],"gg")

        # correlation
        correl <- correlate(decon)
        expect_equal(class(correl[[1]]),"data.frame")
        expect_equal(class(correl[[2]]),"data.frame")
        expect_equal(class(correl[[3]]),"data.frame")
        expect_equal(class(correl[[4]]),"data.frame")
        expect_equal(class(plot_correlate(correl))[1],"gg")

    }
})

