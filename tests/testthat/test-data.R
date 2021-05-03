context("Load data")

test_that("Load data", {

    # load data
    load_ABIS()

    # test sigMatrix_ABIS_S0
    expect_equal(class(sigMatrix_ABIS_S0)[1],"matrix")
    expect_equal(class(plot_similarity(sigMatrix_ABIS_S0))[1],"gg")

    # test sigMatrix_ABIS_S1
    expect_equal(class(sigMatrix_ABIS_S1)[1],"matrix")
    expect_equal(class(plot_similarity(sigMatrix_ABIS_S1))[1],"gg")

    # test sigMatrix_ABIS_S2
    expect_equal(class(sigMatrix_ABIS_S2)[1],"matrix")
    expect_equal(class(plot_similarity(sigMatrix_ABIS_S2))[1],"gg")

    # test sigMatrix_ABIS_S3
    expect_equal(class(sigMatrix_ABIS_S3)[1],"matrix")
    expect_equal(class(plot_similarity(sigMatrix_ABIS_S3))[1],"gg")

    # test bulkRNAseq_ABIS
    expect_equal(class(bulkRNAseq_ABIS)[1],"matrix")
    
    # test groundTruth_ABIS
    expect_equal(class(groundTruth_ABIS)[1],"matrix")
})
