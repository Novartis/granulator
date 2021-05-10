granulator: Rapid benchmarking of methods for *in silico* deconvolution of 
bulk RNA-seq data
================

<!-- badges: start -->
[![Lifecycle:stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
<!-- badges: end -->

## Introduction

Heterogeneity in the cellular composition of bulk RNA-seq data may
prevent or bias the results from differential expression analysis. To
circumvent this limitation, *in silico* deconvolution infers cell type
abundances by modelling gene expression levels as weighted sums of the
cell-type specific expression profiles. Several computational methods
have been developed to estimate cell type proportions from bulk
transcriptomics data, and to account for cell type heterogeneity in the
statistical analysis. The R package *granulator* provides a unified
testing interface to rapidly run and benchmark multiple state-of-the-art
deconvolution methods. We demonstrate its usage on published bulk
RNA-seq data from peripheral blood mononuclear cells.

## Methods

The methods currently implemented in *granulator* are reported in
**Table 1**.

<table>
<colgroup>
<col style="width: 11%" />
<col style="width: 8%" />
<col style="width: 27%" />
<col style="width: 35%" />
<col style="width: 16%" />
</colgroup>
<thead>
<tr class="header">
<th>Name</th>
<th>Function</th>
<th>Method</th>
<th>License</th>
<th>Reference</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>ols</td>
<td>stats::lsfit</td>
<td>Ordinary least squares</td>
<td>free (<a href="https://cran.r-project.org/web/packages/L1pack/L1pack.pdf">GPL-2</a>)</td>
<td></td>
</tr>
<tr class="even">
<td>nnls</td>
<td>nnls::nnls</td>
<td>Non-negative least squares</td>
<td>free (<a href="https://cran.r-project.org/web/packages/nnls/index.html">GPL-2, GPL-3</a>)</td>
<td>reimplemented based on <span class="citation" data-cites="Abbas2009">(Abbas et al. 2009)</span></td>
</tr>
<tr class="odd">
<td>qprogwc</td>
<td>limSolve::lsei</td>
<td>Quadratic programming with non-negativity and sum-to-one constraint</td>
<td>free (<a href="https://cran.r-project.org/web/packages/limSolve/index.html">GPL-2, GPL-3</a>)</td>
<td>reimplemented based on <span class="citation" data-cites="Gong2013">(Gong and Szustakowski 2013)</span></td>
</tr>
<tr class="even">
<td>qprog</td>
<td>limSolve::Solve</td>
<td>Quadratic programming without constraints</td>
<td>free (<a href="https://cran.r-project.org/web/packages/limSolve/index.html">GPL-2, GPL-3</a>)</td>
<td></td>
</tr>
<tr class="odd">
<td>rls</td>
<td>MASS::rlm</td>
<td>Re-weighted least squares</td>
<td>free (<a href="https://cran.r-project.org/web/packages/MASS/index.html">GPL-2, GPL-3</a>)</td>
<td>reimplemented based on <span class="citation" data-cites="Monaco2019">(Monaco et al. 2019)</span></td>
</tr>
<tr class="even">
<td>svr</td>
<td>e1071::svr</td>
<td>Support vector regression</td>
<td>free (<a href="https://cran.r-project.org/web/packages/e1071/index.html">GPL-2, GPL-3</a>)</td>
<td>reimplemented based on <span class="citation" data-cites="Newman2015">(Newman et al. 2015)</span></td>
</tr>
<tr class="odd">
<td>dtangle</td>
<td>dtangle::dtangle</td>
<td>Linear mixing model</td>
<td>free (<a href="https://cran.r-project.org/web/packages/dtangle/index.html">GPL-3</a>)</td>
<td><span class="citation" data-cites="Hunt2018">(Hunt et al. 2018)</span></td>
</tr>
</tbody>
</table>

**Table 1** - Deconvolution methods. List of deconvolution algorithms
available in *granulator*.

## Installation

*granulator* can be installed from Bioconductor using:

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("granulator")

The package can be loaded using:

    library(granulator)

## Data

The datasets included in the package comprises bulk RNA-seq gene
expression data of peripheral blood mononuclear cells (PBMCs) from 12
healthy donors and bulk RNA-seq data of 29 isolated immune cell types
from 4 healthy donors (Monaco et al. 2019), publicly available at NCBI
database under GEO accession number
[GSE107011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011).

## Vignettes

We show how to use *granulator* for the deconvolution of bulk RNA-seq
data from peripheral blood mononuclear cells (PBMCs) into the individual
cellular components and how to assess the quality of the obtained
predictions in the following vignette: [Deconvolution of bulk RNA-seq
data with
granulator](https://bioconductor.org/packages/granulator).

## References

Abbas, Alexander R., Kristen Wolslegel, Dhaya Seshasayee, Zora Modrusan,
and Hilary F. Clark. 2009. “<span class="nocase">Deconvolution of blood
microarray data identifies cellular activation patterns in systemic
lupus erythematosus</span>.” *PLoS ONE* 4 (7).
<https://doi.org/10.1371/journal.pone.0006098>.

Gong, Ting, and Joseph D Szustakowski. 2013. “<span
class="nocase">DeconRNASeq: a statistical framework for deconvolution of
heterogeneous tissue samples based on mRNA-Seq data</span>.”
*Bioinformatics* 29 (8): 1083–85.
<https://doi.org/10.1093/bioinformatics/btt090>.

Hunt, Gregory J, Saskia Freytag, Melanie Bahlo, and Johann A
Gagnon-Bartsch. 2018. “Dtangle: Accurate and Robust Cell Type
Deconvolution.” *Bioinformatics* 35 (12): 2093–99.
<https://doi.org/10.1093/bioinformatics/bty926>.

Monaco, Gianni, Bernett Lee, Weili Xu, Seri Mustafah, You Yi Hwang,
Christophe Carré, Nicolas Burdin, et al. 2019. “<span
class="nocase">RNA-Seq Signatures Normalized by mRNA Abundance Allow
Absolute Deconvolution of Human Immune Cell Types</span>.” *Cell
Reports* 26 (6): 1627–1640.e7.
<https://doi.org/10.1016/j.celrep.2019.01.041>.

Newman, Aaron M., Chih Long Liu, Michael R. Green, Andrew J. Gentles,
Weiguo Feng, Yue Xu, Chuong D. Hoang, Maximilian Diehn, and Ash A.
Alizadeh. 2015. “<span class="nocase">Robust enumeration of cell subsets
from tissue expression profiles</span>.” *Nature Methods* 12 (5):
453–57. <https://doi.org/10.1038/nmeth.3337>.
