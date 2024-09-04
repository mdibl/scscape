# ![nf-core/scscape](docs/images/nf-core-scscape_logo_light.png#gh-light-mode-only) ![nf-core/scscape](docs/images/nf-core-scscape_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/scscape/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/scscape/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/scscape/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/scscape/actions?query=workflow%3A%22nf-core+linting%22)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/scscape/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/scscape)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23scscape-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/scscape)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/scscape** is a bioinformatics pipeline that was built for multi-sample single cell analysis downstream from the generation of count matrices.
The pipeline operates using many functional components derived from the [Seurat](https://satijalab.org/seurat/) R package. Input data is expected to be in the
format of barcodes, features, and matrix files. Output includes Seurat objects that contain QC metrics, identified cell clusters, and dimensionally reduced projections that
encompass the experiments gene expression variability.

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. Gzip all raw input files for consistency
2. Initialize seurat object for each sample
3. Normalize gene expression counts & perform mitochondrial / cell-cycle scoring
4. Detect and remove suspected doublets from each sample
5. Merge - normalize - find variable features - scale data (SCTransform)
6. Run principal component analysis
7. Perform integration to remove technical confounding variables
8. Find k nearest-neighbors & cluster (Louvain)
9. Dimensionally reduce expression variance and plot

## Documentation

The nf-core/scscape pipeline comes with documentation about the pipeline

![scscape workflow](docs/images/SubwayMap.png)

## Usage


>**Note:**
>If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.


## Configuration

First, prepare a sample sheet with your input data that looks as follows:

`Samples.csv`:

```csv
id,data_directory,mt_cc_rm_genes
00_dpa_1,/filtered_feature_bc_matrix/,AuxillaryGeneList.csv
```
Each row represents a samples matrix files (barcodes.tsv, features.tsv, matrix.mtx) and associated genes used in the analysis.

Second, add mitochondrial, S phase, G2 / M phase, removal genes

`AuxillaryGeneList.csv`:

```csv
MTgenes,G2Mgenes,Sgenes,RMgenes
mt-nd1,hmgb2a,mcm5,
mt-nd2,cdk1,pcna,
```

Finally, construct a segmentation file defining the analysis groups for the experiment (ex: treatment, rep, age, sex).

`segmentation.csv`:

```csv
id,00_dpa,04_dpa,all
00_dpa_1,true,false,true
00_dpa_2,true,false,true
04_dpa_1,false,true,true
04_dpa_2,false,true,true
```

>***Make sure id columns match between `segmentation.csv` & `Samples.csv`***

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run nf-core/scscape \
   -profile docker \
   -params-file paramaters.json \
   -c custom.config
```

>**Warning:**
>Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

>**Note:**
>There is the ability to create a `.loupe` file within the configuration options of this pipeline. This file can be used with the [10X Loupe Browser](https://www.10xgenomics.com/support/software/loupe-browser/latest) to interactively explore you single cell experiment. In order to successfully generate the file, you are required by 10X to both read the [10X End User License Agreement](https://www.10xgenomics.com/legal/end-user-software-license-agreement) and accept their terms by setting the `eula_agreement` parameter to `Agree` (in addition to setting `makeLoupe` to `true`).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/scscape/usage) and the [parameter documentation](https://nf-co.re/scscape/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/scscape/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/scscape/output).

## Credits

nf-core/scscape was originally written by Ryan Seaman, Riley Grindle, Joel Graber.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#scscape` channel](https://nfcore.slack.com/channels/scscape) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/scscape for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
