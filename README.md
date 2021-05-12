# ![nf-core/pikavirus](docs/images/nf-core-pikavirus_logo.png)

**A workflow for metagenomics**.

[![GitHub Actions CI Status](https://github.com/nf-core/pikavirus/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/pikavirus/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/pikavirus/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/pikavirus/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/pikavirus.svg)](https://hub.docker.com/r/nfcore/pikavirus)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23pikavirus-4A154B?logo=slack)](https://nfcore.slack.com/channels/pikavirus)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**nf-core/pikavirus** is a bioinformatics best-practise analysis pipeline for metagenomic analysis following a new approach, based on eliminatory k-mer analysis, followed by assembly and posterior contig-binning.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Podman`](https://podman.io/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/pikavirus -profile test,<docker/singularity/podman/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    ```bash
    nextflow run nf-core/pikavirus -profile <docker/singularity/podman/conda> --input '*_R{1,2}.fastq.gz'
    ```

See [usage docs](https://nf-co.re/pikavirus/usage) for all of the available options when running the pipeline.

## Pipeline Summary

By default, the pipeline currently performs the following:


* Sequencing quality control (`FastQC`)
* Trimming of low-quality regions in the reads (`FastP`)
* Identification and elimination of reads from the host, and isolation of viral, bacterial, fungal and unknown reads (`Kraken2`)
* Assembly of unknow reads (`MetaQuast`) and mapping against databases (`Kaiju`) to identify new possible pathogens
* Alignment of viral, bacterial and fungal reads against reference genomes to ensure the presence of certain organisms (`Bowtie2`)

## Documentation

The nf-core/pikavirus pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/pikavirus/usage) and [output](https://nf-co.re/pikavirus/output).

<!-- TODO nf-core: Add a brief overview of what the pipeline does and how it works -->

## Credits

nf-core/pikavirus was originally written by Guillermo Jorge Gorines Cordero, under supervision of the BU-ISCIII team in Madrid, Spain.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#pikavirus` channel](https://nfcore.slack.com/channels/pikavirus) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/pikavirus for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
> ReadCube: [Full Access Link](https://rdcu.be/b1GjZ)

In addition, references of tools and data used in this pipeline are as follows:

> **Improved metagenomic analysis with Kraken 2.**
>
> Derrick E Wood, Jennifer Lu & Ben Langmead.
>
> _Genome biology_ 2019 Nov 28. doi: [10.1186/s13059-019-1891-0](https://doi.org/10.1186/s13059-019-1891-0)
>

>**fastp: an ultra-fast all-in-one FASTQ preprocessor.**
>
> Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu.
> 
> _Bioinformatics_, Volume 34, Issue 17, 01 September 2018, Pages i884–i890,. doi: [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)
>

>**Fast and sensitive taxonomic classification for metagenomics with Kaiju**
>
> Peter Menzel, Kim Lee Ng & Anders Krogh
>
> _Nature Communications_ volume 7, Article number: 11257 (2016). doi [10.1038/ncomms11257](https://doi.org/10.1038/ncomms11257)
>

>**QUAST: quality assessment tool for genome assemblies**
>
> Alexey Gurevich, Vladislav Saveliev, Nikolay Vyahhi & Glenn Tesler 
>
> _Bioinformatics_ Volume 29, Issue 8, 15 April 2013, Pages 1072–1075. doi [10.1093/bioinformatics/btt086](https://doi.org/10.1093/bioinformatics/btt086)
>

>**Bioconda: sustainable and comprehensive software distribution for the life sciences**
>
> Björn Grüning, Ryan Dale, Andreas Sjödin, Brad A. Chapman, Jillian Rowe, Christopher H. Tomkins-Tinch, Renan Valieris, Johannes Köster & The Bioconda Team
>
> _Nature Methods_ volume 15, pages 475–476(2018). doi [10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7)
>

>**Mash: fast genome and metagenome distance estimation using MinHash**
>
> Brian D. Ondov, Todd J. Treangen, Páll Melsted, Adam B. Mallonee, Nicholas H. Bergman, Sergey Koren & Adam M. Phillippy
>
> _Genome Biology_ 17, Article number: 132 (2016). doi [10.1186/s13059-016-0997-x](https://www.doi.org/10.1186/s13059-016-0997-x)
>