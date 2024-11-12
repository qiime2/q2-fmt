# QIIME 2 FMT Book

This tutorial was orignally made for the 2024 [ITN](https://www.itcrtraining.org/) / [QIIME 2](https://qiime2.org) Workshop, taught 27 - 29 August 2024 at the National Institutes of Health. 
If you are interested in the original workshop tutorial, please visit: https://q2-cat-book.readthedocs.io/en/2024.08.27-q2-itn/. 
This tutorial is now the offical documentation for the q2-fmt plugin.

```{admonition} Code of Conduct
Please review the [QIIME 2 Community Code of Conduct](https://forum.qiime2.org/t/qiime-2-community-code-of-conduct/9057).
All workshop instructors agreed to abide by these community expectations and we hope you will too as you become involved in the QIIME 2 community.
```

This tutorial will walk you through all the functionality available in q2-fmt. 
We will start with a metadata file and feature table. 
We will discuss choosing an even sampling depth and generating diversity metrics so that we are prepared to run our data through q2-fmt.
Grab a cup of coffee or tea ☕, and lets jump right in!

## Recommended readings

The following readings will provide background information for this tutorial.
These are presented in priority order.

1. [_Assessing Engraftment Following Fecal Microbiota Transplant_](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11042410/) ({cite:t}`herman_2024`)
1. [_Reconstitution of the gut microbiota of antibiotic-treated patients by autologous fecal microbiota transplant_](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6468978/) ({cite:t}`taur-autofmt-2018`)
1. [_Compilation of longitudinal microbiota data and hospitalome from hematopoietic cell transplantation patients_](https://www.nature.com/articles/s41597-021-00860-8) ({cite:t}`liao-data-2021`)
1. [_Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2_](https://doi.org/10.1038/s41587-019-0209-9) ({cite:t}`bolyen-qiime2-2019`)
1. [_q2-longitudinal: Longitudinal and Paired-Sample Analyses of Microbiome Data_](http://dx.doi.org/10.1128/mSystems.00219-18) ({cite:t}`bokulich-q2long-2018`)
1. [Getting started with QIIME 2](https://gregcaporaso.github.io/q2book/using/getting-started.html)
1. [_Experiences and lessons learned from two virtual, hands-on microbiome bioinformatics workshops_](https://doi.org/10.1371/journal.pcbi.1009056) ({cite:t}`dillon-workshops-2021`)

## Thank you to our funders

This tutorial  was funded in part by NIH National Cancer Institute Informatics Technology for Cancer Research grant [1U24CA248454-01](https://reporter.nih.gov/project-details/9951750) to JGC.

The preparation and dissemination of the data used in this tutorial was funded in part by the National Institutes of Health (NIH) grants U01 AI124275, R01 AI137269 and U54 CA209975 to Joao B. Xavier.

This website is built with MyST Markdown and Jupyter Book, which are supported in part with [funding](https://sloan.org/grant-detail/6620) from the Alfred P. Sloan Foundation.

## Reusing

 <p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><a property="dct:title" rel="cc:attributionURL" href="https://q2-fmt.readthedocs.io/en/latest/intro.html">q2-fmt-book</a> by <a rel="cc:attributionURL dct:creator" property="cc:attributionName" href="https://cap-lab.bio">Caporaso Lab</a> is licensed under <a href="https://creativecommons.org/licenses/by-nc-sa/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY-NC-SA 4.0.<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/nc.svg?ref=chooser-v1" alt=""><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/sa.svg?ref=chooser-v1" alt=""></a></p>

 You are free distribute, remix, adapt, and build upon the material in any medium or format, for noncommercial purposes only, as long as credit is given to the original creator. Adaptations must be shared under the same terms.


## Table of Contents

```{tableofcontents}
```
