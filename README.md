# q2-fmt (fecal microbiota transplant)

![](https://github.com/qiime2/q2-fmt/workflows/ci/badge.svg)

## Demo
[![](https://raw.githubusercontent.com/qiime2/q2-fmt/master/demo/screenshot.png)
**Interactive Link**](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fraw.githubusercontent.com%2Fqiime2%2Fq2-fmt%2Fmaster%2Fdemo%2Fraincloud-baseline0.qzv)

## Installation

Please note that this software is now an alpha release. Conda packages will
be provided in the future. In the interim, use this command in the latest
QIIME 2 environment.

```bash
pip install git+https://github.com/qiime2/q2-fmt.git
```


## Examples
Example data can be generated using the `--example-data` flag on each action
described below. This will create a directory structure to match the examples
(you will need to `cd` into the appropriate directory first).

```bash
qiime fmt --example-data fmt-examples/
cd fmt-examples/engraftment/engraftment-baseline/
```

**Plot Engraftment**

This pipeline will execute the below actions in order to produce engraftment
plots in a single step. To learn more, see the examples below this.
```bash
qiime fmt engraftment \
  --i-diversity-measure div-measure.qza \
  --m-metadata-file md.tsv \
  --p-hypothesis baseline \
  --p-time-column week \
  --p-reference-column InitialDonorSampleID \
  --p-subject-column SubjectID \
  --p-where '[SampleType]="stool"' \
  --p-filter-missing-references \
  --p-against-group 0 \
  --p-p-val-approx asymptotic \
  --o-stats stats.qza \
  --o-raincloud-plot raincloud-plot.qzv
```

---

**Group Timepoints**

First the data must be collected from your study and organized into groups with
some relevant diversity measure. We will use Faith's Phylogenetic Diversity in
these examples, but other measures, including beta diversity metrics, are
supported.

This command will produce two collections of groups: one for subjects at each
timepoint, and another for references/donors and any control groups.
```bash
qiime fmt group-timepoints \
  --i-diversity-measure div-measure.qza \
  --m-metadata-file md.tsv \
  --p-where '[SampleType]="stool"' \
  --p-time-column week \
  --p-reference-column InitialDonorSampleID \
  --p-subject-column SubjectID \
  --p-filter-missing-references \
  --o-timepoint-dists timepoint-dists.qza \
  --o-reference-dists reference-dists.qza
```

**Statistical Tests**

Wilcoxon Signed Rank tests and Mann-Whitney U tests are available and
differentiated by the semantic type of the input to ensure the applicability of
the test (i.e. only matched pairs for Wilcoxon, and independent groups for
Mann-Whitney).

Hypotheses for Wilcoxon Signed Rank:
 - baseline: Compare each timepoint against a reference timpoint
 - consecutive: Compare each timepoint against its next timepoint

Hypotheses for Mann-Whitney U:
 - reference: Compare other groups against a particular reference group
 - all-pairwise: Compare all pairwise groups.

The `against-each` parameter is needed to compare the timepoints against the
references/controls. The `engraftment` pipeline above does this
automatically.

```bash
qiime fmt wilcoxon-srt \
  --i-distribution timepoint-dists.qza \
  --p-hypothesis baseline \
  --p-baseline-group 0 \
  --p-p-val-approx asymptotic \
  --o-stats stats_baseline0.qza
```

**Raincloud Plots**

A statistical test is not required to generate the plots, but if provided, will
produce a table showing the results of the test.
```bash
qiime fmt plot-rainclouds \
  --i-data timepoint-dists.qza \
  --i-stats stats_baseline0.qza \
  --o-visualization raincloud-baseline0.qzv
```

