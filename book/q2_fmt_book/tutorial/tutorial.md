# q2-fmt Tutorial
(exploring-the-data)=
## Exploring the data

```{usage-scope}
---
name: tutorial
---
```

```{usage-selector}
---
default-interface: cli-usage
---
```

(access-metadata)=
### Access and summarize the study metadata

To begin our work with QIIME 2 and the tutorial data we will
start by downloading the metadata, generating a summary, and exploring
that summary.

First, download the metadata.

```{usage}
md_url = 'https://qiime2-workshops.s3.us-west-2.amazonaws.com/itn-aug2024/sample-metadata-v3-q2-fmt.tsv'

sample_metadata = use.init_metadata_from_url('sample_metadata', md_url)
```

Next, we’ll get a view of the study metadata using QIIME 2. This will allow you to assess whether the metadata that QIIME 2 is using is as you expect. You can do this using the tabulate action in QIIME 2’s q2-metadata plugin as follows.

```{usage}
use.action(
    use.UsageAction(plugin_id='metadata', action_id='tabulate'),
    use.UsageInputs(input=sample_metadata),
    use.UsageOutputNames(visualization='metadata_summ_1')
)
```

(access-feature-table)=
### Access and summarize the feature table

The feature table will describe the amplicon sequence variants (ASVs) observed in each sample, and how many times each ASV was observed in each sample. The feature data in this case is the sequence that defines each ASV.

In this tutorial, we're going to work specifically with samples that were
included in the autoFMT randomized trial.

Lets generate and explore a summary of the feature table we will be using.

```{usage}

feature_table_url = 'https://qiime2-workshops.s3.us-west-2.amazonaws.com/itn-aug2024/autofmt-table.qza'

autofmt_table = use.init_artifact_from_url('feature-table', feature_table_url)
```

```{usage}
use.action(
    use.UsageAction(plugin_id='feature_table', action_id='summarize'),
    use.UsageInputs(table=autofmt_table, sample_metadata=sample_metadata),
    use.UsageOutputNames(visualization='autofmt_table_summ'),
)
```
## Selecting an Even Sampling Depth

### Rarefying, rarefaction, and q2-boots

A first step in analyzing our microbiome feature table is to choose an "even sampling depth", or the number of sequences that we should select at random from each of our samples to ensure that all samples are sequenced at equivalent depth or with equivalent effort.
This processes is referred to as rarefying our feature table.
Rarefying feature tables, or sampling them to a user-specified sampling depth and discarding samples with a total frequency that is less than the sampling depth, is a bit of a controversial topic because it throws away some of the data that was collected and is subject to biases as a result.
*Rarefaction* differs in a subtle way from rarefying, in that it involves repeat sampling from the input feature table to a user-specified sampling depth.
These concepts were recently discussed in {cite}`Schloss2024-aq`.
In this tutorial, for the sake of time, we are going to focus our diversity analyses on rarefying (i.e., a single iteration of random sampling).
To perform rarefaction-based diversity analysis with QIIME 2, refer to the [q2-boots](https://q2-boots.readthedocs.io/en/latest/) plugin {cite}`Raspet2024-om`.

### Selecting an even sampling depth

To start our diversity analyses, we first need to determine what even sampling depth (or "rarefaction depth") we want to select for computing our diversity metrics.
Because most diversity metrics are sensitive to different sampling depths across different samples, it is common to randomly subsample the counts from each sample to a specific value.
For example, if you define your sampling depth as 500 sequences per sample, the counts in each sample will be subsampled without replacement so that each sample in the resulting table has a total count of 500.
If the total count for any sample(s) are smaller than this value, those samples will be dropped from the downstream analyses. Choosing this value is tricky.
We recommend making your choice by reviewing the information presented in the feature table summary file.
Choose a value that is as high as possible (so you retain more sequences per sample) while excluding as few samples as possible.

Open up the feature table summary that you previously created with either Galaxy or in your QIIME 2 container and we'll discuss this as a group.

### Alpha rarefaction plots

After choosing an even sampling depth, it's helpful to see if your diversity metrics appear stable at that depth of coverage.
You can do this for alpha diversity using an alpha rarefaction plot.

```{usage}
use.action(
    use.UsageAction(plugin_id='diversity', action_id='alpha_rarefaction'),
    use.UsageInputs(table=autofmt_table, metrics={'observed_features'},
                    metadata=sample_metadata, max_depth=33000),
    use.UsageOutputNames(visualization='obs-features-alpha-rarefaction'))
```

## Computing diversity metrics



The next step that we'll work through is computing a series of common diversity metrics on our feature table.
We'll do this using the `q2-diversity` plugin's `core-metrics` action.
This action is a QIIME 2 `Pipeline` which combines over ten different actions in a single command.

### Core diversity metrics

The `core-metrics` action requires your feature table and your sample metadata as input.
It additionally requires that you provide the sampling depth that this analysis will be performed at.
Determining what value to provide for this parameter is often one of the most confusing steps of an analysis.
Refer to the previous chapter for details.
Here we prioritize retaining samples and so we select a sampling depth of 10,000.

```{usage}
core_metrics_results = use.action(
    use.UsageAction(plugin_id='diversity', action_id='core_metrics'),
    use.UsageInputs(table=autofmt_table,
                    sampling_depth=10000, metadata=sample_metadata),
    use.UsageOutputNames(rarefied_table='rarefied_table',
                            observed_features_vector='observed_features_vector',
                            shannon_vector='shannon_vector',
                            evenness_vector='evenness_vector',
                            jaccard_distance_matrix='jaccard_distance_matrix',
                            bray_curtis_distance_matrix='bray_curtis_distance_matrix',
                            jaccard_pcoa_results='jaccard_pcoa_results',
                            bray_curtis_pcoa_results='bray_curtis_pcoa_results',
                            jaccard_emperor='jaccard_emperor',
                            bray_curtis_emperor='bray_curtis_emperor'),
)
```

As you can see, this command generates many outputs including both QIIME 2 artifacts and visualizations.
We'll work together on a guided exploration of these results.
