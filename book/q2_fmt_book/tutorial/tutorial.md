# Assessing engraftment extent with q2-FMT
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

## Assessing engraftment extent with q2-FMT

```{warning}
This section of this Jupyter Book is still under development 🛠️. Content in this section is subject to change.
```
When investigating FMTs, it is really important to understand the extent to which the recipient microbiome engrafted the donated microbiome. Without assessing engraftment extent, we can never fully understand the clinical outcomes of the study {numref}`engraftment-extent`.

```{figure} _images/engraftment-extent.png
:name: engraftment-extent

Created with BioRender.
```

[Herman et al. (2024)](https://pubmed.ncbi.nlm.nih.gov/38659636/) {cite}`herman_2024` defines three criteria that are important for understanding if a microbiome engrafted following FMT.
These are:
1. Chimeric Asymmetric Community Coalescence
1. Donated Microbiome Indicator Features
1. Temporal Stability.

q2-fmt is a QIIME 2 plugin that was designed to help you investigate all three of these criteria.
For additional infomation, Here is a lecture video of Chloe Herman discussing the importance of assessing engraftment extent![](https://www.youtube.com/watch?v=5vTTDp-t0Bw). 

In this tutorial we will be using q2-fmt to investigating all of these criteron.
## Chimeric Asymmetric Community Coalescence

First, we'll assess the *distance to donor* (using Jaccard distance) to see how the recipient's distance to their donor sample changes following the FMT.
Before continuing, think about what you might expect to see.
Should you expect a recipient's gut microbiome composition to become more or less similar to their donor after FMT?
How quickly would you expect to see that change?
What do you expect would happen a few months or years after the transplant?

```{usage}
use.action(
    use.UsageAction('fmt', 'cc'),
    use.UsageInputs(
        diversity_measure=core_metrics_results.jaccard_distance_matrix,
        metadata=sample_metadata,
        distance_to='donor',
        compare='baseline',
        time_column='timepoints',
        reference_column='DonorSampleID',
        subject_column='PatientID',
        filter_missing_references=True,
        against_group='0',
        p_val_approx='asymptotic'
    ),
    use.UsageOutputNames(
        stats='jaccard_raincloud_stats',
        raincloud_plot='jaccard_raincloud_plot'
    )
)
```

In the resulting *raincloud plot*, we see that following cancer treatment (Timepoint 0-2) the recipients distance to donor is relatively high.This seems to be naturally resolving itself but after FMT intervention(Timepoint 3) the recipient's microbiome looks more similar to the donated microbiome. We see some stability in this but by the last time point the recipient's microbiome mostly looks unique from its donated microbiome.

### Community Richness

While FMT research is still in its infancy, some early studies suggest that long-term, the emergence of "personalize microbiome" is common after FMT, but in some cases features of the donated microbiome are retained.
Let's investigate this idea in an exercise.

(alpha-diversity-exercise)=
```{admonition} Tracking alpha diversity following FMT
:class: note

Above we tracked the *distance to donor* following FMT intervention.
This is a beta diversity approach to assessing Chimeric Asymmetric Community Coalescence that tells us about the overall membership or composition of a microbiome.
We can also directly investigate what happens to community richness, an alpha diversity metric, following FMT intervention. Retaining alpha diversity long term is an indicator of this emergence of "personalize microbiome" as opposed to reverting to their baseline. 

Try to adapt your above command, and perhaps referring to other sections of the tutorial, to compute a raincloud plot that will allow you to investigate this using the *observed features* vector you computed earlier.
When you're done, or if you get stuck, expand the following dropdown box for a command that will generate this plot for you.
```

(alpha-diversity-solution)=
```{admonition} Solution
:class: dropdown, note

```{usage}
use.action(
    use.UsageAction('fmt', 'cc'),
    use.UsageInputs(
        diversity_measure=core_metrics_results.observed_features_vector,
        metadata=sample_metadata,
        distance_to='donor',
        compare='baseline',
        time_column='timepoints',
        reference_column='DonorSampleID',
        subject_column='PatientID',
        filter_missing_references=True,
        against_group='0',
        p_val_approx='asymptotic'
    ),
    use.UsageOutputNames(
        stats='obs_features_raincloud_stats',
        raincloud_plot='obs_features_raincloud_plot'
    )
)
``` 
The community richness at each timepoint seems to be considerably variable and is never sigificantly higher than the baseline. It does look like after FMT intervention the community richness of the recipients is less variable but again it is not significantly higher than baseline. The last timepoint has the lowest community richness but there are only 3 subjects that were sampled that long after FMT. This makes it hard to tell if that is representative of the whole population or if those individuals where sampled that long after FMT because of complications. 

### Distance to Baseline
Developing a personalized microbiome following FMT intervention is expected. Although, there is no expected threshold for the length of time before personalization starts. It is important that the recipients microbiome doesn't return to their baseline composition. `qiime fmt cc` can help us investigate this! 

```{usage}
use.action(
    use.UsageAction('fmt', 'cc'),
    use.UsageInputs(
        diversity_measure=core_metrics_results.jaccard_distance_matrix,
        metadata=sample_metadata,
        distance_to='baseline',
        compare='baseline',
        time_column='timepoints',
        baseline_timepoint='0',
        subject_column='PatientID',
        filter_missing_references=True,
        against_group='1',
        p_val_approx='asymptotic'
    ),
    use.UsageOutputNames(
        stats='baseline_jaccard_stats',
        raincloud_plot='baseline_jaccard_raincloud_plot'
    )
)

```
Now lets look at the distance to their baseline. This will help us identify if the microbiome that emerges after FMT intervention is a unique personalized microbiome or if it is reverting to their baseline. Looking at this raincloud plot it looks like the microbiome always looks distinct from the baseline microbiome and that pattern continues all the way to the last timepoint. This probably indicates that the recipient microbiome following FMT is unique! 


### Proportional Engraftment of Donor Strains(PEDS)

If you are interested in how many microbes from the donor engrafted in the recipient, Proportional Engraftment of Donor Strains helps capture just that. 

The above metrics capture how similar the recipient and donor microbiomes are. However, We want these microbiome to coalesce asymmetrically meaning that we want the donor's features to be more prominment in the recipeint following FMT than baseline features. This metrics investigates this asymmetric colescence and captures how many donated features engrafted. Note that Drop incomplete timepoints will be run before drop incomplete subjects. 

```{usage}
sample_peds_dists, = use.action(
        use.UsageAction('fmt', 'sample_peds'),
        use.UsageInputs(
            table=core_metrics_results.rarefied_table,
            metadata=sample_metadata,
            time_column='timepoints',
            reference_column='DonorSampleID',
            filter_missing_references=True,
            subject_column='PatientID', 
            drop_incomplete_timepoints = ["2","5","6"],
            drop_incomplete_subjects=True
        ),
        use.UsageOutputNames(
            peds_dists='sample_peds_dist'
        )

    )
```
Now, we have our PEDS metrics and we want to visualize them. Lets use `qiime fmt peds-heatmap`. 

```{usage}
use.action(
    use.UsageAction('fmt', 'heatmap'),
    use.UsageInputs(
        data=sample_peds_dists,
    ),
    use.UsageOutputNames(
        visualization='sample-heatmap',

    )
)
```
This can also be run using `qiime fmt peds`. Run `qiime fmt ped --help` for more info! 

Looking at this output, timepoint 0 seems to have very low  proportional engraftment of donor strains. However, we can see that there is a a kind of bi-modal distribution at timepoint 1. Some of the samples have a relatively proportional engraftment of donor strains and some of them have a relatively low proportional engraftment of donor strains . At timepoint 3, we can see that there is less variation between recipients. This lines up pretty well with what we were seeing with the raincloud plot. This makes sense because they both are investigating Community Coalesense! 

Note: you can also use `qiime stats plot-rainclouds` to visualize this data! This will create plots more similar to the ouput of `qiime fmt cc` above. 

### Monte Carlo Simulation of PEDS

Now an important question to ask is: 
"Is the overlap between features due to the FMT or are the similarities between these 2 microbiome by random chance". 
That's where `peds-simulation` comes in! 

PEDS simulations is a Monte Carlo simulation that randomizes the relationships between donors and recipients, to test whether the PEDS score between a recipient and their actual donor is significantly higher than PEDS scores between other recipients paired with random donors.

Note that we are testing this on an equally amount of pre-fmt and post FMT samples and this will likely lead to a more conservative
global test results. 

Alright, Lets take a look 👀.

```{usage}
per_subject_stats, global_stats = use.action(
        use.UsageAction('fmt', 'peds_simulation'),
        use.UsageInputs(
            table=core_metrics_results.rarefied_table,
            metadata=sample_metadata,
            time_column='timepoints',
            reference_column='DonorSampleID',
            filter_missing_references=True,
            subject_column='PatientID', 
            drop_incomplete_timepoints = ["2","5","6"],
            drop_incomplete_subjects=True
        ),
        use.UsageOutputNames(
            per_subject_stats='per-subject-stats',
            global_stats='global-stats'
        )
    )
```
We can now re-vizualize our heatmap and include these new stats that we created. 

```{usage}
use.action(
    use.UsageAction('fmt', 'heatmap'),
    use.UsageInputs(
        data=sample_peds_dists,
        per_subject_stats=per_subject_stats,
        global_stats=global_stats
    ),
    use.UsageOutputNames(
        visualization='peds-stats-heatmap',

    )
)
```
We can see that there are many per-subject stats that where the simulated data(randomly paired recipients and donors) have higher PEDS than the true donor recipient pair but the majority of our comparisons are significant and globally are true pairs are significantly higher than our simulated donor recipient pairs. 

Thats good news! 


### Proportional Persistence of Recipient Strains(PPRS)

Proportional Persistence of Recipient Strains (PPRS) investigates if there are any features from the recipients baseline that stick around after FMT intervention.

It is a very similar investigation to using `fmt cc` with the distance-to parameter set to baseline. This will again help us evaluate if a uniquew personalized microbiome is emerging or if the microbiome is reverting back to baseline.

```{usage}
sample_pprs_dists, = use.action(
        use.UsageAction('fmt', 'sample_pprs'),
        use.UsageInputs(
            table=core_metrics_results.rarefied_table,
            metadata=sample_metadata,
            time_column='timepoints',
            baseline_timepoint='0',
            filter_missing_references=True,
            subject_column='PatientID', 
            drop_incomplete_timepoints = ["2","5","6"],
            drop_incomplete_subjects=True
        ),
        use.UsageOutputNames(
            pprs_dists='sample_pprs_dist'
        )

    )
```
This can also be viewed using `fmt heatmap`!

```{usage}
use.action(
    use.UsageAction('fmt', 'heatmap'),
    use.UsageInputs(
        data=sample_pprs_dists
    ),
    use.UsageOutputNames(
        visualization='pprs-heatmap',

    )
)
```
From this visualization we can see that there is a relatively low presentage of persistant recipient features. Generally that's a good sign of engraftment! 

However, It looks like FMT.0035 has many persistant recipient features by timepoint 4. For an extra challenge go through our previous visualizations and see if you can identify any other signs of low engraftment extent! 

## Donated Microbiome Indicator Features

Another method that researchers commonly use to assess engraftment is donated microbiome indicator features

Currently, we do this by running ANCOMBC comparing the recipient at baseline to their donor. However because ancombc can not be run on dependent samples (i.e. a patient over time), we are not able to compare the baseline recipient to their donor because in this case they are the same patient. 

However, some people investigate microbes that other papers reported as imported. We will investigate if the feature was in the donated microbiome and if the recipient recieved the feature. 

(This review)[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8436969/] reported that the class of *Clostridia* is correlated with positive health outcomes. 

So here we are going to use q2-longitudinal to track the relative abunadnce of *Clostridia* over time. This will help us identify if the donated microbiome had the previously reported feature and if the recipients recieve this feature.

```{usage}

taxonomy_url = 'https://qiime2-workshops.s3.us-west-2.amazonaws.com/itn-aug2024/taxonomy.qza'

taxonomy = use.init_artifact_from_url('taxonomy', taxonomy_url)
```
Lets use this Taxonomy to collapse our table at level three or Class. 

```{usage}
collapsed3_table,  = use.action(
    use.UsageAction('taxa', 'collapse'),
    use.UsageInputs(
        table=autofmt_table,
        taxonomy=taxonomy,
        level=3
    ),
    use.UsageOutputNames(
        collapsed_table='collapsed3-table.qza'
    )
)
```

Let's transform this feature-table into relative frequency feature table. This will allow us to look at relative frequency instead of raw counts that can be misleading. 

```{usage}
collapsed3_table_rf,  = use.action(
    use.UsageAction('feature_table', 'relative_frequency'),
    use.UsageInputs(
        table=collapsed3_table
    ),
    use.UsageOutputNames(
        relative_frequency_table='collapsed3-table-rf.qza'
    )
)
```
Now we are ready to use `q2-longitudinal` to track our previously identified feature

### Tracking Features with `q2-longitudinal`

```{usage}
md_collapsed3_table_rf= use.view_as_metadata(
    'md_collapsed3_table_rf',
     collapsed3_table_rf)

merged_collaped3_rf = use.merge_metadata(
    'merged_collaped3_rf',
    sample_metadata,
    md_collapsed3_table_rf)

use.action(
    use.UsageAction('longitudinal', 'linear_mixed_effects'),
    use.UsageInputs(metadata=merged_collaped3_rf,
                    state_column='day-relative-to-fmt',
                    group_columns='autoFmtGroup',
                    individual_id_column='PatientID',
                    metric='k__Bacteria;p__Firmicutes;c__Clostridia'),
    use.UsageOutputNames(visualization='lme-clostridia-treatmentVScontrol')
)
```
Looking at this it looks like both the control and FMT groups had Clostrida in relatively low abundances. Both Groups seems to have in increase in Clostrida over time, which is good becuase we know thats correlated with postive treament outcomes. However, there is no significant change between groups. This suggests that Clostrida is not a great donor indicator for this group of patients. 

Check back soon for more updates on tracking donated microbiome indicator features using ANCOMBC2 which will allow for repeated sampling! 

### Feature Engraftment

Feature engraftment doesn't "assess engraftment" of a recipient but instead investigates if there are features that engraft across all subjects. This will help researchers understand which features are sucessfully engrafting and help them decide which successfully engrafted features correlate with postive clinical outcome. 

This could be utilized to help researchers decide on specified communities to donate (as opposed the black-box commmunity approach we have currently). This specified communities would be full of microbes that have successfully engraftment accross patients and are associated with postive clincal outcomes.

Since we previously looked at the Clostrida, lets continue look at how our classes engraft accross subjects. 

Let's take a look now and see if we have any features that are sucessfully engrafting accross subjects.

```{usage}
feature_peds_dist, = use.action(
    use.UsageAction('fmt', 'feature_peds'),
    use.UsageInputs(
        table=collapsed3_table,
        metadata=sample_metadata,
        time_column='timepoints',
        reference_column='DonorSampleID',
        subject_column='PatientID',
        filter_missing_references=True
    ),
    use.UsageOutputNames(
        peds_dists='feature_peds_dist'
    )

)
```
And again, Visualize the data using `qiime fmt heatmap`. 
We are using the `level-delimiter` parameter here so that the taxonomic strings will only show the lowest relevant taxonomic level! 

```{usage}
use.action(
    use.UsageAction('fmt', 'heatmap'),
    use.UsageInputs(
        data=feature_peds_dist,
        level_delimiter = ';',
    ),
    use.UsageOutputNames(
        visualization='feature-heatmap',

    )
)

```
Oh look! Clostridia seems to be a feature that consistantly engrafts from our donor! It seems like it doesnt engraft "better" than spontanous recover but it is a consistant engrafter. Firmicutes also seems to pretty consistently engraft but it is important to note that our N for this study is quite small. Bacilli seems to engraft sucessfully too but its also present before the FMT so its hard to tell if thats from the donor or just a common gut 🐛 bug! 

Overall in this study it seems that we have Asymetric Chimeric Community Coalesense. We also have stability in our shift towards the donated microbiome as there are no signs that our subjects are reverting back to baseline. From these factors I would conclude that we have relatively high engraftment extent and clincal findings are probably due to FMT engraftment. I would love to see more donated microbiome indicator species to further support this claim! 
