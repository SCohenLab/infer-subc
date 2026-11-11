# Example input for a volcano plot explorer

These fictional files contain three conditions, six independent samples per
condition, four measurements, and two comparisons against control.

Open `../test.html` in a browser and select **Rhoads, et al., 2025 data** or
**Fictional example** from the **Sample dataset** dropdown. The sample files
download directly from the infer-subc GitHub repository, including when opening
the HTML file directly. The sample URLs are pinned to repository revision
`4195e1dfc418b0babe7b1a9c2c3f0f53f150a4b0` to keep the files consistent;
update that revision in the HTML if the published sample data changes.
An internet connection is required for the sample files and the plotting, CSV,
and workbook libraries. Use **Reload sample** to retry a failed download.

Choose a comparison and click a point in the volcano plot to see individual
cell values. The first feature is displayed initially. Use the metric selector
to reach points that overlap on the volcano plot.

To load local files instead, expand **Load local files instead**. For the
study, select the Table S1 workbook, `20240722_NvA-untr_fp-metrics.csv` for
neuron vs astrocyte, and `20250110_N_uncVcon_ctrls_fp-metrics.csv` for neuron
vs astrocyte-naive neuron, then click **Load study files**. For the fictional
CSV format, select `sample_measurements.csv` and `sample_statistics.csv` and
click **Load selected files**. Local files are processed in the browser.

The study loader is specific to the Table S1 layout and these two CSVs:
it uses `cell_type` (Neuron/Astrocyte) in the first CSV and
`media_condition` (Conditioned/Unconditioned) in the second. The conditioned
group corresponds to the neuron control and the unconditioned group to the
astrocyte-naive neuron. `dataset` and `image_name` identify individual cells.
The plot uses the workbook's mean rank difference as reported, with neuron
control minus the other group, and its reported q-value. An upward triangle
plots a `<` q-value at its bound, so the true -log10(q) is higher; a downward
triangle plots a `>` q-value at its bound, so the true -log10(q) is lower.
Hover to see the original bound. Missing q-values or rank differences are
excluded and counted. Missing, NA, NaN and infinite cell values are excluded
from the bars and dots, with per-group exclusion counts shown on the chart.
Colors use the workbook's significance flags, not recalculated statistics.

## Sample measurements

`sample_measurements.csv` has one row per independent sample:

- `sample_id`: unique identifier for a sample.
- `condition`: experimental condition, matching the names in the statistics file.
- Other columns: numeric measurements, with one column per feature. The column
  name must exactly match the `feature` field in the statistics file.

Use descriptive, unique feature names and consistent units for each feature.
For missing measurements, leave the cell empty, not zero; record sample counts
per feature after excluding missing values. If rows are cells or technical
replicates rather than independent biological samples, include biological sample
identifiers and use the analysis appropriate to the experimental design.

## Statistics

`sample_statistics.csv` has one row for each feature and comparison. A
comparison is identified by `reference_condition` and `comparison_condition`;
each combination of these and `feature` must be unique.

- `n_reference`, `n_comparison`: non-missing sample counts.
- `mean_reference`, `mean_comparison`: means of the corresponding measurements.
- `log2_fold_change`: `log2(mean_comparison / mean_reference)`, used for the
  illustrative effect size. A positive value means the comparison mean is higher.
- `mean_rank_difference`: mean rank of reference minus mean rank of comparison,
  used for the volcano x-axis.
- `u_statistic`: Mann-Whitney U for the comparison condition (the first
  argument to the test).
- `p_value`: two-sided Mann-Whitney U p-value.
- `p_adjusted_bh`: Benjamini-Hochberg adjusted p-value.
- `significant_fdr_0_10`: `true` when the adjusted p-value is at most 0.10.

For these examples, p-values use SciPy's asymptotic Mann-Whitney test with tie
and continuity corrections. Benjamini-Hochberg adjustment is applied **separately
within each comparison**, across all four features. The sample sizes, means,
test statistics and p-values were calculated from the measurements file; numeric
values are rounded for display.

For the fictional CSV format, the explorer uses `-log10(p_adjusted_bh)` for its y-axis and identifies
significant points with a horizontal line at `-log10(0.10) = 1`. Hover text
shows the selected comparison and both p-values. Bars display means calculated
from the loaded samples, not the means in the statistics file. The explorer uses
the supplied statistics without recalculating the tests or FDR adjustment.

The log2 mean ratio is an illustrative effect size, not an output of the
Mann-Whitney test. It is defined here because both means are positive. If your
measurements can be zero or negative, use a scientifically justified effect
size and document its definition; do not silently add an arbitrary pseudocount.
The explorer rejects an adjusted p-value of zero rather than silently replacing
it for plotting. Resolve such values in the source analysis before loading.

For real data, preserve the actual test method, alternative hypothesis,
multiple-testing method and adjustment family used by your existing analysis;
do not recompute adjustments using only displayed or filtered features.
