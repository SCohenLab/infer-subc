# Example input for a volcano plot explorer

These fictional files contain three conditions, six independent samples per
condition, four measurements, and two comparisons against control.

Open `../test.html` in a browser, select the two CSV files, and click
**Load selected files**. Choose a comparison and click a point in the volcano
plot to see its sample values. The first feature is displayed initially.
The files are read locally; Plotly and the CSV parser load from CDNs and
therefore require an internet connection.

Alternatively, run `py -m http.server 8000` from the directory containing
`test.html`, open `http://localhost:8000/test.html`, and click **Load bundled
example**. Browsers do not allow that button to fetch local files when the HTML
is opened directly, so use the file selectors in that case.

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
  volcano x-axis. A positive value means the comparison mean is higher.
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

The explorer uses `-log10(p_adjusted_bh)` for its y-axis and identifies
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
