# Plot Cox regression results on patient-level CTN presence

Visualizes the prognostic value and cellular composition of Cell-type
Triad Niches (CTNs).

The left panel features a lollipop plot where the x-axis represents the
log hazard ratio (HR) from Cox regression models, with lines colored by
Benjamini-Hochberg adjusted P-values (red: FDR \< 0.05; black: 0.05
\\\le\\ FDR \< 0.1).

The right panel features an aligned bar plot illustrating the cell type
composition for each CTN, with core cell types are highlighted with
black boxes.

[ggtext](https://wilkelab.org/ggtext/reference/ggtext.html),
[forcats](https://forcats.tidyverse.org/reference/forcats-package.html),
[cowplot](https://wilkelab.org/cowplot/reference/cowplot-package.html),
and
[patchwork](https://patchwork.data-imaginist.com/reference/patchwork-package.html)
packages are required for the function.

## Usage

``` r
draw_CTN_coxph_logHR(
  coxph_res,
  CTN_merged_celltype_prop,
  top = NULL,
  alpha = 0.1,
  celltype_levs,
  celltype_palette,
  dot_fontsize = 14,
  rename_celltype = NULL,
  line_width = 0.4,
  y_text_size = NULL,
  width_ratio = c(0.5, 0.5)
)
```

## Arguments

- coxph_res:

  A data frame representing the Cox proportional hazards regression
  results on patient-level CTN presence. Output of the
  [`CTN_presence_coxph()`](https://github.com/YufanWu147/SCOPE.CTN/reference/CTN_presence_coxph.md)
  function

- CTN_merged_celltype_prop:

  A data frame representing the cell type proportion for each CTN with
  the following columns: `CTN` (the CTN names), `Celltype` (cell type
  labels), and `prop` (cell type proportion of the CTN).

- top:

  An integer specifying the maximum number of top-ranked CTNs ordered by
  \\\|logHR\|\\ to display on the plot.

- alpha:

  A numeric value specifying the False Discovery Rate (FDR) significance
  threshold for filtering the CTNs.

- celltype_levs:

  A character vector defining the levels or ordering of the cell types.

- celltype_palette:

  Input of the
  [`color_CTN_names()`](https://github.com/YufanWu147/SCOPE.CTN/reference/color_CTN_names.md)
  function. A named vector with colors as values and cell type labels as
  names. Must contain all cell types.

- dot_fontsize:

  Input of the
  [`color_CTN_names()`](https://github.com/YufanWu147/SCOPE.CTN/reference/color_CTN_names.md)
  function. Size of the dots. Default is 16pt.

- rename_celltype:

  A function for renaming cell type labels.

- line_width:

  A numeric value specifying the line width for the boxes highlighting
  the core cell types.

- y_text_size:

  A numeric value specifying the font size of the y-axis labels.

- width_ratio:

  A numeric vector of length 2 specifying the width ratios between the
  log harzard ratio lollipop plot and cell type proportion bar plot.

## Value

A combined
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object containing both the survival lollipop plot and the cell type
composition bar plot.
