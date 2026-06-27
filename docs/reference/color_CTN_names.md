# Separating cell types in CTN labels with colored dots

Adding colored dots before each cell type in cell type triad (CTN)
labels. Outputs can be used in titles and labels `ggplot` with markdown
text enabled using
[ggtext::element_markdown](https://wilkelab.org/ggtext/reference/element_markdown.html).

## Usage

``` r
color_CTN_names(
  CTN,
  celltype_palette = c(B = "#fee440", CD4T = "#3a86ff", CD8T = "#8338ec"),
  dot_fontsize = 16,
  sep = "_"
)
```

## Arguments

- CTN:

  CTN label; e.g. B_CD4T_CD8T.

- celltype_palette:

  A named vector with colors as values and cell type labels as names.
  Must contain all cell types.

- dot_fontsize:

  Size of the dots. Default is 16pt.

- sep:

  A character string separating the cell types in CTN. Note that it
  can't appear in cell type labels.

## Examples

``` r
library(ggplot2)
library(ggtext)

CTN_colored <- color_CTN_names(
  CTN = "B_CD4T_CD8T",
  celltype_palette = c("B" = "#fee440", "CD4T" = "#3a86ff", "CD8T" = "#8338ec"))

# ggplot() +
# labs(title = CTN_colored) +
#  theme(plot.title = element_markdown(hjust = 0.5))
```
