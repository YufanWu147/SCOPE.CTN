# Generate a Thin Plate Spline Basis Matrix

Internal function to converts a spatial neighborhood cell type
proportion matrix into a thin plate regression spline basis matrix.

## Usage

``` r
apply_smoothing(dat, mgcv_df)
```

## Arguments

- dat:

  Neighborhood cell type proportion matrix

- mgcv_df:

  The dimension of the basis used to represent the thin plate regression
  spline smooth term. The default is 15. See
  [mgcv::s](https://rdrr.io/pkg/mgcv/man/s.html) for details.

## Value

A numeric matrix containing the smoothed cell type proportions, with row
names preserved from the original input.
