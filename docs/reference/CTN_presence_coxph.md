# Cox regression on patient survival based on CTN presence

Evaluates the prognostic ability of Celltype Triad Niches (CTNs) by
fitting a Cox proportional hazards regression model using patient-level
CTN presence as the main predictor.

## Usage

``` r
CTN_presence_coxph(
  Full_CTN_table,
  CTN_ls,
  clinical_df,
  time = "OS",
  event = "Ev.O",
  patient_id = "Patient_ID",
  confounder = NULL,
  min.n.patients = 3,
  min.n.cells = 100
)
```

## Arguments

- Full_CTN_table:

  A data frame containing clustering labels for all cell-type triads

- CTN_ls:

  A character vector containing the names of the CTNs to evaluate.

- clinical_df:

  A data frame containing image-level clinical data. It must include
  columns for patient ID (`patient_id`), survival time (`time`), status
  event (`event`), and any specified confounders (`confounder`).

- time:

  A character string specifying the survival time column name. See
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) for
  details.

- event:

  A character string specifying the survival status indicator column
  name, normally 0=alive, 1=dead. See
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) for
  details.

- patient_id:

  A character string specifying the patient ID column name.

- confounder:

  A character string or vector specifying column names of confounders to
  be adjusted for in the model. Defaults to `NULL`.

- min.n.patients:

  An integer specifying the minimum patient threshold. Only CTNs that
  are both present present and absent in \\\ge\\ `min.n.patients` will
  be evaluated. Defaults to 3.

- min.n.cells:

  An integer specifying the minimum cell count threshold required to
  define a CTN as present on a tissue section. Defaults to 100.

## Value

A data frame containing the survival analysis results for each CTN:

- CTN:

  The name of the CTN.

- logHR:

  Log harzard ratio.

- HR:

  Hazard ratio.

- pval:

  P-value

- c_index:

  Concordance-index (c-index)

- FDR:

  The Benjamini-Hochberg adjusted p-value

## Details

Patient-level presence of a CTN is defined as having \\\ge\\
`min.n.cells` cells within that niche in at least one tissue section.
