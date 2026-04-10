# Post-hoc Knockoffs

This repository contains the official code accompanying the paper:

> **Choosing the nominal level post-hoc with knockoffs using e-values**  
> Lasse Fischer and Konstantinos Sechidis  
> [arXiv:2511.11166](https://arxiv.org/abs/2511.11166)

## Overview

The knockoff filter is a powerful framework for controlled variable selection with false discovery rate (FDR) guarantees. While widely used in high-dimensional settings, its performance can deteriorate in low-dimensional scenarios with sparse signals. A key limitation is that the knockoff filter requires a minimum number of selections (at least 1/α) determined by the nominal FDR level.

In this work, we introduce a **post-hoc adjustment using e-values** that allows the nominal FDR level to be adapted *after* running the knockoff procedure. This leads to two practical benefits:

- **If the standard knockoff filter makes no discoveries:** we can safely increase the nominal level to potentially recover meaningful signals.
- **If the standard knockoff filter makes discoveries:** we can often reduce the nominal level to improve precision.

These improvements come at no cost — the post-hoc results are always at least as informative as those from the original knockoff procedure. We also extend the method to **derandomised knockoff** techniques for both FDR and PFER control.

## Repository Structure

```
posthoc_knockoffs/
├── src/
│   ├── utils.R              # Core implementation of all methods
│   ├── sysdata.rda          # Precomputed data for the case study
│   ├── Section_5_3_1/       # Simulations: §5.3.1 (Figures 2–3)
│   ├── Section_5_3_2/       # Simulations: §5.3.2 (Figure 4)
│   ├── Section_5_3_3/       # Simulations: §5.3.3 (Figure 5)
│   ├── Section_5_3_4/       # Simulations: Suppl. §S2 (Figure S1)
│   └── Section_6_1/         # Case study: Suppl. §S3 (Figure S2)
├── figures/                  # Pre-generated PDF figures from the paper
├── LICENSE
├── posthoc_knockoffs.Rproj
└── README.md
```

### Core Methods (`src/utils.R`)

The file `utils.R` implements the following methods:

**Standard knockoff procedures:**

- `KO_FDR(W, FDR_nominal)` — Original BC-knockoff filter for FDR control (Candès et al., 2018).
- `KO_select_PFER(W, PFER_nominal)` — Knockoff filter for per-family error rate (PFER) control (Janson and Su, 2016).

**Post-hoc knockoff procedures (our contributions):**

- `KO_FDR_posthoc(W, alpha_nominal)` — Post-hoc knockoff filter for FDR control (**Algorithm 1** in the paper). Returns a data-dependent α and the rejection set.
- `KO_FDR_derand_posthoc(W, FDR_nominal, alpha_ebh)` — derandomised post-hoc knockoffs for FDR control (**Algorithm 2**). Takes a matrix of knockoff statistics from multiple runs.
- `KO_PFER_derand_posthoc(W, PFER_nominal)` — derandomised post-hoc knockoffs for PFER control (**Algorithm S1** in the Supplementary Material). Adaptively selects the threshold η from the data.
- `KO_FDR_derand_posthoc_given_rejection_set(W, FDR_nominal, rejection_set)` — Computes the smallest post-hoc α for a user-specified rejection set (**Proposition 4**).

**derandomised knockoff baselines:**

- `KO_FDR_derand(W, FDR_nominal, alpha_ebh)` — Original derandomised knockoffs for FDR control (Ren and Barber, 2024).
- `KO_PFER_derand(W, PFER_nominal, eta)` — Original derandomised knockoffs for PFER control (Ren et al., 2023).

**Simulation utilities:**

- `gene_X(X_type, n, p, ...)` — Generates design matrices with various correlation structures (AR, block, IID, etc.) used in the simulation studies.
- `evenly_spaced_indices(p, p_relevant)` — Returns evenly spaced indices for placing non-null variables in the coefficient vector β.
- `plot_selection_heatmap(R_mat, feature_names, ...)` — Produces variable selection heatmaps (as in Supplementary Figure S2).

### Simulation Scripts (`src/Section_*`)

Each subfolder corresponds to a subsection of the paper and contains the R script(s) needed to reproduce the corresponding figures:

| Subfolder | Paper Section | Description | Figures |
|-----------|--------------|-------------|---------|
| `Section_5_3_1/` | §5.3.1 | Post-hoc KO (Alg. 1) vs. original knockoffs (Candès et al., 2018) | Figures 2, 3 |
| `Section_5_3_2/` | §5.3.2 | Post-hoc KO (Alg. 1) vs. calibration knockoffs (Luo et al., 2025) | Figure 4 |
| `Section_5_3_3/` | §5.3.3 | derandomised post-hoc KO for FDR (Alg. 2) vs. Ren & Barber (2024) | Figure 5 |
| `Section_5_3_4/` | Suppl. §S2 | derandomised post-hoc KO for PFER (Alg. S1) vs. Ren et al. (2023) | Figure S1 |
| `Section_6_1/` | Suppl. §S3 | Prognostic variable selection on synthetic clinical trial data (benchtm) | Figure S2 |

**To reproduce a figure:**

1. Navigate to the relevant subfolder under `src/`.
2. Open the R script and read the **header comments** — they describe the full process for generating the figure, including parameters and configuration.
3. Run the script in R.

### Pre-generated Figures (`figures/`)

The `figures/` directory contains all pre-generated PDF figures included in the paper and Supplementary Material.

## Getting Started

### Prerequisites

- **R** (version ≥ 4.0 recommended)
- Required R packages (install from CRAN):

```r
install.packages(c("knockoff", "ggplot2", "dplyr", "tidyverse", "pheatmap"))
```

- The [`knockofftools`](https://github.com/Novartis/knockofftools) package (install from GitHub):

```r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("Novartis/knockofftools")
```

- For the case study in Supplementary Section S3, the [`benchtm`](https://github.com/openpharma/benchtm) package is also needed:

```r
devtools::install_github("openpharma/benchtm")
```

> **Note:** Check the header of each script for any additional package dependencies specific to that analysis.

### Running the Code

1. **Clone the repository:**

   ```bash
   git clone https://github.com/sechidis/posthoc_knockoffs.git
   cd posthoc_knockoffs
   ```

2. **Open the R project** in RStudio (`posthoc_knockoffs.Rproj`), or set the repository root as your working directory:

   ```r
   setwd("/path/to/posthoc_knockoffs")
   ```

3. **Source the utility functions** to use the methods in your own analysis:

   ```r
   source("src/utils.R")
   ```

4. **Reproduce a figure** by navigating to the relevant subfolder under `src/` and running the R script. For example, to reproduce the comparison with original knockoffs (Figures 2–3):

   ```r
   # See header comments in the script for full instructions
   source("src/Section_5_3_1/<script_name>.R")
   ```

5. **Run the simulated clinical trial case study** (Supplementary Figure S2):

   ```r
   source("src/Section_6_1/<script_name>.R")
   ```

## Citation

If you use this code in your research, please cite our paper:

```bibtex
@misc{fischer2026choosingnominallevelposthoc,
      title={Choosing the nominal level post-hoc with knockoffs using e-values}, 
      author={Lasse Fischer and Konstantinos Sechidis},
      year={2026},
      eprint={2511.11166},
      archivePrefix={arXiv},
      primaryClass={stat.ME},
      url={https://arxiv.org/abs/2511.11166}, 
}
```

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Contact

- **Lasse Fischer** — [fischer1@uni-bremen.de](mailto:fischer1@uni-bremen.de), Competence Center for Clinical Trials Bremen, University of Bremen
- **Konstantinos Sechidis** — [kostas.sechidis@novartis.com](mailto:kostas.sechidis@novartis.com), Advanced Methodology and Data Science, Novartis Pharma AG
