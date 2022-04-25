# CellClock-rep :: Analysis and inference on Cell Cycle kinetics: renewal rate and phase distribution
Computational methods to analyze statistics of experimental cell frequencies in distinct cell-cycle phases and infer cell-cycle progression in the esophageal epithelium of adult mice from synchronized S-phase cells from dual EdU-BrdU chase experiments.

This repository contains code that was originally used in the following manuscript:
  > Abby E, Dentro SC, Hall MWJ, et al. (2022) Notch1 mutation drives clonal expansion in normal esophageal epithelium but impairs tumor growth. _under revision_

### Overview
The repository divides into code to perform the statistical analysis of experimental data in WT and _Notch1_ mutant animals (R language) and stochastic simulation-based inference on cell cycle kinetics, renewal rate and changes in cell-cycle phase distributions (Matlab).

The simulation algorithms to infer WT and mutant cell behavior are based on the **Single Progenitor (SP) model** paradigm of cell renewal in murine squamous epithelium. It is assumed that both WT and mutant cell dynamics respond to a single population of progenitor cells renewing with a given division rate (_Lambda_) defining a relatively short cell-cycle period, and division daughters committing to stochastic fates. The distribution of individual cell-cycle periods inferred for WT progenitors in Piedrafita et al (2020) is here used to simulate a random time a progenitor cell stays in each of the cell-cycle phases according to the relative predominance with which that phase is found across the entire overall basal cell population.

The Monte-Carlo simulator of progenitor cell-cycle phase dynamics is extended to incorporate differentiation and stratification processes, allowing to study renewal and stratification dynamics concomittantly.

### Main scripts
- **StatisticalAnalysis-CellCyclePhases.R** : main R script to reproduce the statistical analysis of data from the EdU-BrdU chase experiment. It reads the original data in _Summary-EdU-BrdU-data.xlsx_, builds a multinomial logistic regression model and shows the contingency analysis of cell-cycle phase distributions and statistics on cell renewal and stratification estimates.
- **Inference-CellCyclePhases.m** : Matlab script used for cell-cycle phase inference, implementing a simple simulation algorithm of individual progenitor cell progression through different cell cycle phases. It infers cell-cycle phase distributions for EdU-labelled populations at different times post-EdU labelling, assuming all basal cells behave as progenitors and considering as input the inferred distribution of cell cycle periods for normal progenitor cells in the esophagus (Piedrafita et al, 2020) to model the time in different phase as a fraction of a random total cycle time period.
- **Inference-CellCyclePhases-RenewalStratif.m** : Matlab script used for cell-cycle phase inference. This version considers a more realistic description of cellular dynamics including stochastic division fate outcomes and basal differentiating progeny that transiently contributes to G0/G1 frequencies. Apart from providing inference on cell-cycle phase distributions over time, it allows to simulate and fit cell renewal and stratification rate time courses.

### Dependencies
- **MCsimulator-dynamics-EdU-SP-total.m** : Matlab function implementing a Monte Carlo simulator of EdU-labelled cell clone dynamics under the SP model with balanced fates (time-dependent division and shedding rates are assumed realistically). As time progresses, simulations record both individual cell-cycle state (whether progenitor in G1, S, G2 or M phase or differentiating in G0) and location (whether basal or suprabasal) starting from an initial population of cells synchronized in S-phase (EdU stained).
- **Summary-EdU-BrdU-data.xlsx** : Experimental data on numbers of cells in different cell-cycle phases (stained with different cell-cycle phase markers) from the EdU-BrdU chase experiment. (Imported by the R script)
- **freq_CellCyclePhases.xlsx** : Summary of cell frequencies and summary statistics from the EdU-BrdU chase experiment. (Imported by the Matlab script)

- `Plot-Rscripts` folder : Set of functions used for generating R plots called from _StatisticalAnalysis-CellCyclePhases.R_ (mostly ggplot).


### Requirements
- R version 4.0.3 (2020-10-10)
- Matlab R2020b
- (R packages loaded in StatisticalAnalysis_CellCyclePhases.R)

