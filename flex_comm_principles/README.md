# Principles of flexible neural communication
This directory contains scripts used in "Assessment of oscillatory mechanisms underlying flexible neural communication in the human brain" 2025 by Varun Madan Mohan, Thomas F Varley, Robin F H Cash, Caio Seguin, and Andrew Zalesky.

### Pre-requisites

#### Software / packages
```MATLAB2022b``` with the [Signal Processing Toolbox](https://au.mathworks.com/products/signal.html) for spectral analyses.
[Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/home) - contains function ```threshold_proportional.m``` used to threshold the structural connectivity matrix.
[parc_plotter](https://github.com/faskowit/parc_plotter) - used in ```visualise.m``` to project correlation coefficients onto the cortical surface to visualise heterogeneous neural oscillatory dependence.

#### Data
1. A .mat file containing Source-localised MEG time series for each subject (variable name ```main_data```), and an Nregion x Nregion matrix of inter-regional delays (units of timesteps - e.g. for a sampling rate of 1000Hz, a latency of 10ms between regions _i_ and _j_ would be _delay(i,j)=10_ - variable name ```delay```) (.mat) ```<SUB_IDX>_resting.mat``` (SUB_IDX numbered from 1 to Nsub)
2. Group-normative structural connectivity matrix (thresholded to 15%, inter-hemispheric connections pruned and binarised in code) (.txt) - ```group_SC.txt```

### Sub-directory descriptions:

```observed/``` contains scripts to estimate the correlation between communication inferred using EWC and target power/PLV. (Note: Intersite Phase Clustering (ISPC) which is used in the main text, is an alternative term for the Phase Locking Value (PLV)).

```cyclic_surr/``` contains scripts that carry out similar estimations as in ```observed/```, but for cyclically permuted surrogates

```surr_correction.m``` compares the observed and surrogate results, performs a surrogate-derived false discovery rate corrected significance threshold on the results

```eventiden.m``` finds the timepoints of significant "communication events", and is used for the estimation of EWC and neural oscillatory measures in both observed and surrogate analyses.

```visualise.m``` loads surrogate-corrected results from ```surr_correction.m``` and projects correlation values onto the cortical surface using ```parc_plotter```.

```masterLoop.sh``` computes the observed communication principles for _Nsub_ subjects, runs _M_ surrogates per subject, performs a surrogate-based correction of observed EWC-power/ISPC relationships, and projects the mean correlation coefficients onto the cortical surface for visualisation.
