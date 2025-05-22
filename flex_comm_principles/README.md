# Principles of flexible neural communication
This directory contains scripts used in "Assessment of oscillatory mechanisms underlying flexible neural communication in the human brain" 2025 by Varun Madan Mohan, Thomas F Varley, Robin F H Cash, Caio Seguin, and Andrew Zalesky.

Sub-directory descriptions:

```observed/``` contains scripts to estimate the correlation between communication inferred using EWC and target power/PLV.

```cyclic_surr/``` contains scripts that carry out similar estimations as in ```observed/```, but for cyclically permuted surrogates

```surr_correction.m``` compares the observed and surrogate results, performs a surrogate-derived false discovery rate corrected significance threshold on the results

```eventiden.m``` finds the timepoints of significant "communication events", and is used for the estimation of EWC and neural oscillatory measures in both observed and surrogate analyses.
