# Event-marked Windowed Communication (EWC)
This repo contains the codes used in carrying out the analyses in "Event-marked Windowed Communication: Inferring activity propagation patterns from neural time series". The Event-marked Windowed Communication (EWC) is an implementation to gauge directed interactions from regional time-series, which can then be used to infer communication between the regions. The EWC can, in principle, be estimated using any symmetric measure of statistical dependence - We use Partial Correlation (PC), Conditional Mutual Information (cMI), and bivariate Transfer Entropy (TE). 
cMI and TE were estimated using the [Java Information Dynamics Toolkit (JIDT)](jlizier.github.io/jidt/).

MATLAB version: R2022b \\
Python packages used:
```numpy, scipy, osl, pandas```

Codes are organised into 3 directories: 

```LSM/``` - contains codes to model an simple network-motif with Linear stochastic model dynamics, and Poisson firing sources. Corresponds to the section "Asymmetric signalling over a network motif".

```empirical/``` - contains codes to gauge computational tractability, and empirical results. Corresponds to sections "Computational tractability of the EWC protocol" and "Inferring whole-brain interaction patterns from MEG recordings"

```functions/``` - contains all the functions used in the codes. Ensure that the path to this directory is included in all codes.
