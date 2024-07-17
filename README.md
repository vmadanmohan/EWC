# Event-marked Windowed Communication (EWC)
This repo contains the codes used in carrying out the analyses in "An empirical assessment of neural-oscillatory mechanisms underlying flexible communication". The Event-marked Windowed Communication (EWC) is an implementation to gauge directed interactions from regional time-series, which can then be used to infer communication between the regions. The EWC can in principle be estimated using any symmetric measure of statistical dependence - We use Partial Correlation (PC), Conditional Mutual Information (cMI), and bivariate Transfer Entropy (TE). 
cMI and TE were estimated using the [Java Information Dynamics Toolkit (JIDT)](jlizier.github.io/jidt/).

Python packages used:
```numpy, scipy, osl, pandas```
MATLAB version: R2022b

### Workflow for 1 subject
```emp_comm.py``` $\rightarrow$ (```surr_comm.py```) $\times$ M $\rightarrow$ ```surr_corr.py``` $\rightarrow$ ```symmetry.py, flexibility.py``` $\rightarrow$ ```EWC_SC_analysis.m, pow_coh_analysis.m```

## File overview

```subIDs.txt``` - contains the IDs of the HCP subjects. Change to suit your data.

```mapping.txt``` - the ordering of the ROIs post-source localisation in Brainstorm. Change this mapping as required to suit your data. Re-order the matrices using ```matrestruct()``` (for square matrices) or ```rrestruct()``` (for the principle matrices) in ```miscfunc.py``` prior to cortical projection.

### ```LSM/```

### ```empirical/```

### ```functions/```
