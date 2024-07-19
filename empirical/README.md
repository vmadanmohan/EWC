# Empirical analysis
This directory contains codes to infer inter-ROI communication from source-localised MEG recordings using the partial correlation (EWC) and bivariate transfer entropy (Full).

### Pre-requisites
Source localised MEG recordings in ```.csv``` format, ```.txt``` file containing inter-regional conduction delays (in seconds) stored as an $N \times N$ matrix, and ```.txt``` file containing the Start and End timepoints of bad segments (in seconds)
Ensure that the delay and badsegment files conform to the formats given in ```example_delay.txt``` and ```example_badseg.txt```.
## File overview
```cleanup.py``` - takes in source localised time series (.csv), inter-regional delays (.txt) and bad segment intervals (.txt). The main data is orthogonalised, epoched into 10s segments, and epochs containing bad segments are removed. The data is then trimmed to remove possible processing artifacts and is then stored in a ```.mat``` file as a $Time \times ROIs$ matrix.

```PC_EWC_empirical.m``` - takes in the ```.mat``` file containing the orthogonalised regional time series, and estimates PC-EWC (using the function in ```EWC/functions```).

```PC_TE_comparison.m``` - Computes the hemisphere-wise correlation between the communication patterns inferred using PC-EWC and TE-Full, for all subjects.

```TE_Full_empirical.m``` - takes in the ```.mat``` file containing the orthogonalised regional time series, and estimates TE-Full (using the function in ```EWC/functions```).

```mapping.txt``` - the ordering of the ROIs post-source localisation in Brainstorm. Change this mapping as required to suit your data. Re-order the matrices using ```matrestruct``` (for square matrices) in ```EWC/functions```.

```subIDs.txt``` - contains the IDs of the HCP subjects. Change to suit your data.
