# Functions
Functions used to estimate activity propagation using the EWC and conventional (Full) implementations. Ensure that path to this directory is included in all analysis codes.

### File overview

```CondMutInfo.m``` - computes the conditional Mutual Information (cMI) between all pairs of regions, using the EWC implementation (cMI-EWC).

```MIFull.m``` - computes the Mutual Information (MI) between all pairs of regions, using the conventional implementation (MI-Full).

```PearsonEWC.m``` - computes the Partial Correlation (PC) between all pairs of regions, using the EWC implementation (PC-EWC).

```TransferEnt.m``` - computes the bivariate Transfer Entropy (TE) between all pairs of regions, using the EWC implementation (TE-EWC).

```TransferEntFull.m``` - computes the bivariate Transfer Entropy (TE) between all pairs of regions, using the conventional implementation (TE-Full).

```matrestruct.m``` - used to restructure $N \times N$ matrices as per some pre-defined mapping.
