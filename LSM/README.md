# Communication over a network motif (with Linear Stochastic Model (LSM) dynamics)
### File overview
```PC_EWC_noisevar.m, cMI_EWC_noisevar.m, TE_Full_noisevar.m``` - simulates a 4-node network with linear stochastic model (LSM) dynamics, and 3 Poisson-process driven sources (one node isolated), and varying noise amplitudes. The communication between nodes is then inferred using partial correlation (PC-EWC), conditional Mutual Information (cMI-EWC), and Transfer Entropy (TE-Full) respectively.

```PC_EWC_timedelayvar.m, TE_Full_timedelayvar.m``` - similar to ```_noisevar.m```, but with a 3-node motif, 2 Poisson sources, fixed noise amplitude and varying timedelays between nodes 2 and 3. Communication inferred as PC-EWC and TE-Full respectively.

```PC_EWC_firingratevar.m, TE_Full_firingratevar.m``` - similar to ```_timedelayvar.m```, but with fixed noise amplitude and timedelays, and varying source firing rates. Communication inferred as PC-EWC and TE-Full respectively.
