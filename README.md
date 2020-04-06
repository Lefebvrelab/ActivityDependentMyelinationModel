# ActivityDependentMyelinationModel

This repository contains code and data for the paper

> Noori et al. (PNAS 2020) *Activity-dependent myelination: a glial mechanism of oscillatory self-organization in large-scale networks*

## Model

The myelin plasticity model consists of a set of stochastic differential equations. This system describing oscillatory neural activity and synchronization behaviour as a function of myelin-related plasticity changes that evolve on the time scale of hours to tens of hours. 

The model is implemented in C++. Figure visualizations are done in Python. 


## Repository organization

The `code` and `data` contain, well, code and data. 

```
code/
data/
  connectivity/
  sim_results
```

Connectivity data (connectivity weights, tract lengths) used by the model are in the `connectivity` folder. Code for simulations and plotting is in the `code` folder. Additionally, the code folder includes python and bash scripts for compiling the C++ code to executable binaries, and submitting as a series of batch job to a compute grid. 


The [figures notebook](https://github.com/Lefebvrelab/ActivityDependentMyelinationModel/blob/master/code/figures.ipynb) loads in the simulation results in the `sim_resuls` folder and produces plots, s shown in the figures of the paper. 

