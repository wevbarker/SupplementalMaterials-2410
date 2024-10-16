[![arXiv](https://img.shields.io/badge/arXiv-2311.11790-b31b1b.svg)](https://arxiv.org/abs/2410.#####)

# Supplemental Materials to arXiv:2410.#####

Will Barker, Benjamin Gladwyn and Sebastian Zell 

## How to use this supplement 

The user is assumed to have read [arXiv:2410.#####](https://arxiv.org/abs/2410.#####) in its entirety before using this supplement.

### **Step 1:** Mukhanov-Sasaki integration and Press-Schechter formalism

This supplement contains analyses of the four models in Table 1 in Section 4 of the manuscript, corresponding to four different PBH masses and, for each mass, to three variants of the Press-Schechter formalism. Each PBH mass is associated with the string `1e#g`, where `#` ranges from six to nine, and the syntax encodes the PBH mass in grams. In the names of the main _Python_ files, this string may be followed by `Lower`, nothing, or `Upper`, representing a lower bound, representative value and upper bound on stochastic GW production. These three variants are defined by specific choices of the critical threshold $\delta_c$ and window function $W(x)$ used in the Press-Schechter formalism: the details are explained in Section 3.3 of the manuscript, and the implementation can also be verified by examination of `PBH_MS.py`. Each analysis is begun by running the corresponding _Python_ script, which performs the automated parts of the tuning procedure presented in Figure 3 of Section 3 of the manuscript. For example, to run the lower bound on GW production for PBHs of mass $10^6$ grams, we use the following:
```console, bash
[user@system SupplementalMaterials-2410]$ python3 1e6gLower.py
```
All such model scripts call the files `PBH_MS.py` and `radau.py`. The output files for this model are:
- `1e6g_lower_outputs.txt`: the spectral index, tensor-to-scalar ratio and e-folds to the end of inflation
- `1e6g_lower_PBHAbundance.txt`: dark matter fraction of PBHs as a function of PBH mass
- `1e6g_lower_PBHAbundance.pdf`: plot of the above 
- `1e6g_lower_PowerSpectrum.txt`: primordial power spectrum as a function of wavenumber
- `1e6g_lower_PowerSpectrum.pdf`: plot of the above

The computationally intensive part of the analysis is the integration of the Mukhanov-Sasaki equation. To accommodate the variable availability of computational resources, the script `PBH_MS.py` uses subprocesses to sample the power spectrum at a variable `M = 1000` points in $k$-space. Note that `M` can be adjusted: the pre-set value is suitable for researchers whose funding supports the basic HPC requirements of precision cosmology (we do not include job-scheduling scripts in this repository, since they are institution-specific). Consistent with the [_open science_ principle](https://horizoneuropencpportal.eu/sites/default/files/2023-04/task-3.6-open_science_brief.pdf), arbitrarily smaller values of `M` can be used for citizen-science projects on a personal computer (and for development purposes). Unavoidably, the accuracy of the resulting power spectrum and PBH abundance fractions will suffer at lower `M`, and the data products may no longer be science-grade.

### **Step 2:** Stochastic gravitational wave spectra

Once the model scripts have been run, the stochastic GW spectra can be computed. For expedience, the integration is performed using a _Wolfram Language_ script, simply to avoid numerical instabilities which we encountered when running our first-pass _Python_ implementation (which we do not include). We posit that these errors would not be difficult to resolve with further development, allowing for a fully open-source pipeline. The _Mathematica_ notebook `PBH_GW.nb` records the outcome of this integration procedure, as produced by the one-line script: 
```mathematica
In[1]:= Get@FileNameJoin[{NotebookDirectory[],"PBH_GW.m"}];
```
The source file `PBH_GW.m` can equally well be run from the command line, using:
```console, bash
[user@system SupplementalMaterials-2410]$ math -run < PBH_GW.m
```
Within `PBH_GW.m`, the global parameters `$Compute=True` and `$TheProcessorCount=100` reflect the intention to actually _compute_ the spectra in an HPC environment (for which the command line is more appropriate). As with the main _Python_ model scripts, the integrals can be run on a personal computer with the variable `$TheProcessorCount` set between (e.g.) four and eight. In this case, the line
```mathematica
k1=10^Range[Log10[lowF/(conv)],Log10[highF/(conv)],(Log10[highF/(conv)]-Log10[lowF/(conv)])/1000];
```
may need ammending, so that fewer frequencies are sampled and the computation can take place in a reasonable wallclock time. By setting `$Compute=False`, the script will instead _load_ the pre-computed spectra from binaries such as `1e6g_lower_GW.mx` and plot them (evidently, the notebook is more suitable in this case). Note that the script requires the _xPlain_ package, an open-source contribution to _xAct_ which can be found at [this GitHub repository](https://github.com/wevbarker/xPlain). Note also that the binary files are, in principle, architecture-dependent, and may not be compatible with all systems. For this reason we also include files such as `1e6g_lower_GW.csv`, which contain the same data in a human-readable format.

### Miscellaneous fitting scripts

The remaining files in this repository pertain to the GW forecasts from fitted power spectra. We provide only the data associated with these spectra, rather than the fitting scripts that generate them (these can be made available upon request). The GW forecasts are computed by `PBH_GW_fitted.m` and `PBH_GW_fitted.nb`, which have an analogous structure to the GW forecast scripts for the models.

### System context

The data in this repository were generated on combinations of Cascade Lake, Ice Lake, Sapphire Rapids and Threadripper CPUs, running either Rocky Linux 8 (a rebuild of RHEL8) and Arch Linux.
