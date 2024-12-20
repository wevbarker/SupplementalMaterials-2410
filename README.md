[![arXiv](https://img.shields.io/badge/arXiv-2410.11948-b31b1b.svg)](https://arxiv.org/abs/2410.11948)

# Supplemental Materials to arXiv:2410.11948

Will Barker, Benjamin Gladwyn and Sebastian Zell 

## How to use this supplement 

The user is assumed to have read [arXiv:2410.11948](https://arxiv.org/abs/2410.11948) in its entirety before using this supplement.

## Software requirements

The scripts in this repository were run using the following software versions:
- _Python_ 3.12.6 of which libraries
    - _matplotlib_
    - _numpy_ 2.0.1
    - _scipy_ 1.14.1
- _Mathematica_ 14.0.0.0
- _xAct_ 1.2.0 of which packages
    - _xPlain_ 0.0.0-developer
Information is available on how to install _Python_ and _xAct_ from the developers' websites. Note that _xPlain_ is available at [this GitHub repository](https://github.com/wevbarker/xPlain). The Python libraries can be installed with 
```console, bash
[user@system SupplementalMaterials-2410]$ pip install matplotlib
```
, and analogously for numpy and scipy.

### System context

The data in this repository were generated on combinations of Cascade Lake, Ice Lake, Sapphire Rapids and Threadripper CPUs, running either Rocky Linux 8 (a rebuild of RHEL8) or Arch Linux.

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

The computationally intensive part of the analysis is the integration of the Mukhanov-Sasaki equation. To accommodate the variable availability of computational resources, the script `PBH_MS.py` uses subprocesses to sample the power spectrum at a variable `M` points in $k$-space, where `M` is introduced at line 28:
```python
        self.M = 600
```
Note that `M` can be adjusted: the pre-set value of `M = 600` is suitable for researchers whose funding supports the basic HPC requirements of precision cosmology (we do not include job-scheduling scripts in this repository, since they are institution-specific). We have confirmed that our final results (i.e. the conclusions we draw from the GW forecast plots) have stabilised by at least `M = 600`, and remain unchanged up to the maximum tested value of `M = 2000`. Consistent with the [_open science_ principle](https://horizoneuropencpportal.eu/sites/default/files/2023-04/task-3.6-open_science_brief.pdf), arbitrarily smaller values of `M` can be used for citizen-science projects on a personal computer (and for development purposes). Unavoidably, the accuracy of the resulting power spectrum and PBH abundance fractions will suffer at lower `M`, and the data products may no longer be science-grade.

### **Step 2:** Stochastic gravitational wave spectra

Once the model scripts have been run, the stochastic GW spectra can be computed. We provide code for the calculation in _MATLAB_, _Wolfram Language_ and _Python_. Performing the integration with _Wolfram Language_ or _Python_ results in numerical instabilities, with Python typically performing worse. We posit that these errors would not be difficult to resolve with further development, allowing for a fully open-source pipeline. 

#### _MATLAB_ Implementation ####

The _MATLAB_ script `MatlabGWCalculation.m` calculates the GW spectra for the 4 models and the bounds. The script `MatlabGWCalculation_FittedSpectra.m` calculates the GW spectra for the fitted power spectra. 

#### _Mathematica_ Implementation ####

The _Mathematica_ notebook `PBH_GW.nb` records the outcome of the GW integration procedure, as produced by the one-line script: 
```mathematica
In[1]:= Get@FileNameJoin[{NotebookDirectory[],"PBH_GW.m"}];
```
The source file `PBH_GW.m` can equally well be run from the command line, using:
```console, bash
[user@system SupplementalMaterials-2410]$ math -run < PBH_GW.m
```
Within `PBH_GW.m`, the global parameters `$Compute=True` and `$TheProcessorCount=100` (lines 8 and 9) reflect the intention to actually _compute_ the spectra in an HPC environment (for which the command line is more appropriate). The output of this process is a binary such as `1e6g_lower_GW.mx`. Note that you must have a Mathematica license which supports the number of parallel processes you are trying to run, otherwise you will encounter errors such as `No valid password found`. As with the main _Python_ model scripts, the integrals can be run on a personal computer with the variable `$TheProcessorCount` set between (e.g.) four and eight. In this case, attention is drawn to line 64:
```mathematica
k1=10^Range[Log10[lowF/(conv)],Log10[highF/(conv)],(Log10[highF/(conv)]-Log10[lowF/(conv)])/1000];
```
The value `1000` may need ammending, so that fewer frequencies are sampled and the computation can take place in a reasonable wallclock time. By setting `$Compute=False`, the script will instead _load_ the pre-computed spectra from binaries such as `1e6g_lower_GW.mx` and plot them (evidently, the notebook is more suitable in this case). Note that the script requires the _xPlain_ package, an open-source contribution to _xAct_ which can be found at [this GitHub repository](https://github.com/wevbarker/xPlain). Note also that the binary files are, in principle, architecture-dependent, and may not be compatible with all systems. For this reason we also include files such as `1e6g_lower_GW.csv`, which contain the same data in a human-readable format, and which are produced from the binaries by running with `$Compute=False`. The remaining files pertain to the GW forecasts from fitted power spectra. We provide only the data associated with these spectra, rather than the fitting scripts that generate them (these can be made available upon request). The GW forecasts are computed by `PBH_GW_fitted.m` and `PBH_GW_fitted.nb`, which have an analogous structure to the GW forecast scripts for the models.

#### _Python_ Implementation ####

The script `PythonGWCalcluation.py` demonstrates the code used to calculate GW signatures in Python. The use of `scipy.integrate.nquad` results in a numerical instabilities at low frequencies.
