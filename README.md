# Supplemental Materials to arXiv:2410.#####

W. Barker, B. Gladwyn and S. Zell 

## How to use this supplement 

The reader is assumed to have read arXiv:2410.##### in its entirety before using this supplement.

### Mukhanov-Sasaki integration and Press-Schechter formalism

This supplement contains for analyses for Model 2 in Table 1 in Section 4 of the manuscript, applied to four different PBH masses and, for each mass, to three variants of the Press-Schechter formalism. Each PBH mass is associated with the string `1e#g`, where `#` ranges from six to nine, and the syntax encodes the PBH mass in grams. This string may be followed by `_lower_`, `_center_` or `_upper_`, representing a lower bound, representative value and upper bound on stochastic GW production. These three variants are defined by specific choice of the critical threshold $\delta_c$ and window function $W(x)$ used in the Press-Schechter formalism: the details are explained in Section 3.3 of the manuscript, and the implementation can also be found in `PBH_MS.py`. Each analysis is begun by running the corresponing _Python_ script. For example, to run the lower bound on GW production for PBHs of mass $10^6$ grams, we use the following:
```console, bash
[user@system SupplementalMaterials-2410]$ python3 1e6g_lower.py
```
All such model scripts call the files `PBH_MS.py` and `radau.py`. The output files for this model are:
- `1e6g_outputs.txt`: the spectral index, tensor-to-scalar ratio and e-folds to the end of inflation
- `1e6g_PBHAbundance.txt`: dark matter fraction of PBHs as a function of PBH mass
- `1e6g_PBHAbundance.pdf`: dark matter fraction of PBHs as a function of PBH mass
- `1e6g_PowerSpectrum.txt`: primordial power spectrum as a function of wavenumber
- `1e6g_PowerSpectrum.pdf`: primordial power spectrum as a function of wavenumber

The computationally intensive part of the analysis is the integration of the Mukhanov-Sasaki equation. To accommodate the variable availability of computational resources, the script `PBH_MS.py` uses subprocesses to sample the power spectrum at a variable `M = 1000` points in $k$-space. Note that `M` can be adjusted: the pre-set value is suitable for researchers whose funding supports the basic HPC requirements of precision cosmology (we do not include job-scheduling scripts in this repository, since they are institution-specific). Arbitrarily smaller values of `M` can be used for citizen-science projects on a personal computer (and for development purposes), however the accuracy of the resulting power spectrum and PBH abundance fractions will suffer as a result.

### Stochastic gravitational wave spectra

Once the model scripts have been run, the stochastic GW spectra can be computed. The integration is performed using _Mathematica_, to avoid errors encountered when using _Python_: we posit that these errors would not be difficult to resolve with further development, allowing for a fully open-source pipeline. The _Mathematica_ notebook `PBH_GW.nb` records the outcome of this integration procedure, as produced by the one-line script: 
```mathematica
Get@FileNameJoin[{NotebookDirectory[],"PBH_GW.m"}];
```
The source file `PBH_GW.m` can equally well be run from the command line, using:
```console, bash
[user@system SupplementalMaterials-2410]$ math -run < PBH_GW.m
```
Within `PBH_GW.m`, the global parameters `$Compute=True` and `$TheProcessorCount=100` reflects the intention to actually compute the spectra in an HPC environment (for which the command line is more appropriate). By setting `$Compute=False`, the script will instead load the pre-computed spectra from binaries such as `1e6g_lower_GW.mx` and plot them (evidently, the notebook is more suitable in this case). Note that the script requires the _xPlain_ package 
