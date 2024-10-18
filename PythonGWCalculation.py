import numpy as np
from scipy.integrate import nquad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import os
print(os.getcwd())


# Open the power spectrum file
data = np.genfromtxt('1e8g_PowerSpectrum.txt', delimiter = ' ')
k_list = data[:,0]
spectrum = data[:,1]

# Interpolate the power spectrum
logPInterpolation = interp1d(np.log10(k_list), np.log10(spectrum), kind='linear', fill_value='extrapolate')

def P(k):
    return 10 ** logPInterpolation(np.log10(k))

# Constants
cg = 0.4
OmegaR0 = 2.47e-5

# Functions
def Ic(d, s):
    return -36 * np.pi * (((s**2 + d**2 - 2)**2 / (s**2 - d**2)**3) * (np.heaviside(s-1,1)))

def Is(d, s):
    return -36 * ((s**2 + d**2 - 2) / (s**2 - d**2)**2) * \
           ((s**2 + d**2 - 2) / (s**2 - d**2) * np.log((1 - d**2) / abs(s**2 - 1)) + 2)

# Integrand for GW signature
def Integrand2(k, d, s):
    term1 = ((d**2 - 1/3) * (s**2 - 1/3) / (s**2 - d**2))**2
    term2 = P(k * np.sqrt(3) * (s + d) / 2)
    term3 = P(k * np.sqrt(3) * (s - d) / 2)
    term4 = Is(d, s)**2 + Ic(d, s)**2
    return term1 * term2 * term3 * term4

# GW signature calculation
# Define h2OmegaGW2 as a function of k
def h2OmegaGW2(k):
    integrand = lambda d, s: Integrand2(k, d, s)
    result, error = nquad(integrand, [(0, 1/np.sqrt(3)), (1/np.sqrt(3), np.inf)], opts={'limit': 100,'epsabs':1e-22})
    print(f"Calculated wavenumber: {k}, result is h2Omega = {result}, error = {error}")
    return cg * OmegaR0 / 36 * result

# Works from 1e-1  opts={'limit': 50, 'epsrel': 1e-6,'epsabs':1e-20}
# Works to 7e-2 'epsrel':1e-22,'epsabs':1e-23
# Works to 7e-2 'epsrel':1e-22,'epsabs':1e-24

# Frequency region to plot GW signature in
lowF = 1e-2
highF = 1e6
conv = 1.546e-15
k1 = np.logspace(np.log10(lowF / conv), np.log10(highF / conv), 100)

# Calculate h2OmegaGW2 for each k1
h2OmegaGW2Sol = [h2OmegaGW2(k) for k in k1]

# Plotting
plt.figure()
plt.loglog(conv * k1, h2OmegaGW2Sol, '-')
plt.xlabel('f/Hz')
plt.ylabel('h2OmegaGW2(k)')
plt.title('Log-Log Plot of OmegaGW2(k)')
plt.grid(True)
plt.show()
