import matplotlib.pyplot as plt
from scipy import optimize
import numpy as np
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.integrate import solve_ivp, odeint, quad
from scipy.optimize import minimize
# Will: extra imports
import subprocess
import time
import os
import shutil
# Will: end of extra imports

class CosmologicalModel:
    def __init__(self, c2, c3, phi0, lam, h0, maxN, name):
        self.c2 = c2
        self.c3 = c3
        self.phi0 = phi0
        self.lam = lam
        self.h0 = h0
        self.maxN = maxN
        self.xi = None # To be calculated
        self.peakPowerList = []
        self.evolution_result = None
        self.k_list = None
        self.spectrum = None
        self.name = name
        self.M = 6

    # Define the potential V(phi)
    def V_phi(self, phi, xi, lam):
        phi0 = self.phi0
        c2 = self.c2
        c3 = self.c3
        return lam / (24 * (1 + xi * phi**2)**2) * (
        3 * phi**4 + xi**2 * phi0**4 * phi**4 - 8 * (1 + c3) * phi0 * phi**3 +
        2 * (1 + c2) * (3 + xi * phi0**2) * phi0**2 * phi**2)

    # Define the Weyl transformation and canonical normalization
    def hofphi(self,phi, xi):
        return np.sqrt((1 + 6 * xi) / xi) * np.arcsinh(phi * np.sqrt(xi * (1 + 6 * xi))) - \
            np.sqrt(6) * np.arctanh((np.sqrt(6) * xi * phi) / np.sqrt(1 + xi * (1 + 6 * xi) * phi**2))

    # Function to find the root of
    def root_function(self,phi, h, xi):
        return self.hofphi(phi, xi) - h

    # Function to find phi for an arbitrary value of h
    # Will: start of modified lines 
    def find_phi_for_h(self,h, xi):
        for phi_initial_guess in np.linspace(1.0, 12.0, 100):
            solution = optimize.root(lambda phi, h: self.root_function(phi, h, xi), phi_initial_guess, args=(h,))
            if solution.success:
                return solution.x[0]
                break
            else:
                continue
    # Will: end of modified lines

    # Define the potential V(h) in terms of h
    def V_h(self,h, xi, lam):
        phi = self.find_phi_for_h(h, xi)
        return self.V_phi(phi, xi, lam)

    def h_derivatives(self,x, Ne, V_h_interp, V_h_grad_interp):
        return [x[1], -3 * x[1] + x[1]**3 / 2 - (3 - (x[1]**2 / 2)) * V_h_grad_interp(x[0])/V_h_interp(x[0])]

    def calculate_a(self):
        print("The values for a are:")
        print(f'a2 = {self.lam / 24 * 2 * (1 + self.c2) * (3 + self.xi * self.phi0**2) * self.phi0**2}')
        print(f'a3 = {self.lam / 24 * (-8) * (1 + self.c3) * self.phi0}')
        print(f'a4 = {self.lam / 24 * (3 + self.xi**2 * self.phi0**4)}')

    def CalculateDNe(self,xi=None):
        if self.xi == None and xi == None:
            raise ValueError("xi has not been set - please calculate or provide a value")
            return
        
        if xi == None:
            xi = self.xi

        # Generate values for h between 0 and 2
        h_values = np.linspace(0, 14, 1000)
        V_h_values = [self.V_h(h, xi, self.lam) for h in h_values]

        # Create the interpolation function for V(h)
        V_h_interp = UnivariateSpline(h_values, V_h_values, k=4,s=0)

        # Generate values for h between 0 and 2 for plotting the gradient
        V_h_grad = V_h_interp.derivative()
        grad_V_h_values = V_h_grad(h_values)

        # Solve the background evolution
        NeArray = np.arange(0,self.maxN,1e-3)

        H0 = np.sqrt(V_h_interp(self.h0) / 3)
        dh = - V_h_grad(self.h0) / (3 * H0)
         
        # Define an event function to stop when h crosses zero
        def h_crosses_zero_event(Ne, x):
            return x[0]  # This triggers when h crosses zero

        # Event properties
        h_crosses_zero_event.terminal = True  # Stop the solver when the event is triggered
        h_crosses_zero_event.direction = -1   # Only trigger when h is decreasing (crossing zero from positive to negative)

        # Solve the system using solve_ivp with the event
        sol = solve_ivp(lambda Ne, x: self.h_derivatives(x, Ne, V_h_interp, V_h_grad),
            (NeArray[0], NeArray[-1]),
            [self.h0, dh],
            t_eval=NeArray,
            method='DOP853',  # or 'DOP853' for better performance
            atol=1e-12,     # Absolute tolerance
            rtol=1e-12,      # Relative tolerance,
            events=h_crosses_zero_event)

        # Access the solution up to the stopping point
        h_sol = sol.y[0]  # First variable (h_sol)
        dhdNe_sol = sol.y[1]  # Second variable (dhdNe_sol)
        NeArray = sol.t
         
        # Check if the event was triggered 
        #if sol.t_events[0].size > 0:
        #    print(f"Stopped early because h crossed zero at Ne = {sol.t_events[0][0]}")
        
        if sol.t_events[0].size == 0:
            #print("Inflaton stuck")
            return

        # Identify when/if h goes below 0.4 so that it has left the USR phase
        try:
            Ne_h_lt_0dot4 = NeArray[np.where(h_sol < 0.4)[0][0]]
        except:
            print('An error ocurred')

        #Ne_h_lt_0 = NeArray[np.where(h_sol < 0)[0][0]]

        epsilon = dhdNe_sol**2/2
        NeUpper = NeArray[np.where(epsilon > 1)[0][0]]

        # Create a linear interpolation of epsilon
        #epsilon_interp = interp1d(NeArray, epsilon, kind='linear')
        epsilon_interp = UnivariateSpline(NeArray, epsilon, k=1,s=0)

        # Find the root using brentq within the interval [0, NeUpper]
        Nend = optimize.brentq(lambda Ne: epsilon_interp(Ne)-1.0, Ne_h_lt_0dot4, NeArray[-1])

        #print(f"The value of Ne where epsilon > 1 for the first time is: {Nend}")

        # Compute H and P
        H = np.sqrt(2 * V_h_interp(h_sol) / (6.0 - dhdNe_sol**2))
        H_interp = interp1d(NeArray, H, kind='linear')

        P = H**2 / (8 * np.pi**2 * epsilon)
        
        # Interpolate the power spectrum P
        P_interp = interp1d(NeArray, P, kind='linear')
 
        # Find the root using brentq within the interval [0, NeArray[-1]]
        As = 1e-10 * np.exp(3.044)
        NeCMB = optimize.brentq(lambda Ne: P_interp(Ne)-As, 0, Nend)
        h_interp = UnivariateSpline(NeArray, h_sol, k=1,s=0)
        #print(f"h at CMB: {h_interp(NeCMB)}")

        #print(f"The value of Ne where P equals As is: {NeCMB}")

        # Define the wavenumbers calibrated to the pivot scale at 0.05
        k = 0.05 * H * np.exp(NeArray-NeCMB) / H_interp(NeCMB)
        k_interp = interp1d(NeArray, k, kind='linear')

        ### Determine the maximum Power
        # Restrict the range to NeCMB and Nend
        Ne_restricted = NeArray[(NeArray >= NeCMB) & (NeArray <= Nend)]
        P_restricted = P[(NeArray >= NeCMB) & (NeArray <= Nend)]
        k_restricted = k[(NeArray >= NeCMB) & (NeArray <= Nend)]

        # Find the maximum value of P in the restricted range
        P_max = np.max(P_restricted)
        Ne_max = Ne_restricted[np.argmax(P_restricted)]
 
        dNe = Nend - NeCMB
        #print(f'{dNe} efolds between CMB and end of inflation')

        return NeArray, epsilon_interp, h_interp, dhdNe_sol, k_restricted, P_restricted, k_interp, H_interp, NeCMB, dNe, P_max, P_interp

    def optimise_xi(self, xismall, xibig, maxIteration):
        for i in range(0,maxIteration):
            xitemp = (xibig+xismall)/2
            evolution_result = self.CalculateDNe(xitemp)
            if evolution_result == None:
                #print('Inflaton stuck, decreasing xi')
                xibig = xitemp
            else:
                #print('Increasing xi')
                NeArray, epsilon_interp, h_interp, dhdNe_sol, k_restricted, P_restricted, k_interp, H_interp, NeCMB, dNe, P_max, P_interp = evolution_result
                xismall = xitemp
                self.peakPowerList.append((xitemp,P_max))

    #def identify_xi(self, chosenPower,plot=False):
    def identify_xi(self, chosenPower,plot=True):
        # Convert peakPowerList to numpy array for easier manipulation
        peakPowerArray = np.array(self.peakPowerList)

        # Extract xi and peak power values
        xi_values = peakPowerArray[:, 0]
        peak_power_values = peakPowerArray[:, 1]

        # Create logarithmic interpolation of the peak power
        log_xi_values = np.log10(xi_values)
        log_peak_power_values = np.log10(peak_power_values)
        log_peak_power_interp = interp1d(log_xi_values, log_peak_power_values, kind='linear')

        # Define the peak power function for a given xi
        def peakPower(xi1):
            return 10 ** log_peak_power_interp(np.log10(xi1))

        xi_plot = np.linspace(np.min(xi_values), np.max(xi_values), 1000)

        # Find the root xi1 where the peak power equals the chosen power
        xi_min = np.min(xi_values)
        xi_max = np.max(xi_values)
        
        try:
            result = optimize.brentq(lambda xi1: peakPower(xi1) - chosenPower, xi_min, xi_max)
            self.xi = result

            if plot==True:
                plt.figure(figsize=(10, 6))
                plt.axvline(self.xi, color='g', linestyle='--', label=f'Optimal xi = {self.xi:.5f}')
                plt.loglog(xi_values, peak_power_values, 'o', label='Peak Power Data')
                plt.loglog(xi_plot, peakPower(xi_plot), '-', label='Interpolated Curve')
                plt.axhline(chosenPower, color='r', linestyle='--', label=f'Chosen Power = {chosenPower}') 
                plt.xlabel('xi')
                plt.ylabel('Peak in Power Spectrum')
                plt.title('Peak in Power Spectrum vs xi')
                plt.legend()
                plt.grid(False)
                #plt.show()
            print(f"The optimal value of xi is: {self.xi}")

        except:
            print("Value not found")

    def optimise_lam_visual(self,lamMin,lamMax,efoldsDesired,step):
        def objectiveFunction(lam):
            self.lam = lam
            evolution_result = self.CalculateDNe(self.xi)
            if evolution_result is None:
                return np.inf  # If the result is None, return a large value
            NeArray, epsilon_interp, h_interp, dhdNe_sol, k_restricted, P_restricted, k_interp, H_interp, NeCMB, dNe, P_max, P_interp = evolution_result
            return (dNe - efoldsDesired) ** 2

        # Sample the objective function at values in the interval
        q = np.arange(lamMin, lamMax, step)
        values = [objectiveFunction(elem) for elem in q]

        # Plot the objective function
        plt.figure(figsize=(10, 6))
        plt.plot(q, values, label='Objective Function', marker='o')
        plt.xlabel('Lambda')
        plt.ylabel('Objective Function')
        plt.title('Objective Function vs Lambda')
        plt.legend()
        plt.grid(True)
        #plt.show()

        # Find the lambda that minimizes the objective function
        minIndex = np.argmin(values)
        optimalLambda = q[minIndex]
        minValue = values[minIndex]

        print(f"The optimal value of lambda is: {optimalLambda}")
        print(f"The minimum value of the objective function is: {minValue}")

        self.lam = optimalLambda

    def optimise_lam_procedure(self,lamMin,lamMax,efoldsDesired):
        def objectiveFunction(lam):
            self.lam = lam
            evolution_result = self.CalculateDNe(self.xi)
            if evolution_result is None:
                return np.inf  # If the result is None, return a large value
            NeArray, epsilon_interp, h_interp, dhdNe_sol, k_restricted, P_restricted, k_interp, H_interp, NeCMB, dNe, P_max, P_interp = evolution_result
            return (dNe - efoldsDesired) ** 2

        # Define the bounds for lambda
        bounds = [(lamMin, lamMax)]

        # Perform the minimization using the 'Nelder-Mead' method as it's robust for non-linear optimization
        minimizationResult = minimize(objectiveFunction, x0=(lamMin + lamMax) / 2, bounds=bounds, method='Nelder-Mead', options={'maxiter': 100, 'disp': False})

        # Extract the optimal lambda value from the minimization result
        self.lam = minimizationResult.x[0]

        print(f"The optimal value of lambda is: {self.lam}")

    # Will: delete the output file if it exists
    def delete_outputs(self):
        os.remove(self.name+'_outputs.txt') if os.path.exists(self.name+'_outputs.txt') else None

    def calculate_outputs(self,savefile=True):
        evolution_result = self.CalculateDNe()
        self.evolution_result = evolution_result

        NeArray, epsilon_interp, h_interp, dhdNe_sol, k_restricted, P_restricted, k_interp, H_interp, NeCMB, dNe, P_max, P_interp = evolution_result
        Nend = NeCMB+dNe

        # Restrict NeArray between NeCMB and Nend
        mask = (NeArray >= NeCMB) & (NeArray <= Nend)
        NeArray_restricted = NeArray[mask]

        # Calculate eta for restricted NeArray
        eta = epsilon_interp(NeArray_restricted) - epsilon_interp.derivative()(NeArray_restricted) / (2 * epsilon_interp(NeArray_restricted))
        eta_interp = UnivariateSpline(NeArray_restricted, eta, k=1, s=0)

        ns = 1 + 2 * eta_interp(NeCMB) - 6 * epsilon_interp(NeCMB)
        r = 16 * epsilon_interp(NeCMB)

        print(f'ns = {ns} , Measured: 0.9649+-0.0042')
        print(f'r = {r} , r<0.056')
        print(f'e-folds = {dNe}')

        if savefile == True:
            # Open a file with the name stored in self.name
            if not os.path.isfile(f'{self.name}_outputs.txt'):
                with open(f'{self.name}_outputs.txt', 'w') as file:
                    # Write the formatted outputs to the file
                    file.write(f'ns = {ns} , Measured: 0.9649+-0.0042\n')
                    file.write(f'r = {r} , r<0.056\n')
                    file.write(f'e-folds = {dNe}\n')
            else:
                with open(f'{self.name}_outputs.txt', 'a') as file:
                    # Write the formatted outputs to the file
                    file.write(f'ns = {ns} , Measured: 0.9649+-0.0042\n')
                    file.write(f'r = {r} , r<0.056\n')
                    file.write(f'e-folds = {dNe}\n')

    def calculate_MS_spectrum(self, plot=False, savefile=False):
        NeArray, epsilon_interp, h_interp, dhdNe_sol, k_restricted, P_restricted, k_interp, H_interp, NeCMB, dNe, P_max, P_interp = self.evolution_result
        
        # Calculate eta for restricted NeArray
        Nend = NeCMB+dNe
        mask = (NeArray >= NeCMB) & (NeArray <= Nend)
        NeArray_restricted = NeArray[mask]
        eta = epsilon_interp(NeArray_restricted) - epsilon_interp.derivative()(NeArray_restricted) / (2 * epsilon_interp(NeArray_restricted))
        eta_interp = UnivariateSpline(NeArray_restricted, eta, k=1, s=0)

        # Constants and initial conditions
        Mpc = 2.4e18  # GeV
        H0_Mpc = 1.56e38  # in units Mpc^-1

        # Function to find Ne for an arbitrary value of k
        def find_Ne_for_k(k, k_interp):
            solution = optimize.brentq(
                lambda Ne: k_interp(Ne) - k, 
                NeCMB - 5, Nend
            )
            return solution

        def Func(Ne, k):
            return (k/k_interp(Ne))**2 + (1 + epsilon_interp(Ne)-eta_interp(Ne))*(eta_interp(Ne)-2) - (epsilon_interp.derivative()(Ne)-eta_interp.derivative()(Ne))

        a = 0.05 * np.exp(NeArray - NeCMB) / H_interp(NeCMB)
        z = a * dhdNe_sol
        z_interp = interp1d(NeArray, z, kind='linear')

        # Mukhanov Sasaki Equation
        #M = 600  # Points on graph

        # Starting value of k   k_interp(NeCMB)
        k_list = np.logspace(np.log10(k_interp(NeCMB)), np.log10(k_interp(Nend)), self.M)

        # Will: extra lines here 
        self.data_dir = self.name+'WideData'
        shutil.rmtree(self.data_dir, ignore_errors=True)
        os.mkdir(self.data_dir)
        # Will: end of extra lines

        # Evolve from
        N_start = find_Ne_for_k(k_interp(NeCMB)/20, k_interp)
        N = np.linspace(N_start, Nend, 100000)
        g = (1-epsilon_interp(N))/2

        # 10000 up to 1e15

        # Will: extra lines here
        def findnewsol(self,i,N,w,g,Nin,Nout,Ruin,dRuin):
            if not os.path.isfile(self.data_dir+'/N.npy'):
                np.save(self.data_dir+'/N',N)
                np.save(self.data_dir+'/g',g)
            np.save(self.data_dir+'/lims_'+str(i),(Nin,Nout,Ruin,dRuin))
            np.save(self.data_dir+'/w_'+str(i),w)
        # Will: end of extra lines

        # Will: extra lines here
        for i, k in enumerate(k_list):
            Nin = find_Ne_for_k(k/20, k_interp)
            Nout = Nin + 11 if (Nin + 11 < Nend) else Nend
            scaling = 1e4
            Ruin = 1 / scaling
            dRuin = 0
            Iuin = 0
            dIuin = -20 / scaling
            w = np.sqrt(Func(N, k)+0j)
            findnewsol(self,i,N,w,g,Nin,Nout,Ruin,dIuin*1j)
        # Will: end of extra lines

        # Will: extra lines here
        process_list = []
        for i, k in enumerate(k_list):
            command_string='python3 radau.py '+str(i)+' "'+self.data_dir+'"'
            #print(command_string)
            process_list.append(subprocess.Popen(command_string,shell=True))

        for process in process_list:
            process.wait()
        # Will: end of extra lines

        # Will: extra lines here
        old_sol = None
        spectrum = np.zeros_like(k_list)
        self.data_dir = self.name+'WideData'
        for i, k in enumerate(k_list):
            Nin = find_Ne_for_k(k/20, k_interp)
            Nout = Nin + 11 if (Nin + 11 < Nend) else Nend
            scaling = 1e4
            try:
                sol_=np.load(self.data_dir+'/p_'+str(i)+'.npy')
                old_sol=sol_
            except:
                sol_=old_sol
            MSPower = (scaling**2 / (2 * k)) * (k**3 / (2 * np.pi**2)) * (np.abs(sol_)**2 / (z_interp(Nout))**2)
            spectrum[i] = MSPower
        # Will: end of extra lines

        self.k_list = k_list
        self.spectrum = spectrum    

        # Will: extra lines here
        plot = True
        savefile = True
        # Will: end of extra lines

        # Plot the power spectrum
        if plot == True:
            plt.figure()
            plt.loglog(k_list, spectrum, label="Mukhanov-Sasaki")
            plt.plot(k_interp(NeArray_restricted), P_interp(NeArray_restricted), label='slow-roll')
            plt.xlabel(r'$\mathrm{Wavenumber}\ k\ \mathrm{[Mpc}^{-1}\mathrm{]}$', fontsize=12)
            plt.ylabel(r'$\mathrm{Comoving\ Curvature\ Spectrum}\ \mathcal{P}_\mathcal{R}(k)$', fontsize=12)
            plt.title(r'Comoving Curvature Spectrum', fontsize=14)
            plt.tight_layout()  # Adjust layout to prevent overlap
            plt.legend(loc='upper left')
            plt.xlim(1e-6, 1e21)
            plt.ylim(1e-15, 1)
            #plt.show()
        
        if savefile == True:
            # Plotting the figure
            plt.figure()
            plt.loglog(k_list, spectrum, label="Mukhanov-Sasaki")
            plt.plot(k_interp(NeArray_restricted), P_interp(NeArray_restricted), label='slow-roll')
            plt.xlabel(r'$\mathrm{Wavenumber}\ k\ \mathrm{[Mpc}^{-1}\mathrm{]}$', fontsize=12)
            plt.ylabel(r'$\mathrm{Comoving\ Curvature\ Spectrum}\ \mathcal{P}_\mathcal{R}(k)$', fontsize=12)
            plt.title(r'Comoving Curvature Spectrum', fontsize=14)
            plt.tight_layout()  # Adjust layout to prevent overlap
            plt.legend(loc='upper left')
            plt.xlim(1e-6, 1e21)
            plt.ylim(1e-15, 1)

            # Saving the figure
            plt.savefig(f'{self.name}_PowerSpectrum.pdf', bbox_inches='tight')

            # Saving the k_list and spectrum data to a txt file without a header
            data = np.column_stack((k_list, spectrum))
            np.savetxt(f'{self.name}_PowerSpectrum.txt', data, delimiter=' ')
            


    def calculate_PBH_abundance(self, lower_range = 1e-6, upper_range = 1, bound=None,plot=False,savefile=False):
        # Interpolate the power spectrum
        logPInterpolation = interp1d(np.log10(self.k_list), np.log10(self.spectrum), kind='linear', fill_value='extrapolate')

        def P(k):
            return 10 ** logPInterpolation(np.log10(k))

        # Constants and functions
        def kpbh(PBHMass):
            return 2 * 10**14 * np.sqrt(10**17 / PBHMass)

        def mpbh(k):
            return 10**17 / (k / 2e14)**2

        if bound == None:
            def W(p):
                return np.exp(-p**2 / 2)  # Window function
            deltac = 0.45

        elif bound == "upper":
            def W(p):
                return np.exp(-p**2 / 2)
            deltac = 0.6

        elif bound == "lower":
            def W(p):
                return 3*(np.sin(p)-p * np.cos(p))/p**3  # Window function
            deltac = 0.4

        # Define the sigmaSquared function
        def sigmaSquared(PBHMass):
            # Define the integrand function
            kpbh_val = kpbh(PBHMass)
            def integrand(q):
                return P(q) * (q**3 / kpbh_val**4) * W(q / kpbh_val)**2

            # Perform the numerical integration from 0 to infinity
            result, error = quad(integrand, 0, 10*kpbh_val,limit=100,epsrel = 1.49e-6)
            return (16 / 81) * result


        def fPBH(PBHMass):
            sig2 = sigmaSquared(PBHMass)

            def integrand(delta):
                return np.exp(-delta**2 / (2 * sig2))

            prefactor = np.sqrt(10**18 / PBHMass) * (1.25 * 10**15 / np.sqrt(2 * np.pi * sig2))
            result, _ = quad(integrand, deltac, np.inf)
            return prefactor * result

        lowM = 1e5
        highM = 1e10
        M1 = np.logspace(np.log10(lowM), np.log10(highM), 300)
        fPBHValues = np.array([fPBH(m) for m in M1])

        if savefile == True:
            plt.figure()
            plt.loglog(M1, fPBHValues)
            plt.xlabel('Mass / g')
            plt.ylabel('f_PBH')
            plt.ylim(lower_range,upper_range)
            plt.title('Log-Log Plot of f_PBH vs Mass')
            plt.savefig(f'{self.name}_PBHAbundance.pdf', bbox_inches='tight')

            # Saving the k_list and spectrum data to a txt file without a header
            data = np.column_stack((M1, fPBHValues))
            np.savetxt(f'{self.name}_PBHAbundance.txt', data, delimiter=' ')

        elif plot==True:
            plt.loglog(M1, fPBHValues)
            plt.xlabel('Mass / g')
            plt.ylabel('f_PBH')
            plt.ylim(lower_range,upper_range)
            plt.title('Log-Log Plot of f_PBH vs Mass')
            #plt.show()
        else:
            return np.max(fPBHValues)

    def tune_power(self, lamMin,lamMax, efoldsDesired, power_min, power_max,bound=None):

        def objectiveFunction(chosenPower):
            self.identify_xi(chosenPower)
            self.optimise_lam_procedure(lamMin,lamMax,efoldsDesired)
            self.calculate_outputs()
            self.calculate_MS_spectrum()
            PBH_abundance = self.calculate_PBH_abundance(bound=bound)
            return np.log(PBH_abundance)**2

        bounds = [(power_min, power_max)]
        minimizationResult = minimize(objectiveFunction, x0=(power_min+power_max) / 2, bounds=bounds, method='Nelder-Mead', options={'maxiter': 10, 'disp': False})

        # Extract the optimal lambda value from the minimization result
        optimum_power =  minimizationResult.x[0]

        # Final result
        self.identify_xi(optimum_power)
        self.optimise_lam_procedure(lamMin,lamMax,efoldsDesired)
        self.calculate_outputs(savefile=True)
        self.calculate_MS_spectrum(plot=False, savefile=True)
        self.calculate_PBH_abundance(bound=bound,plot=False,savefile=True) 
