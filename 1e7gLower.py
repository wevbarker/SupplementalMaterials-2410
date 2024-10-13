from PBH_MS import CosmologicalModel

# Define model 2
c2 = 0.03
c3 = 0.075
phi0 = 0.59

# Paramaters to optimise xi
h0 = 9 # Initial value for inflaton
maxN = 89 # Maximum number of e-folds to calculate up to
xismall = 0.5 # Lower bound on xi
xibig = 1 # Upper bound on xi
maxIteration = 20 # Maximum number of iterations to find optimal xi
lam = 6e-10 # An initial guess for lambda

# For tuning lambda
efoldsDesired = 51.4
lamMin = 1.68e-9
lamMax = 1.76e-9

# Create an instance of the model
model = CosmologicalModel(c2, c3, phi0, lam, h0, maxN, "1e7g_lower")

# Optimize xi values
model.optimise_xi(xismall, xibig, maxIteration)

# Choose the power
chosenPower = 0.00000001515
power_min = 0.000000013
power_max = 0.000000016

model.delete_outputs()
model.tune_power(lamMin,lamMax,efoldsDesired,power_min,power_max,bound="lower")