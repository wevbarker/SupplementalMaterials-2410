from PBH_MS import CosmologicalModel

# Define model 2
c2 = 0.03
c3 = 0.075
phi0 = 0.82

# Paramaters to optimise xi
h0 = 9.5 # Initial value for inflaton
maxN = 75 # Maximum number of e-folds to calculate up to
xismall = 0.24 # Lower bound on xi
xibig = 0.3 # Upper bound on xi
maxIteration = 20 # Maximum number of iterations to find optimal xi
lam = 3e-10 # An initial guess for lambda

# For tuning lambda
efoldsDesired = 51.4
lamMin = 2.4e-10
lamMax = 2.8e-10

# Create an instance of the model
model = CosmologicalModel(c2, c3, phi0, lam, h0, maxN, "1e8g")

# Optimize xi values
model.optimise_xi(xismall, xibig, maxIteration)

# Choose the power
chosenPower = 0.000000215
power_min = 0.0000002
power_max = 0.00000023

model.delete_outputs()
model.tune_power(lamMin,lamMax,efoldsDesired,power_min,power_max)
