from PBH_MS import CosmologicalModel

# Define model 2
c2 = 0.03
c3 = 0.075
phi0 = 0.97

# Paramaters to optimise xi
h0 = 11 # Initial value for inflaton
maxN = 90 # Maximum number of e-folds to calculate up to
xismall = 0.05 # Lower bound on xi
xibig = 0.13 # Upper bound on xi
maxIteration = 24 # Maximum number of iterations to find optimal xi
lam = 1e-10 # An initial guess for lambda

# For tuning lambda
efoldsDesired = 51.4
lamMin = 6.2e-11
lamMax = 7.2e-11

# Create an instance of the model
model = CosmologicalModel(c2, c3, phi0, lam, h0, maxN, "1e9g")

# Optimize xi values
model.optimise_xi(xismall, xibig, maxIteration)

# Choose the power
chosenPower = 0.000001045
power_min = 0.00000099
power_max = 0.0000011

model.delete_outputs()
model.tune_power(lamMin,lamMax,efoldsDesired,power_min,power_max)
