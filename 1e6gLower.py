from PBH_MS import CosmologicalModel

# Define model 2
c2 = 0.05
c3 = 0.075
phi0 = 0.01

# Paramaters to optimise xi
h0 = 8.5 # Initial value for inflaton
maxN = 97 # Maximum number of e-folds to calculate up to
xismall = 3250 # Lower bound on xi
xibig = 3260 # Upper bound on xi
maxIteration = 40 # Maximum number of iterations to find optimal xi
lam = 4e-5 # An initial guess for lambda

# For tuning lambda
efoldsDesired = 51.4
lamMin = 2.8e-2
lamMax = 4e-2

# Create an instance of the model
model = CosmologicalModel(c2, c3, phi0, lam, h0, maxN, "1e6g_lower")

# Optimize xi values
model.optimise_xi(xismall, xibig, maxIteration)

# Choose the power
chosenPower = 0.00000000282
power_min = 0.0000000026
power_max = 0.0000000029

model.delete_outputs()
model.tune_power(lamMin,lamMax,efoldsDesired,power_min,power_max,bound="lower")
