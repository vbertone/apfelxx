import apfelpy as ap
import numpy as np

# Test scale
Mu = 2.

# Reference value of the strong coupling and heavy-quark
# thresholds.
AlphaQCDRef = 0.118
MuQCDRef    = 91.1876
QuarkThresholds = [0, 0, 0, np.sqrt(2), 4.5, 175]

# Iniatialize the running of the coupling at all available
# perturbative orders.
asLO    = ap.AlphaQCD(AlphaQCDRef, MuQCDRef, QuarkThresholds, 0)
asNLO   = ap.AlphaQCD(AlphaQCDRef, MuQCDRef, QuarkThresholds, 1)
asNNLO  = ap.AlphaQCD(AlphaQCDRef, MuQCDRef, QuarkThresholds, 2)
asNNNLO = ap.AlphaQCD(AlphaQCDRef, MuQCDRef, QuarkThresholds, 3)

# Compute and print values at Mu.
print("\nNumerical evolution of the strong coupling:")
print("LO:    alpha_s(Mu = ", Mu, " GeV) = ", asLO.Evaluate(Mu))
print("NLO:   alpha_s(Mu = ", Mu, " GeV) = ", asNLO.Evaluate(Mu) , " (NLO/LO     = ", 100 * asNLO.Evaluate(Mu) / asLO.Evaluate(Mu), "%)")
print("NNLO:  alpha_s(Mu = ", Mu, " GeV) = ", asNNLO.Evaluate(Mu), " (NNLO/NLO   = ", 100 * asNNLO.Evaluate(Mu) / asNLO.Evaluate(Mu), "%)")
print("NNNLO: alpha_s(Mu = ", Mu, " GeV) = ", asNNNLO.Evaluate(Mu), " (NNNLO/NNLO = ", 100 * asNNNLO.Evaluate(Mu) / asNNLO.Evaluate(Mu), "%)")

asLOg    = ap.AlphaQCDg(AlphaQCDRef, MuQCDRef, QuarkThresholds, 0)
asNLOg   = ap.AlphaQCDg(AlphaQCDRef, MuQCDRef, QuarkThresholds, 1)
asNNLOg  = ap.AlphaQCDg(AlphaQCDRef, MuQCDRef, QuarkThresholds, 2)
asNNNLOg = ap.AlphaQCDg(AlphaQCDRef, MuQCDRef, QuarkThresholds, 3)

# Compute and print values at Mu.
print("\nAnalytic evolution of the strong coupling:")
print("LO:    alpha_s(Mu = ", Mu, " GeV) = ", asLOg.Evaluate(Mu))
print("NLO:   alpha_s(Mu = ", Mu, " GeV) = ", asNLOg.Evaluate(Mu) , " (NLO/LO     = ", 100 * asNLOg.Evaluate(Mu) / asLOg.Evaluate(Mu), "%)")
print("NNLO:  alpha_s(Mu = ", Mu, " GeV) = ", asNNLOg.Evaluate(Mu), " (NNLO/NLO   = ", 100 * asNNLOg.Evaluate(Mu) / asNLOg.Evaluate(Mu), "%)")
print("NNNLO: alpha_s(Mu = ", Mu, " GeV) = ", asNNNLOg.Evaluate(Mu), " (NNNLO/NNLO = ", 100 * asNNNLOg.Evaluate(Mu) / asNNLOg.Evaluate(Mu), "%)")

# Reference value of the QED coupling and heavy-quark
# thresholds.
AlphaQEDRef = 1. / 128.
MuQEDRef    = 91.2
LeptThresholds = [0, 0, 1.777]

# Iniatialize the running of the QED coupling at all available
# perturbative orders.
aLO  = ap.AlphaQED(AlphaQEDRef, MuQEDRef, QuarkThresholds, LeptThresholds, 0)
aNLO = ap.AlphaQED(AlphaQEDRef, MuQEDRef, QuarkThresholds, LeptThresholds, 1)

# Compute and print values at Mu.
Mu = 1e10
print("\nNumeric evolution of the electromagnetic coupling:")
print("LO:    alpha_em(Mu = ", Mu, " GeV) = ", aLO.Evaluate(Mu))
print("NLO:   alpha_em(Mu = ", Mu, " GeV) = ", aNLO.Evaluate(Mu) , " (NLO/LO = ", 100 * aNLO.Evaluate(Mu) / aLO.Evaluate(Mu), "%)\n")
