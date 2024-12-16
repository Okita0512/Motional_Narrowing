import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator, tick_params
fig, ax = plt.subplots()
plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams["font.family"] = "Helvetica"

def Jomg(w):
    lam = 30 / 1000
    gam = 24.8 / 1000
    return 2 * lam * gam * w / (w**2 + gam**2)

def Bose(w):
    beta = 1052.6 / (27.2114)
    return 1. / (np.exp(beta * w) - 1)

Gammac = 8.83 * 5
Nmol = 2

# ====== theoretical linewidth ======
Gammac = Gammac / 1000      # convert meV to eV
Gammae = 0.0765
wc = np.linspace(1.0, 3.0, 50)
w0 = 2.0 + 0.03
Delta = wc - w0
gc = 0.0681
Rabi = gc * np.sqrt(Nmol)

# Hopfield coefficients
C_square = 0.5 * (1 + Delta / np.sqrt(Delta**2 + 4 * Rabi**2))
X_square = 0.5 * (1 - Delta / np.sqrt(Delta**2 + 4 * Rabi**2))

# light-matter mixing angle, defined on [0, \pi / 2]
theta = 0.5 * np.arctan(2 * Rabi / Delta)
for i in range(len(theta)):
    if(theta[i] < 0.0):
        theta[i] = (np.pi / 2) - np.abs(theta[i])

deltaE_UP_D = 0.5 * np.sqrt(Delta**2 + 4 * Rabi**2) + 0.5 * Delta
Gamma_UP_D = ((Nmol - 1) / Nmol) * (1 + np.cos(2 * theta)) * Jomg(deltaE_UP_D) * (1 + Bose(deltaE_UP_D))
Gamma_UP_D_Hopfield = (2 * (Nmol - 1) / Nmol) * C_square * Jomg(deltaE_UP_D) * (1 + Bose(deltaE_UP_D))

deltaE_UP_LP = np.sqrt(Delta**2 + 4 * Rabi**2)
Gamma_UP_LP = (1 / (2 * Nmol)) * (np.sin(2 * theta))**2 * Jomg(deltaE_UP_LP) * (1 + Bose(deltaE_UP_LP))
Gamma_UP_LP_Hopfield = (2 / Nmol) * C_square * X_square * Jomg(deltaE_UP_LP) * (1 + Bose(deltaE_UP_LP))

plt.plot(wc, Gamma_UP_D, "-", color = 'red', label = 'UP -> D: Mixing Angle')
plt.plot(wc, Gamma_UP_D_Hopfield, "o", color = 'red', label = 'UP -> D: Hopfield')

plt.plot(wc, Gamma_UP_LP, "-", color = 'blue', label = 'UP -> LP: Mixing Angle')
plt.plot(wc, Gamma_UP_LP_Hopfield, "o", color = 'blue', label = 'UP -> LP: Hopfield')

plt.legend(frameon = False)
plt.show()