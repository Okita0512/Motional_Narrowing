from scipy import integrate
import numpy as np

# ==============================================================================================
#                                       Global Parameters     
# ==============================================================================================
conv = 27.211397                            # 1 a.u. = 27.211397 eV
fs_to_au = 41.341                           # 1 fs = 41.341 a.u.
cm_to_au = 4.556335e-06                     # 1 cm^-1 = 4.556335e-06 a.u.
au_to_K = 3.1577464e+05                     # 1 au = 3.1577464e+05 K
kcal_to_au = 1.5936e-03                     # 1 kcal/mol = 1.5936e-3 a.u.

beta = 1052.6

def coth(x):                                # mathematical function, cot(x)
    return 1 / np.tanh(x)

def Bose(x):
    return 1 / (np.exp(beta * x) - 1)

def Drude(x):                               # the molecular bath spectral density function, J_v(w)
    lam = 0.03 / conv
    gam = 0.0248 / conv
    return (2 * lam * gam * x / (x**2 + gam**2)) * coth(beta * x / 2)

"""
sigma^2 = (1 / pi) \int_0^{\infty} dw J_v (w) coth(beta w / 2)
"""

# to get the variance
w0 = 2.0 / conv
Rij = 1.0 # 1.0 / np.sqrt(2 * w0)
wi = np.linspace(1e-15, 0.0248 / conv, 10000000)     # for intergration. Better to be larger (at least 10^3)
y = Drude(wi)
sigma_2 = integrate.trapz(y, wi)
sigma_2 = Rij**2 * sigma_2 / (np.pi)
print("sigma_2 = ", np.sqrt(sigma_2) * conv, '\t eV')

