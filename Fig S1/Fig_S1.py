import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator, tick_params
fig, ax = plt.subplots()
plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams["font.family"] = "Helvetica"

conv = 27.211397                            # 1 a.u. = 27.211397 eV
fs_to_au = 41.341                           # 1 fs = 41.341 a.u.
cm_to_au = 4.556335e-06                     # 1 cm^-1 = 4.556335e-06 a.u.
au_to_K = 3.1577464e+05                     # 1 au = 3.1577464e+05 K
kcal_to_au = 1.5936e-03                     # 1 kcal/mol = 1.5936e-3 a.u.

def Jomg(w):
    lam = 30 / 1000
    gam = 24.8 / 1000
    return 2 * lam * gam * w / (w**2 + gam**2)

def Bose(w):
    beta = 1052.6 / (27.2114)
    return 1. / (np.exp(beta * w) - 1)

Gammac = 8.83 * 5
Gammac = Gammac / 1000      # convert meV to eV
Gammae = 0.0765
wc = np.linspace(1.0, 3.0, 50)
w0 = 2.0 + 0.03
Delta = wc - w0
gc = 0.0681

# ==============================================================================================
#                                       Global Parameters     
# ==============================================================================================

lw = 3.0
legendsize = 48         # size for legend
font_legend = {'family':'Times New Roman', 'weight': 'roman', 'size': 22}
# axis label size
lsize = 30             
txtsize = 12
# tick length
lmajortick = 15
lminortick = 5
transparency = .4
y1, y2 = - 0.5, 10.0

unitlen = 7
fig = plt.figure(figsize=(4.0 * unitlen, 0.8 * unitlen), dpi = 128)
plt.subplots_adjust(wspace = 0.2)

# ==============================================================================================
#                           Fig s1a 
# ==============================================================================================
plt.subplot(1, 4, 1)

Nmol = 1
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

deltaE_LP_D = 0.5 * np.sqrt(Delta**2 + 4 * Rabi**2) + 0.5 * Delta
Gamma_LP_D = ((Nmol - 1) / Nmol) * (1 + np.cos(2 * theta)) * Jomg(deltaE_UP_D) * (Bose(deltaE_UP_D))
Gamma_LP_D_Hopfield = (2 * (Nmol - 1) / Nmol) * C_square * Jomg(deltaE_UP_D) * (Bose(deltaE_UP_D))

deltaE_LP_UP = np.sqrt(Delta**2 + 4 * Rabi**2)
Gamma_LP_UP = (1 / (2 * Nmol)) * (np.sin(2 * theta))**2 * Jomg(deltaE_UP_LP) * (Bose(deltaE_UP_LP))
Gamma_LP_UP_Hopfield = (2 / Nmol) * C_square * X_square * Jomg(deltaE_UP_LP) * (Bose(deltaE_UP_LP))

plt.plot(wc, Gamma_UP_D * 1000, "o-", color = 'red', label = r'$\Gamma_{+\to \{\mathrm{D}_k\}}$')#: Mixing Angle')
# plt.plot(wc, Gamma_UP_D_Hopfield * 1000, "o", color = 'red')#, label = 'UP -> D: Hopfield')

plt.plot(wc, Gamma_UP_LP * 1000, "o-", color = 'blue', label = r'$\Gamma_{+\to -}$')#: Mixing Angle')
# plt.plot(wc, Gamma_UP_LP_Hopfield * 1000, "o", color = 'blue')#, label = 'UP -> LP: Hopfield')

plt.plot(wc, Gamma_LP_D * 1000, "o-", color = 'orange', label = r'$\Gamma_{-\to \{\mathrm{D}_k\}}$') #: Mixing Angle')
# plt.plot(wc, Gamma_LP_D_Hopfield * 1000, "o", color = 'orange')#, label = 'LP -> D: Hopfield')

plt.plot(wc, Gamma_LP_UP * 1000, "o-", color = 'violet', label = r'$\Gamma_{-\to +}$')#: Mixing Angle')
# plt.plot(wc, Gamma_LP_UP_Hopfield * 1000, "o", color = 'violet')#, label = 'LP -> UP: Hopfield')

# ==============================================================================================

# RHS y-axis
ax = plt.gca()
x_major_locator = MultipleLocator(0.5)
x_minor_locator = MultipleLocator(0.1)
y_major_locator = MultipleLocator(5)
y_minor_locator = MultipleLocator(1)
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.tick_params(which = 'major', length = 15, pad = 10)
ax.tick_params(which = 'minor', length = 5)

x1_label = ax.get_xticklabels()
[x1_label_temp.set_fontname('Times New Roman') for x1_label_temp in x1_label]
y1_label = ax.get_yticklabels()
[y1_label_temp.set_fontname('Times New Roman') for y1_label_temp in y1_label]

plt.tick_params(which = 'both', direction = 'in', labelsize = 30)
plt.xlim(1.0, 3.0)
plt.ylim(y1, y2)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 15)
ax2.tick_params(which = 'minor', length = 5)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(y1, y2)

ax.set_xlabel(r'$\omega_\mathrm{c}\ (\mathrm{eV})$', size = 32)
ax.set_ylabel(r'Linewidth $(\mathrm{meV})$', size = 32)
ax.legend(loc = 'upper left', frameon = False, prop = font_legend)
ax.set_title(r"$N = 1$", size = 48)
plt.legend(title = '(a)', loc = "upper right", frameon = False, title_fontsize = legendsize)

# ==============================================================================================
#                           Fig s1a 
# ==============================================================================================
plt.subplot(1, 4, 2)

Nmol = 2
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

deltaE_LP_D = 0.5 * np.sqrt(Delta**2 + 4 * Rabi**2) + 0.5 * Delta
Gamma_LP_D = ((Nmol - 1) / Nmol) * (1 + np.cos(2 * theta)) * Jomg(deltaE_UP_D) * (Bose(deltaE_UP_D))
Gamma_LP_D_Hopfield = (2 * (Nmol - 1) / Nmol) * C_square * Jomg(deltaE_UP_D) * (Bose(deltaE_UP_D))

deltaE_LP_UP = np.sqrt(Delta**2 + 4 * Rabi**2)
Gamma_LP_UP = (1 / (2 * Nmol)) * (np.sin(2 * theta))**2 * Jomg(deltaE_UP_LP) * (Bose(deltaE_UP_LP))
Gamma_LP_UP_Hopfield = (2 / Nmol) * C_square * X_square * Jomg(deltaE_UP_LP) * (Bose(deltaE_UP_LP))

plt.plot(wc, Gamma_UP_D * 1000, "-", color = 'red', label = 'UP -> D')#: Mixing Angle')
plt.plot(wc, Gamma_UP_D_Hopfield * 1000, "o", color = 'red')#, label = 'UP -> D: Hopfield')

plt.plot(wc, Gamma_UP_LP * 1000, "-", color = 'blue', label = 'UP -> LP')#: Mixing Angle')
plt.plot(wc, Gamma_UP_LP_Hopfield * 1000, "o", color = 'blue')#, label = 'UP -> LP: Hopfield')

plt.plot(wc, Gamma_LP_D * 1000, "-", color = 'orange', label = 'LP -> D') #: Mixing Angle')
plt.plot(wc, Gamma_LP_D_Hopfield * 1000, "o", color = 'orange')#, label = 'LP -> D: Hopfield')

plt.plot(wc, Gamma_LP_UP * 1000, "-", color = 'violet', label = 'LP -> UP')#: Mixing Angle')
plt.plot(wc, Gamma_LP_UP_Hopfield * 1000, "o", color = 'violet')#, label = 'LP -> UP: Hopfield')

# ==============================================================================================

# RHS y-axis
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.tick_params(which = 'major', length = 15, pad = 10)
ax.tick_params(which = 'minor', length = 5)

x1_label = ax.get_xticklabels()
[x1_label_temp.set_fontname('Times New Roman') for x1_label_temp in x1_label]
y1_label = ax.get_yticklabels()
[y1_label_temp.set_fontname('Times New Roman') for y1_label_temp in y1_label]

plt.tick_params(which = 'both', direction = 'in', labelsize = 30)
plt.xlim(1.0, 3.0)
plt.ylim(y1, y2)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 15)
ax2.tick_params(which = 'minor', length = 5)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(y1, y2)

ax.set_xlabel(r'$\omega_\mathrm{c}\ (\mathrm{eV})$', size = 32)
# ax.set_ylabel(r'Linewidth $(\mathrm{meV})$', size = 32)
# ax.legend(loc = 'upper left', frameon = False, prop = font_legend)
ax.set_title(r"$N = 2$", size = 48)
plt.legend(title = '(b)', loc = "upper right", frameon = False, title_fontsize = legendsize)

# ==============================================================================================
#                           Fig s1a 
# ==============================================================================================
plt.subplot(1, 4, 3)

Nmol = 4
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

deltaE_LP_D = 0.5 * np.sqrt(Delta**2 + 4 * Rabi**2) + 0.5 * Delta
Gamma_LP_D = ((Nmol - 1) / Nmol) * (1 + np.cos(2 * theta)) * Jomg(deltaE_UP_D) * (Bose(deltaE_UP_D))
Gamma_LP_D_Hopfield = (2 * (Nmol - 1) / Nmol) * C_square * Jomg(deltaE_UP_D) * (Bose(deltaE_UP_D))

deltaE_LP_UP = np.sqrt(Delta**2 + 4 * Rabi**2)
Gamma_LP_UP = (1 / (2 * Nmol)) * (np.sin(2 * theta))**2 * Jomg(deltaE_UP_LP) * (Bose(deltaE_UP_LP))
Gamma_LP_UP_Hopfield = (2 / Nmol) * C_square * X_square * Jomg(deltaE_UP_LP) * (Bose(deltaE_UP_LP))

plt.plot(wc, Gamma_UP_D * 1000, "-", color = 'red', label = 'UP -> D')#: Mixing Angle')
plt.plot(wc, Gamma_UP_D_Hopfield * 1000, "o", color = 'red')#, label = 'UP -> D: Hopfield')

plt.plot(wc, Gamma_UP_LP * 1000, "-", color = 'blue', label = 'UP -> LP')#: Mixing Angle')
plt.plot(wc, Gamma_UP_LP_Hopfield * 1000, "o", color = 'blue')#, label = 'UP -> LP: Hopfield')

plt.plot(wc, Gamma_LP_D * 1000, "-", color = 'orange', label = 'LP -> D') #: Mixing Angle')
plt.plot(wc, Gamma_LP_D_Hopfield * 1000, "o", color = 'orange')#, label = 'LP -> D: Hopfield')

plt.plot(wc, Gamma_LP_UP * 1000, "-", color = 'violet', label = 'LP -> UP')#: Mixing Angle')
plt.plot(wc, Gamma_LP_UP_Hopfield * 1000, "o", color = 'violet')#, label = 'LP -> UP: Hopfield')

# ==============================================================================================

# RHS y-axis
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.tick_params(which = 'major', length = 15, pad = 10)
ax.tick_params(which = 'minor', length = 5)

x1_label = ax.get_xticklabels()
[x1_label_temp.set_fontname('Times New Roman') for x1_label_temp in x1_label]
y1_label = ax.get_yticklabels()
[y1_label_temp.set_fontname('Times New Roman') for y1_label_temp in y1_label]

plt.tick_params(which = 'both', direction = 'in', labelsize = 30)
plt.xlim(1.0, 3.0)
plt.ylim(y1, y2)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 15)
ax2.tick_params(which = 'minor', length = 5)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(y1, y2)

ax.set_xlabel(r'$\omega_\mathrm{c}\ (\mathrm{eV})$', size = 32)
# ax.set_ylabel(r'Linewidth $(\mathrm{meV})$', size = 32)
# ax.legend(loc = 'upper left', frameon = False, prop = font_legend)
ax.set_title(r"$N = 4$", size = 48)
plt.legend(title = '(c)', loc = "upper right", frameon = False, title_fontsize = legendsize)

# ==============================================================================================
#                           Fig s1a 
# ==============================================================================================
plt.subplot(1, 4, 4)

Nmol = 8
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

deltaE_LP_D = 0.5 * np.sqrt(Delta**2 + 4 * Rabi**2) + 0.5 * Delta
Gamma_LP_D = ((Nmol - 1) / Nmol) * (1 + np.cos(2 * theta)) * Jomg(deltaE_UP_D) * (Bose(deltaE_UP_D))
Gamma_LP_D_Hopfield = (2 * (Nmol - 1) / Nmol) * C_square * Jomg(deltaE_UP_D) * (Bose(deltaE_UP_D))

deltaE_LP_UP = np.sqrt(Delta**2 + 4 * Rabi**2)
Gamma_LP_UP = (1 / (2 * Nmol)) * (np.sin(2 * theta))**2 * Jomg(deltaE_UP_LP) * (Bose(deltaE_UP_LP))
Gamma_LP_UP_Hopfield = (2 / Nmol) * C_square * X_square * Jomg(deltaE_UP_LP) * (Bose(deltaE_UP_LP))

plt.plot(wc, Gamma_UP_D * 1000, "-", color = 'red', label = 'UP -> D')#: Mixing Angle')
plt.plot(wc, Gamma_UP_D_Hopfield * 1000, "o", color = 'red')#, label = 'UP -> D: Hopfield')

plt.plot(wc, Gamma_UP_LP * 1000, "-", color = 'blue', label = 'UP -> LP')#: Mixing Angle')
plt.plot(wc, Gamma_UP_LP_Hopfield * 1000, "o", color = 'blue')#, label = 'UP -> LP: Hopfield')

plt.plot(wc, Gamma_LP_D * 1000, "-", color = 'orange', label = 'LP -> D') #: Mixing Angle')
plt.plot(wc, Gamma_LP_D_Hopfield * 1000, "o", color = 'orange')#, label = 'LP -> D: Hopfield')

plt.plot(wc, Gamma_LP_UP * 1000, "-", color = 'violet', label = 'LP -> UP')#: Mixing Angle')
plt.plot(wc, Gamma_LP_UP_Hopfield * 1000, "o", color = 'violet')#, label = 'LP -> UP: Hopfield')

# ==============================================================================================

# RHS y-axis
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.tick_params(which = 'major', length = 15, pad = 10)
ax.tick_params(which = 'minor', length = 5)

x1_label = ax.get_xticklabels()
[x1_label_temp.set_fontname('Times New Roman') for x1_label_temp in x1_label]
y1_label = ax.get_yticklabels()
[y1_label_temp.set_fontname('Times New Roman') for y1_label_temp in y1_label]

plt.tick_params(which = 'both', direction = 'in', labelsize = 30)
plt.xlim(1.0, 3.0)
plt.ylim(y1, y2)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 15)
ax2.tick_params(which = 'minor', length = 5)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(y1, y2)

ax.set_xlabel(r'$\omega_\mathrm{c}\ (\mathrm{eV})$', size = 32)
# ax.set_ylabel(r'Linewidth $(\mathrm{meV})$', size = 32)
# ax.legend(loc = 'upper left', frameon = False, prop = font_legend)
ax.set_title(r"$N = 8$", size = 48)
plt.legend(title = '(d)', loc = "upper right", frameon = False, title_fontsize = legendsize)


plt.savefig("Fig_S1.pdf", bbox_inches='tight')