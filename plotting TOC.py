import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator, tick_params
fig, ax = plt.subplots()
plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams["font.family"] = "Helvetica"

# ==============================================================================================
#                                       Global Parameters     
# ==============================================================================================
lw = 3.0
legendsize = 48         # size for legend
font_legend = {'family':'Times New Roman', 'weight': 'roman', 'size': 12}
# axis label size
lsize = 30             
txtsize = 32
# tick length
lmajortick = 15
lminortick = 5
legend_x, legend_y = - 0.2, 1.03
transparency = .4
y1, y2 = 0.018, 0.08

unitlen = 7
fig = plt.figure(figsize=(1.0 * unitlen, 0.75 * unitlen), dpi = 128)
plt.subplots_adjust(hspace = 0.25, wspace = 0.5)

# whether plot UP or not
pltup = True

def Jomg(w):
    lam = 30 / 1000
    gam = 24.8 / 1000
    return 2 * lam * gam * w / (w**2 + gam**2)

def Bose(w):
    beta = 1052.6 / (27.2114)
    return 1. / (np.exp(beta * w) - 1)

Gammac = 8.83 * 5
Nmol = 8
data = np.loadtxt("width_N=%s_Gammac=%smev.txt" % (Nmol, Gammac), dtype = float)

# plot under detuning
plt.plot(data[:, 0], data[:, 1], 'o', markersize = 10, linewidth = 1, markerfacecolor = 'white', color = 'blue', label = "LP (HEOM)")
plt.plot(data[:, 0], data[:, 2], 'o', markersize = 10, linewidth = 1, markerfacecolor = 'white', color = 'red', label = "UP (HEOM)")

Gammac = Gammac / 1000      # convert meV to eV
Gammae = 0.0765
wc = np.linspace(1.0, 3.0, 1000)
w0 = 2.0 + 0.03
Delta = wc - w0
gc = 0.0681
Rabi = gc * np.sqrt(Nmol)

# Hopfield coefficients
X_square = 0.5 * (1 + Delta / np.sqrt(Delta**2 + 4 * Rabi**2))
C_square = 0.5 * (1 - Delta / np.sqrt(Delta**2 + 4 * Rabi**2))

# light-matter mixing angle
theta = 0.5 * np.arctan(2 * Rabi / Delta)
for i in range(len(theta)):
    if(theta[i] < 0.0):
        theta[i] = np.pi / 2 - np.abs(theta[i])
fN = (np.sin(2 * theta))**2 / Nmol + (np.cos(2 * theta))**2

deltaE_UP_D = 0.5 * np.sqrt(Delta**2 + 4 * Rabi**2) + 0.5 * Delta
# Gamma_UP_D = ((Nmol - 1) / Nmol) * (1 + np.cos(2 * theta)) * Jomg(deltaE_UP_D) * (1 + Bose(deltaE_UP_D))
Gamma_UP_D = ((Nmol - 1) / Nmol) * 2 * X_square * Jomg(deltaE_UP_D) * (1 + Bose(deltaE_UP_D))
deltaE_UP_LP = np.sqrt(Delta**2 + 4 * Rabi**2)
Gamma_UP_LP = (1 / (2 * Nmol)) * (np.sin(2 * theta))**2 * Jomg(deltaE_UP_LP) * (1 + Bose(deltaE_UP_LP))

# Corrected polariton line width
Gamma_UP = (X_square * Gammac) + (C_square**2 / (X_square + C_square**2 / Nmol) * Gammae / Nmol) + Gamma_UP_D + Gamma_UP_LP
Gamma_LP = (C_square * Gammac) + (X_square**2 / (C_square + X_square**2 / Nmol) * Gammae / Nmol)

Gamma_bare_UP = (X_square * Gammac) + (C_square * Gammae)
Gamma_bare_LP = (C_square * Gammac) + (X_square * Gammae)

plt.plot(wc, Gamma_bare_LP, "--", linewidth = lw, color = 'k', alpha = transparency)
plt.plot(wc, Gamma_bare_UP, "--", linewidth = lw, color = 'k', alpha = transparency)
plt.plot(wc, Gamma_LP, "-", linewidth = 2.5, color = 'blue', label = "LP (Theory)")

if pltup:
    plt.plot(wc, Gamma_UP, "-", linewidth = 2.5, color = 'red', label = "UP (Theory)")

# ==============================================================================================

# RHS y-axis
ax = plt.gca()
x_major_locator = MultipleLocator(0.5)
x_minor_locator = MultipleLocator(0.1)
y_major_locator = MultipleLocator(0.02)
y_minor_locator = MultipleLocator(0.01)
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
ax.set_ylabel(r'Linewidth $(\mathrm{eV})$', size = 32)
ax.legend(loc = 'lower left', frameon = False, prop = font_legend)
# ax.set_title(r"$N = 8$", size = 48)
# plt.legend(title = '(g)', bbox_to_anchor = (legend_x, legend_y), frameon = False, title_fontsize = legendsize)



plt.savefig("Fig_TOC.png", bbox_inches='tight')
