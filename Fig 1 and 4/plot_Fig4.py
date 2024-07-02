import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator, tick_params
fig, ax = plt.subplots()
plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams["font.family"] = "Helvetica"

# ==============================================================================================
#                                       Global Parameters     
# ==============================================================================================

conv = 27.211397                            # 1 a.u. = 27.211397 eV
fs_to_au = 41.341374575751                  # 1 fs = 41.341 a.u.
cm_to_au = 4.556335e-06                     # 1 cm^-1 = 4.556335e-06 a.u.
au_to_K = 3.1577464e+05                     # 1 au = 3.1577464e+05 K
kcal_to_au = 1.5936e-03                     # 1 kcal/mol = 1.5936e-3 a.u.

lw = 3.0
legendsize = 48         # size for legend
font_legend = {'family':'Times New Roman', 'weight': 'roman', 'size': 18}
# axis label size
lsize = 30             
txtsize = 32
# tick length
lmajortick = 15
lminortick = 5
legend_x, legend_y = - 0.115, 1.03
transparency = .4

unitlen = 7
fig = plt.figure(figsize=(1.0 * unitlen, 0.8 * unitlen), dpi = 128)

# ==============================================================================================
#                           Fig 2a 
# ==============================================================================================

# determine the normalization factor
data2 = np.loadtxt("./N=4/resp1st_im.w", dtype = float)
norm = np.max(data2)

#=============

data = np.loadtxt("./N=4/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=4/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = r"Lossless")

data = np.loadtxt("./N=4_Gammac=8.83mev/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=4_Gammac=8.83mev/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = r"$\Gamma_\mathrm{c}=$ 8.83 meV")

data = np.loadtxt("./N=4_Gammac=17.66mev/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=4_Gammac=17.66mev/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = r"17.66 meV")

data = np.loadtxt("./N=4_Gammac=44.15mev/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=4_Gammac=44.15mev/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = r"44.15 meV")

data = np.loadtxt("./N=4_Gammac=88.3mev/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=4_Gammac=88.3mev/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = r"88.3 meV")

data = np.loadtxt("./N=4_Gammac=176.6mev/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=4_Gammac=176.6mev/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = r"176.6 meV")

data = np.loadtxt("./N=4_Gammac=441.5mev/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=4_Gammac=441.5mev/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = r"441.5 meV")

#=============

# RHS y-axis
ax = plt.gca()
x_major_locator = MultipleLocator(0.1)
x_minor_locator = MultipleLocator(0.02)
y_major_locator = MultipleLocator(20)
y_minor_locator = MultipleLocator(10)
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
plt.xlim(1.75, 2.3)
plt.ylim(-5, 105)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 15)
ax2.tick_params(which = 'minor', length = 5)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(-5, 105)

ax.set_xlabel(r'$\omega_\mathrm{c}\ (\mathrm{eV})$', size = 32)
ax.set_ylabel(r'Intensity (%)', size = 32)
ax.legend(loc = 'upper center', bbox_to_anchor = (0.55, 1.0), frameon = False, prop = font_legend, ncol = 1)

# plt.legend(title = '(a)', bbox_to_anchor = (legend_x, legend_y), frameon = False, title_fontsize = legendsize)




plt.savefig("Fig_4.pdf", bbox_inches='tight')

