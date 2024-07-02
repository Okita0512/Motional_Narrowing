import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator, tick_params
fig, ax = plt.subplots()
from numpy import polyfit, poly1d
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
legend_x, legend_y = - 0.135, 1.03
transparency = .4

unitlen = 7
fig = plt.figure(figsize=(3.0 * unitlen, 1.8 * unitlen), dpi = 128)
plt.subplots_adjust(wspace = 0.375, hspace = 0.26)

# ==============================================================================================
#                           Fig 2a 
# ==============================================================================================

plt.subplot(2, 3, 1)

# determine the normalization factor
data2 = np.loadtxt("./N=8/resp1st_im.w", dtype = float)
norm = np.max(data2) / 8

#=============

data = np.loadtxt("./N=1_out/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=1_out/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = "Molecule")

data = np.loadtxt("./N=1/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=1/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = r"$N$ = 1")

data = np.loadtxt("./N=2/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=2/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / (2 * norm), '-', linewidth = lw, label = r"$N$ = 2")

data = np.loadtxt("./N=4/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=4/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / (4 * norm), '-', linewidth = lw, label = r"$N$ = 4")

data = np.loadtxt("./N=8/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=8/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / (8 * norm), '-', linewidth = lw, label = r"$N$ = 8")

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
ax.legend(loc = 'upper center', frameon = False, prop = font_legend)

plt.legend(title = '(a)', bbox_to_anchor = (legend_x, legend_y), frameon = False, title_fontsize = legendsize)

# ==============================================================================================
#                           Fig 2b
# ==============================================================================================

plt.subplot(2, 3, 2)

# determine the normalization factor
data2 = np.loadtxt("./N=8/resp1st_im.w", dtype = float)
norm = np.max(data2) / 8

#=============

data = np.loadtxt("./N=1_out/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=1_out/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = "Molecule")

data = np.loadtxt("./N=1_fixing_Rabi_8/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=1_fixing_Rabi_8/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = r"$N$ = 1")

data = np.loadtxt("./N=2_fixing_Rabi_8/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=2_fixing_Rabi_8/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / (2 * norm), '-', linewidth = lw, label = r"$N$ = 2")

data = np.loadtxt("./N=4_fixing_Rabi_8/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=4_fixing_Rabi_8/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / (4 * norm), '-', linewidth = lw, label = r"$N$ = 4")

data = np.loadtxt("./N=8/resp1st.w1", dtype = float)
data2 = np.loadtxt("./N=8/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / (8 * norm), '-', linewidth = lw, label = r"$N$ = 8")

#=============

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
ax.legend(loc = 'upper center', frameon = False, prop = font_legend)

plt.legend(title = '(b)', bbox_to_anchor = (legend_x, legend_y), frameon = False, title_fontsize = legendsize)

# ==============================================================================================
#                           Fig 2c
# ==============================================================================================

plt.subplot(2, 3, 3)

def rsquare(x, y, degree):
    results = {}

    coeffs = polyfit(x, y, degree)

     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y) / len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot

    return results

data = np.loadtxt("./width_scaling_fast.txt", dtype = float)

plt.plot(1. / data[:, 0], 1000 * data[:, 2], 'o', markersize = 10, markerfacecolor = 'white', color = 'red', label = r"$|+\rangle$")
plt.plot(1. / data[:, 0], 1000 * data[:, 1], 'o', markersize = 10, markerfacecolor = 'white', color = 'blue', label = r"$|-\rangle$")

coeff_LP = polyfit(1. / data[:, 0], 1000 * data[:, 1], 1)
coeff_UP = polyfit(1. / data[:, 0], 1000 * data[:, 2], 1)

rsquare_1 = rsquare(1. / data[:, 0], 1000 * data[:, 1], degree = 1)
rsquare_2 = rsquare(1. / data[:, 0], 1000 * data[:, 2], degree = 1)
print("1/N fit for LP: ", rsquare_1)
print("1/N fit for UP: ", rsquare_2)

plt.text(0.3, 8, r"$R^2$ $\approx$ %.4f" % rsquare_1['determination'], size = 15, color = 'blue')
plt.text(0.1, 15, r"$R^2$ $\approx$ %.4f" % rsquare_2['determination'], size = 15, color = 'red')

N0 = np.linspace(1, 10000)
plt.plot(1. / N0, coeff_UP[0] / N0 + coeff_UP[1], '--', color = 'red')
plt.plot(1. / N0, coeff_LP[0] / N0 + coeff_LP[1], '--', color = 'blue')

#=============

# RHS y-axis
ax = plt.gca()
x_major_locator = MultipleLocator(0.2)
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
plt.xlim(0.0, 0.6)
plt.ylim(0, 20)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 15)
ax2.tick_params(which = 'minor', length = 5)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(0, 20)

ax.set_xlabel(r'$1 / N$', size = 32)
ax.set_ylabel(r'Width (meV)', size = 32)
ax.legend(loc = 'lower right', frameon = False, prop = font_legend)

plt.legend(title = '(c)', bbox_to_anchor = (legend_x, legend_y), frameon = False, title_fontsize = legendsize)

# ==============================================================================================
#                           Fig 2a 
# ==============================================================================================

plt.subplot(2, 3, 4)

# determine the normalization factor
data2 = np.loadtxt("./scaling_fixing_gc_slow/8/resp1st_im.w", dtype = float)
norm = np.max(data2) / 8

#=============

data = np.loadtxt("./scaling_fixing_gc_slow/outside/resp1st.w1", dtype = float)
data2 = np.loadtxt("./scaling_fixing_gc_slow/outside/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = "Molecule")

data = np.loadtxt("./scaling_fixing_gc_slow/1/resp1st.w1", dtype = float)
data2 = np.loadtxt("./scaling_fixing_gc_slow/1/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = r"$N$ = 1")

data = np.loadtxt("./scaling_fixing_gc_slow/2/resp1st.w1", dtype = float)
data2 = np.loadtxt("./scaling_fixing_gc_slow/2/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / (2 * norm), '-', linewidth = lw, label = r"$N$ = 2")

data = np.loadtxt("./scaling_fixing_gc_slow/4/resp1st.w1", dtype = float)
data2 = np.loadtxt("./scaling_fixing_gc_slow/4/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / (4 * norm), '-', linewidth = lw, label = r"$N$ = 4")

data = np.loadtxt("./scaling_fixing_gc_slow/8/resp1st.w1", dtype = float)
data2 = np.loadtxt("./scaling_fixing_gc_slow/8/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / (8 * norm), '-', linewidth = lw, label = r"$N$ = 8")

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
ax.legend(loc = 'upper center', frameon = False, prop = font_legend)

plt.legend(title = '(d)', bbox_to_anchor = (legend_x, legend_y), frameon = False, title_fontsize = legendsize)

# ==============================================================================================
#                           Fig 2b
# ==============================================================================================

plt.subplot(2, 3, 5)

# determine the normalization factor
# data2 = np.loadtxt("./scaling_fixing_gc_slow/N=8/resp1st_im.w", dtype = float)
# norm = np.max(data2) / 8

#=============

data = np.loadtxt("./scaling_fixing_gc_slow/outside/resp1st.w1", dtype = float)
data2 = np.loadtxt("./scaling_fixing_gc_slow/outside/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = "Molecule")

data = np.loadtxt("./scaling_slow/1/resp1st.w1", dtype = float)
data2 = np.loadtxt("./scaling_slow/1/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / norm, '-', linewidth = lw, label = r"$N$ = 1")

data = np.loadtxt("./scaling_slow/2/resp1st.w1", dtype = float)
data2 = np.loadtxt("./scaling_slow/2/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / (2 * norm), '-', linewidth = lw, label = r"$N$ = 2")

data = np.loadtxt("./scaling_slow/4/resp1st.w1", dtype = float)
data2 = np.loadtxt("./scaling_slow/4/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / (4 * norm), '-', linewidth = lw, label = r"$N$ = 4")

data = np.loadtxt("./scaling_slow/8/resp1st.w1", dtype = float)
data2 = np.loadtxt("./scaling_slow/8/resp1st_im.w", dtype = float)
plt.plot(data * conv, 100 * data2 / (8 * norm), '-', linewidth = lw, label = r"$N$ = 8")

#=============

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
ax.legend(loc = 'upper center', frameon = False, prop = font_legend)

plt.legend(title = '(e)', bbox_to_anchor = (legend_x, legend_y), frameon = False, title_fontsize = legendsize)

# ==============================================================================================
#                           Fig 2c
# ==============================================================================================

plt.subplot(2, 3, 6)

def rsquare(x, y, degree):
    results = {}

    coeffs = polyfit(x, y, degree)

     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y) / len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot

    return results

data = np.loadtxt("./width_scaling_slow.txt", dtype = float)

# plt.plot(1. / data[:, 0], data[:, 1], 'o', markersize = 10, color = 'blue', label = "LP")
# plt.plot(1. / data[:, 0], data[:, 2], 'o', markersize = 10, color = 'red', label = "UP")

plt.plot(1. / np.sqrt(data[:, 0]), 1000 * data[:, 2], 'o', markersize = 10, markerfacecolor = 'white', color = 'red', label = r"$|+\rangle$")
plt.plot(1. / np.sqrt(data[:, 0]), 1000 * data[:, 1], 'o', markersize = 10, markerfacecolor = 'white', color = 'blue', label = r"$|-\rangle$")

coeff_LP = polyfit(1. / np.sqrt(data[:, 0]), 1000 * data[:, 1], 1)
coeff_UP = polyfit(1. / np.sqrt(data[:, 0]), 1000 * data[:, 2], 1)

rsquare_1 = rsquare(1. / np.sqrt(data[:, 0]), 1000 * data[:, 1], degree = 1)
rsquare_2 = rsquare(1. / np.sqrt(data[:, 0]), 1000 * data[:, 2], degree = 1)
print("1/sqrt(N) fit for LP: ", rsquare_1)
print("1/sqrt(N) fit for UP: ", rsquare_2)

# rsquare_3 = rsquare(1. / data[:, 0], data[:, 1], degree = 1)
# rsquare_4 = rsquare(1. / data[:, 0], data[:, 2], degree = 1)
# print("1/N fit for LP: ", rsquare_3)
# print("1/N fit for UP: ", rsquare_4)

N0 = np.linspace(1, 10000)
plt.plot(1. / N0, coeff_UP[0] / N0 + coeff_UP[1], '--', color = 'red')
plt.plot(1. / N0, coeff_LP[0] / N0 + coeff_LP[1], '--', color = 'blue')

plt.text(0.7, 26, r"$R^2$ $\approx$ %.4f" % rsquare_1['determination'], size = 15, color = 'blue')
plt.text(0.28, 32, r"$R^2$ $\approx$ %.4f" % rsquare_2['determination'], size = 15, color = 'red')

#=============

# RHS y-axis
ax = plt.gca()
x_major_locator = MultipleLocator(0.2)
x_minor_locator = MultipleLocator(0.1)
y_major_locator = MultipleLocator(10)
y_minor_locator = MultipleLocator(2)
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
plt.xlim(0.0, 1.05)
plt.ylim(0, 50)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 15)
ax2.tick_params(which = 'minor', length = 5)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(0, 50)

ax.set_xlabel(r'$1 / \sqrt{N}$', size = 32)
ax.set_ylabel(r'Width (meV)', size = 32)
ax.legend(loc = 'lower right', frameon = False, prop = font_legend)

plt.legend(title = '(f)', bbox_to_anchor = (legend_x, legend_y), frameon = False, title_fontsize = legendsize)



plt.savefig("Fig_1.pdf", bbox_inches='tight')

