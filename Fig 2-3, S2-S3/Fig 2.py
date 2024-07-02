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
font_legend = {'family':'Times New Roman', 'weight': 'roman', 'size': 18}
# axis label size
lsize = 30             
txtsize = 12
# tick length
lmajortick = 15
lminortick = 5
legend_x, legend_y = - 0.19, 0.78
transparency = .4
y1, y2 = - 0.5, 10

unitlen = 7
fig = plt.figure(figsize=(1.0 * unitlen, 0.8 * unitlen), dpi = 128)
plt.subplots_adjust(hspace = 0.3)

# Define the Lorentzian function
def lorentzian(x, amp, center, width):
    return (amp / np.pi) * (width/2) / ((width/2)**2 + (x - center)**2)

# Define the sum of two Lorentzian functions
def double_lorentzian(x, amp1, center1, width1, amp2, center2, width2):
    return lorentzian(x, amp1, center1, width1) + lorentzian(x, amp2, center2, width2)

UP_width = []
LP_width = []
omegac = [1.8,1.85,1.9,1.95,2.0,2.05,2.1,2.15,2.2]

# filename = "N=1_Gammac=8.83mev"
# filename = "N=1_Gammac=44.15mev"
# filename = "N=1_Gammac=88.3mev"
# filename = "N=2_Gammac=8.83mev"
# filename = "N=2_Gammac=44.15mev"
# filename = "N=2_Gammac=88.3mev"
# filename = "N=4_Gammac=8.83mev"
filename = "N=4_Gammac=44.15mev"
# filename = "N=4_Gammac=88.3mev"
# filename = "N=8_Gammac=8.83mev"
# filename = "N=8_Gammac=44.15mev"
# filename = "N=8_Gammac=88.3mev"

count = 0
for wc in omegac:

    # Generate some sample data
    x_data = np.loadtxt('./%s/%s/resp1st.w1' % (filename, wc))
    y_data = np.loadtxt('./%s/%s/resp1st_im.w' % (filename, wc))

    x_data = x_data * 27.2114
    y_data = y_data / max(y_data)

    plt.plot(x_data, y_data + count, "-", linewidth = lw)
    wc = wc + 0.03
    plt.text(1.53, 0.15 + count, "$\omega_\mathrm{c}=$%.2f eV" %wc, size = txtsize)
    count += 1.0 + 0.1
    
# ==============================================================================================

# RHS y-axis
ax = plt.gca()
x_major_locator = MultipleLocator(0.2)
x_minor_locator = MultipleLocator(0.1)
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.tick_params('x', which = 'major', length = 15, pad = 10)
ax.tick_params('x', which = 'minor', length = 5)
ax.tick_params('y', which = 'major', length = 0)
ax.tick_params('y', which = 'minor', length = 0)

x1_label = ax.get_xticklabels()
[x1_label_temp.set_fontname('Times New Roman') for x1_label_temp in x1_label]
# y1_label = ax.get_yticklabels()
# [y1_label_temp.set_fontname('Times New Roman') for y1_label_temp in y1_label]

plt.tick_params('x', which = 'both', direction = 'in', labelsize = 30)
plt.tick_params('y', which = 'both', direction = 'in', labelsize = 0)
plt.xlim(1.5, 2.5)
plt.ylim(y1, y2)

ax.set_xlabel(r'$\omega_\mathrm{c}\ (\mathrm{eV})$', size = 32)

# plt.legend(title = '(b)', bbox_to_anchor = (legend_x, legend_y), frameon = False, title_fontsize = legendsize)




plt.savefig("Fig 2.pdf", bbox_inches='tight')