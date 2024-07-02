import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the Lorentzian function
def lorentzian(x, amp, center, width):
    return (amp / np.pi) * (width/2) / ((width/2)**2 + (x - center)**2)

# Define the sum of two Lorentzian functions
def double_lorentzian(x, amp1, center1, width1, amp2, center2, width2):
    return lorentzian(x, amp1, center1, width1) + lorentzian(x, amp2, center2, width2)

UP_width = []
LP_width = []
Nmol = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

for N in Nmol:

    # Generate some sample data
    x_data = np.loadtxt('./scaling_fast/%s/resp1st.w1' % (N))
    y_data = np.loadtxt('./scaling_fast/%s/resp1st_im.w' % (N))

    x_data = x_data * 27.2114
    y_data = y_data / max(y_data)

    # Set the rough positions of the peaks and the minimum separation
    w0 = 2.0 + 0.03
    wc = w0
    gc = 0.0681 * np.sqrt(N)
    tolerance = 0.00
    rough_position1 = 0.5 * (wc + w0) - 0.5 * np.sqrt((wc - w0)**2 + 4 * gc**2) - tolerance
    rough_position2 = 0.5 * (wc + w0) + 0.5 * np.sqrt((wc - w0)**2 + 4 * gc**2) + tolerance

    # Adjust the initial guess to ensure minimum separation
    initial_guess = [
        0.5,  # Amplitude 1
        rough_position1,  # Center 1
        0.05,  # Width 1
        0.5,  # Amplitude 2
        rough_position2,  # Center 2
        0.05  # Width 2
    ]

    # Fit the data using curve_fit
    params, covariance = curve_fit(double_lorentzian, x_data, y_data, p0 = initial_guess, bounds = (0, [1.0, 3.5, 0.2, 1.0, 3.5, 0.2]))

    if params[1] > params[4]:
        UPw = np.abs(params[1])
        LPw = np.abs(params[4])

    elif params[1] < params[4]: 
        UPw = np.abs(params[4])
        LPw = np.abs(params[1])
    
    difference_array1 = np.absolute(x_data - UPw)
    index_UP = difference_array1.argmin()
    difference_array2 = np.absolute(x_data - LPw)
    index_LP = difference_array2.argmin()
    print(index_LP, index_UP)

    plt.vlines([x_data[index_UP]], 0, y_data[index_UP], linestyles = 'dashed', colors = ["red"])
    plt.vlines([x_data[index_LP]], 0, y_data[index_LP], linestyles = 'dashed', colors = ["g"])

    difference_array3 = np.absolute(y_data[index_UP:] - 0.5 * y_data[index_UP])
    UP_width_Right = x_data[index_UP + difference_array3.argmin()]

    difference_array4 = np.absolute(y_data[int(0.5 * (index_LP + index_UP)) : index_UP] - 0.5 * y_data[index_UP])
    UP_width_Left = x_data[int(0.5 * (index_LP + index_UP)) + difference_array4.argmin()]

    difference_array5 = np.absolute(y_data[index_LP : int(0.5 * (index_LP + index_UP))] - 0.5 * y_data[index_LP])
    LP_width_Right = x_data[index_LP + difference_array5.argmin()]

    difference_array6 = np.absolute(y_data[:index_LP] - 0.5 * y_data[index_LP])
    LP_width_Left = x_data[difference_array6.argmin()]

    print(UP_width_Left, UP_width_Right)
    print(LP_width_Left, LP_width_Right)
    UP_width.append(UP_width_Right - UP_width_Left)
    LP_width.append(LP_width_Right - LP_width_Left)

    plt.hlines([0.5 * y_data[index_UP]], UP_width_Left, UP_width_Right, linestyles = 'dashed', colors = ["red"])
    plt.hlines([0.5 * y_data[index_LP]], LP_width_Left, LP_width_Right, linestyles = 'dashed', colors = ["g"])

    plt.scatter(x_data, y_data, label='Original')
#    plt.plot(x_data, double_lorentzian(x_data, *params), 'r-', label='Fit')
    plt.legend()
    plt.show()

bunch = np.zeros((len(Nmol), 3), dtype = float)
bunch[:, 0] = Nmol
bunch[:, 1] = LP_width
bunch[:, 2] = UP_width

np.savetxt("width_scaling_fast.txt", bunch)