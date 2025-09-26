 # SImulation of fluid flow

import numpy as np

#Polynomial definition
poly_p = np.array([-4, 7, -3, 9])

#polynomial derivative
p_der = np.polyder(poly_p)

#Confiramtion of correct polynomial derivative
print('Derivative of polynomial =', p_der )

#Analytical result at x = 0
p_der_eval = np.polyval(p_der, 0.0)
print('Theorotical Deriative =', p_der_eval)

#Numerical calculations using FDM (forward)
x_0 = 0.0
h = np.float64(0.25)
forward_difference = (np.polyval(poly_p, x_0 + h) - np.polyval(poly_p,x_0)) / h
print('Forward Derivative =', forward_difference)

#Numerical calculations using FDM (backward)
backward_difference = (np.polyval(poly_p, x_0) - np.polyval(poly_p, x_0 - h)) / h
print('Backward Derivative =', backward_difference)

#Numerical calculation using FDM (Central difference)
Central_difference = (np.polyval(poly_p, x_0 + h) - np.polyval(poly_p, x_0 - h)) / (2*h)
print('Central difference =', Central_difference)
