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


