# Getting the numpy module
import numpy as np

# Getting matplot library
import matplotlib.pyplot as plt

#Number of grid points
N = 101

#Domain size
L = 1

#Corresponding grid space
h = np.float64(L / (N - 1))

# Thermal conductivity
k = 0.1;

#Area
A = 0.001;

#Iterations
iterations = 0

# Initialise the temperture field
T = np.zeros(N)
T[N-1] = 1.

#Initialising the iterated temperature
T_new = np.zeros(N)
T_new[N-1] = 1.

#Error related variables
epsilon = 1.E-05
numerical_error = 1

#Checking the error tolerance
while numerical_error > epsilon:
    for i in range(1,N-1):
        a_E = np.float64(k*A/h)
        a_W = np.float64(k*A/h)
        if i==0:
            a_W = 0
        a_P = a_W + a_E
        T_new[i] = (a_E*T[i+1] + a_W*T[i-1]) / a_P
    
    #recalculating the error'
    numerical_error = np.sum(np.abs(T_new-T))
            
    #iterations advancement
    iterations = iterations + 1
    T = T_new.copy()
        
           
#Plotting the results
plt.figure()

#Defining the position vector from the indices
x_dom = np.arange(N) * h

#Plotting the variation with customisation
plt.plot(x_dom, T, 'gx--', linewidth=2, markersize = 10)

#Displaying the grid lines
plt.grid(True, color = 'k')

#Labelling and providing a title
plt.xlabel("Position", size = 20)
plt.ylabel("Temperature", size = 20)        
plt.title("T(x)")

#Show that plot on screen'
plt.show()    