
import numpy as np
import matplotlib.pyplot as plt

# Number of grid posints
N = 11

#Length of the domain
L = 1

#Corresponding grid spacing
h = np.float64(L/(N-1))

#Iteration number
iterations = 0

#Initialising the temperature
T = np.zeros(N) # Initialising the array with zeros having N grid points
T[N-1] = 1.
#Initialising the iterated temperature
T_new = np.zeros(N)
T_new[N-1] = 1.

#Error related variable
epsilon = 1.E-8
numerical_error = 1

# Checking the error tolerance
while numerical_error > epsilon:
    for i in range(1,N-1):
        T_new[i] = 0.5*(T[i-1] + T[i+1])
        
        #Resetting the numerical error and recalculate 
        numerical_error = 0
    for i in range(1, N-1):
        numerical_error = numerical_error + abs(T[i] - T_new[i])
        
        #Iterations advancement and reassignment
    iterations = iterations + 1
    T = T_new.copy()
        
#Plotting the results

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
        
        
        

