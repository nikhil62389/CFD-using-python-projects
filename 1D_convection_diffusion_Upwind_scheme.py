import numpy as np
import matplotlib.pyplot as plt

N = 51

L = 1.0

h = np.float64(L/(N-1))

#velocity
u = np.float64(50)

#Density
rho = np.float64(1.0)

#Diffusivity
gamma = np.float64(0.1)

#Peclet number
Pe = (rho*u*h) / gamma
print("Peclet number is", Pe)

#Iteration number
iterations = 0

#Initialising the temperature field
T = np.zeros(N)
T[0] = 100.
T[N-1] = 500.

#Initialising the new temp field
T_new = np.zeros(N)
T_new[0] = 100.
T_new[N-1] = 500.

numerical_error = np.float64(1.0)
epsilon = 1.E-08

plt.figure(10)

while numerical_error > epsilon:
    for i in range(1,N-1):
        a_W = gamma/h + max(rho*u,0)
        a_E = gamma/h - max(0,-rho*u)
        a_P = a_W + a_E
        T_new[i] = (a_E*T[i+1] + a_W*T[i-1])/ a_P
        
        #Resetting the numerical error
    numerical_error = 0
    for i in range(1,N-1):
        numerical_error = numerical_error + abs(T[i] - T_new[i])
        
    iterations = iterations + 1
    T = T_new.copy()
    
    if iterations%500 == 0:
        plt.figure(10)
        plt.semilogy(iterations, numerical_error,'ko')
        plt.pause(0.01)
        
x_dom = np.arange(N) * h

plt.figure(11)
plt.plot(x_dom, T, 'gx--', linewidth = 2, markersize = 10)

plt.grid(True, color = 'k')

plt.xlabel("Position", size = 20)
plt.ylabel("Temperature", size = 20)
plt.title("T(x)")

plt.show()        
        

        
        
        
        
        
    


