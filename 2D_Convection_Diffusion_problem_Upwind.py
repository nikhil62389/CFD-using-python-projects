import numpy as np
import matplotlib.pyplot as plt

N = 21

L = 1

h = np.float64(L/(N - 1))
iterations = 0

gamma = np.float64(5.0);
rho = np.float64(1.25);

u = np.float64(2.0);
v= np.float64(2.0);

Pe = u*h/gamma
print("Peclet number is", Pe)

T = np.zeros((N,N))
T[0,:] = 1.
T[:,0] = 1.

T_new = np.zeros((N,N))
T_new[0,:] = 1.
T_new[:,0] = 1.

epsilon = 1E-08
numerical_error = 1

plt.figure(10)

while numerical_error > epsilon:
    for i in range (1,N-1):
        for j in range(1,N-1):
            a_E = gamma/h + max(-rho*u,0)
            a_W = gamma/h + max(rho*u,0)
            a_N = gamma/h + max(-rho*v,0)
            a_S = gamma/h + max(rho*u,0)
            a_P = a_E + a_W + a_N + a_S
            T_new[i,j] = (a_E*T[i,j+1]+a_W*T[i,j-1]+a_N*T[i-1,j]+a_S*T[i+1,j]) / a_P
                          
    numerical_error = 0
    for i in range(1,N-1):
        for j in range(1,N-1):
            numerical_error = numerical_error + abs(T[i,j] - T_new[i,j])
    iterations = iterations + 1
    T = T_new.copy()
    
    if iterations%250 ==0:
        plt.figure(10)
        plt.semilogy(iterations,numerical_error, 'ko')
        plt.pause(0.01)
        
x_dom = np.arange(N) * h
y_dom = L - np.arange(N)*h
[X,Y] = np.meshgrid(x_dom,y_dom)
            
plt.figure(11)
plt.contourf(X,Y,T,12)

plt.colorbar(orientation = 'vertical')
plt.title("T(x,y)")
plt.show()
       