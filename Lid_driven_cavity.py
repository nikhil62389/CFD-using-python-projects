import numpy as np
import matplotlib.pyplot as plt

# %% Domain and Parameters
N = 51
L = 1.0
h = L / (N - 1)
Re = 100
nu = 1.0 / Re

alpha = 0.8
alpha_p = 0.8

# %% Initialize Variables (Fixed tuple syntax)
u = np.zeros((N + 1, N))
v = np.zeros((N, N + 1))
p = np.ones((N + 1, N + 1))

u_star = np.zeros_like(u)
v_star = np.zeros_like(v)
d_e = np.zeros_like(u)
d_n = np.zeros_like(v)
pc = np.zeros((N + 1, N + 1))
b = np.zeros((N + 1, N + 1))

# Initial Lid Velocity
u[0, :] = 2.0 

error = 1.0
iterations = 0
error_req = 1e-7

# %% SIMPLE Algorithm Loop
while error > error_req:
    # 1. x-momentum (u_star)
    for i in range(1, N):
        for j in range(1, N - 1):
            u_E = 0.5 * (u[i, j] + u[i, j + 1])
            u_W = 0.5 * (u[i, j] + u[i, j - 1])
            v_N = 0.5 * (v[i - 1, j] + v[i - 1, j + 1])
            v_S = 0.5 * (v[i, j] + v[i, j + 1])
            
            a_E = (-0.5 * u_E * h) + nu
            a_W = (0.5 * u_W * h) + nu
            a_N = (-0.5 * v_N * h) + nu
            a_S = (0.5 * v_S * h) + nu
            a_e = (0.5 * u_E * h) - (0.5 * u_W * h) + (0.5 * v_N * h) - (0.5 * v_S * h) + (4 * nu)
            
            d_e[i, j] = -h / a_e
            u_star[i, j] = (a_E * u[i, j + 1] + a_W * u[i, j - 1] + a_N * u[i - 1, j] + a_S * u[i + 1, j]) / a_e + \
                           d_e[i, j] * (p[i, j + 1] - p[i, j])
    
    # u_star Boundaries
    u_star[0, :] = 2.0 - u_star[1, :]
    u_star[N, :] = -u_star[N - 1, :]
    u_star[1:N, 0] = 0.0
    u_star[1:N, N - 1] = 0.0 
    
    # 2. y-momentum (v_star)
    for i in range(1, N - 1):
        for j in range(1, N):
            u_E = 0.5 * (u[i, j] + u[i + 1, j])
            u_W = 0.5 * (u[i, j - 1] + u[i + 1, j - 1])
            v_N = 0.5 * (v[i, j] + v[i - 1, j])
            v_S = 0.5 * (v[i, j] + v[i + 1, j])
            
            a_E = -0.5 * u_E * h + nu
            a_W = 0.5 * u_W * h + nu
            a_N = -0.5 * v_N * h + nu
            a_S = 0.5 * v_S * h + nu
            a_n = 0.5 * u_E * h - 0.5 * u_W * h + 0.5 * v_N * h - 0.5 * v_S * h + 4 * nu
   
            d_n[i, j] = -h / a_n
            v_star[i, j] = (a_E * v[i, j + 1] + a_W * v[i, j - 1] + a_N * v[i - 1, j] + a_S * v[i + 1, j]) / a_n + \
                           d_n[i, j] * (p[i, j] - p[i + 1, j])
    
    # v_star Boundaries
    v_star[:, 0] = -v_star[:, 1]        
    v_star[:, N] = -v_star[:, N - 1]      
    v_star[0, 1:N] = 0.0                
    v_star[N - 1, 1:N] = 0.0

    # 3. Pressure Correction (Poisson Solver)
    pc[:, :] = 0.0 
    for p_iter in range(100):
        for i in range(1, N):
            for j in range(1, N):
                a_E_pc = -d_e[i, j] * h
                a_W_pc = -d_e[i, j - 1] * h
                a_N_pc = -d_n[i - 1, j] * h
                a_S_pc = -d_n[i, j] * h
                a_P_pc = a_E_pc + a_W_pc + a_N_pc + a_S_pc
                
                # Continuity source term
                b[i, j] = -(u_star[i, j] - u_star[i, j - 1]) * h + (v_star[i, j] - v_star[i - 1, j]) * h
                pc[i, j] = (a_E_pc * pc[i, j + 1] + a_W_pc * pc[i, j - 1] + a_N_pc * pc[i - 1, j] + a_S_pc * pc[i + 1, j] + b[i, j]) / a_P_pc
        
        # pc BCs
        pc[0, :] = pc[1, :]
        pc[N, :] = pc[N - 1, :]
        pc[:, 0] = pc[:, 1]
        pc[:, N] = pc[:, N - 1]

    # 4. Correct Fields
    p = p + alpha_p * pc
    
    for i in range(1, N):
        for j in range(1, N - 1):
            u[i, j] = u_star[i, j] + alpha * d_e[i, j] * (pc[i, j + 1] - pc[i, j])

    for i in range(1, N - 1):
        for j in range(1, N):
            v[i, j] = v_star[i, j] + alpha * d_n[i, j] * (pc[i, j] - pc[i + 1, j])

    # Final Boundary Updates
    u[0, :] = 2.0 - u[1, :]
    u[N, :] = -u[N - 1, :]
    u[1:N, 0] = 0.0
    u[1:N, N - 1] = 0.0
    v[:, 0] = -v[:, 1]
    v[:, N] = -v[:, N - 1]
    v[0, 1:N] = 0.0
    v[N - 1, 1:N] = 0.0

    # Residual calculation
    error = np.sum(np.abs(b[1:N, 1:N]))
    
    if iterations % 500 == 0:
        print(f"Iteration: {iterations}, Error: {error:.4e}")
    
    iterations += 1

# %% Post-Processing
u_final = np.zeros((N, N))
v_final = np.zeros((N, N))

for i in range(N):
    for j in range(N):
        u_final[i, j] = 0.5 * (u[i, j] + u[i + 1, j])
        v_final[i, j] = 0.5 * (v[i, j] + v[i, j + 1])

# %% Visualization
x_dom = np.linspace(0, L, N)
y_dom = np.linspace(0, L, N) # Standard 0 to 1
X, Y = np.meshgrid(x_dom, y_dom)

# Flip u_final vertically so lid is at top (index 0 is y=1)
plt.figure(figsize=(7, 5))
cp = plt.contourf(X, Y, np.flipud(u_final), 21, cmap='jet')
plt.colorbar(cp)
plt.title('U-Velocity Contours')
plt.show()

plt.figure(figsize=(7, 5))
plt.quiver(X, Y, u_final, v_final, color='black', scale=15)
plt.gca().invert_yaxis() # Match CFD convention
plt.axis('equal')
plt.show()