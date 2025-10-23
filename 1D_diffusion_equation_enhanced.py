import numpy as np
import matplotlib.pyplot as plt
import time

def solve_1d_heat_equation(N = 11, L=1.0, T_left = 0.0, T_right = 1.0, epsilon = 1e-08, max_iterations = 10000, method = 'gauss-seidal'):
    

