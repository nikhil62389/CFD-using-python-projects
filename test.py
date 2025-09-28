 # SImulation of fluid flow

import numpy as np
import matplotlib.pyplot as plt

def plot_results(poly, poly_der, x_0, h):
    #Visualisation the polynomial, its true derivative and numerical approximation
    
    #create a range of x values for a smoooth plot
    x = np.linspace(x_0 - 2, x_0 + 2, 400) # linspace to draw a smooth curve, we need a lot of points. This creates an array x with 400 evenly spaced points, starting from 2 units before x_0 and ending 2 units after x_0.
    y = np.polyval(poly, x)
    
    #Calculate the tangent line for true derivative at x_0
    tangent_y = np.polyval(poly_der, x_0)*(x - x_0) + np.polyval(poly,x_0)
    
    #Plotting
    plt.figure(figsize=(12, 8))
    plt.plot(x, y, label = 'Polynomial p(x)', color = 'blue')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    #analyse_derivatives
def analyse_derivative(poly_coeffs, x_0, h):
    # Calculates and compares the nuumerical derivatives and analytical derivatives of the polynomial.
    
    #1.Analytical (True) derivatives
   p_der_coeffs = np.polyder(poly_coeffs)
   p_der_analytical = np.polyval(p_der_coeffs, x_0)
   
   #2.Numerical derivatives
   p_x0 = np.polyval(poly_coeffs, x_0)
   p_x0_plus_h = np.polyval(poly_coeffs, x_0 + h)
   p_x0_minus_h = np.polyval(poly_coeffs, x_0 - h)
   
   forward_diff = (p_x0_plus_h - p_x0) / h
   backward_diff = (p_x0_minus_h - p_x0) / h
   central_diff = (p_x0_plus_h - p_x0_minus_h) / (2*h)
   
   #3.Error analysis
   error_forward = abs(forward_diff - p_der_analytical)
   error_backward = abs(backward_diff - p_der_analytical)
   error_central = abs(central_diff - p_der_analytical)
   
   #4.Visualisation
   plot_results(poly_coeffs, p_der_coeffs, x_0, h)
       
   # - -Main Execution - -    
       
 if __name__ == "__main__":
 
 #actual array of the polynomial is created 
     poly_p = np.array([-4, 7, -3, 9])
 
 x_val = 0.0
 h_val = 0.25
 
 analyse_derivative(poly_p, x_val, h_val)
 
 analyse_derivative(poly_p, x_val, 0.01)
 print("test")      
    
    
    

