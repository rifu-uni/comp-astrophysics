#!/usr/bin/env python
# coding: utf-8

# # Ordinary differential equations

# ## Analytical solution

# Solve the following ODE:
# 
# $$
# y'' + 3y' + \frac{25}{16}y = \cos x
# $$
# 
# with initial conditons:
# 
# $$
# y(0) = 5, \quad y'(0) = 3
# $$

# ### Step-by-step solution
# 
# **1. Homogeneous solution**
# 
# First, we find the roots of the characteristic equation:
# $$
# r^2 + 3r + \frac{25}{16} = 0
# $$
# 
# Using the quadratic formula:
# $$
# r = \frac{-3 \pm \sqrt{9 - 4 \cdot \frac{25}{16}}}{2} = \frac{-3 \pm \sqrt{9 - \frac{25}{4}}}{2} = \frac{-3 \pm \sqrt{\frac{11}{4}}}{2} = -\frac{3}{2} \pm \frac{\sqrt{11}}{4}
# $$
# 
# Thus, the homogeneous solution is:
# $$
# y_h(x) = c_1 e^{\left(-\frac{3}{2} + \frac{\sqrt{11}}{4}\right)x} + c_2 e^{\left(-\frac{3}{2} - \frac{\sqrt{11}}{4}\right)x}
# $$
# 
# **2. Particular solution**
# 
# We guess a particular solution of the form $y_p(x) = A \cos x + B \sin x$.
# Its derivatives are:
# $$
# y_p' = -A \sin x + B \cos x \\
# y_p'' = -A \cos x - B \sin x
# $$
# 
# Substituting into the ODE:
# $$
# (-A \cos x - B \sin x) + 3(-A \sin x + B \cos x) + \frac{25}{16}(A \cos x + B \sin x) = \cos x
# $$
# 
# Grouping sine and cosine terms:
# $$
# \left( \frac{9}{16}A + 3B \right) \cos x + \left( -3A + \frac{9}{16}B \right) \sin x = \cos x
# $$
# 
# This gives a system of equations:
# $$
# \begin{cases} 
# \frac{9}{16}A + 3B = 1 \\
# -3A + \frac{9}{16}B = 0 
# \end{cases}
# $$
# 
# From the second equation, $B = \frac{16}{3}A$. Substituting into the first:
# $$
# \frac{9}{16}A + 3\left(\frac{16}{3}A\right) = 1 \implies \frac{9}{16}A + 16A = 1 \implies \frac{265}{16}A = 1
# $$
# So $A = \frac{16}{265}$ and $B = \frac{16}{3} \cdot \frac{16}{265} = \frac{256}{795}$.
# 
# The particular solution is:
# $$
# y_p(x) = \frac{16}{265} \cos x + \frac{256}{795} \sin x
# $$
# 
# **3. General solution**
# 
# The general solution is $y(x) = y_h(x) + y_p(x)$:
# $$
# y(x) = c_1 e^{\left(-\frac{3}{2} + \frac{\sqrt{11}}{4}\right)x} + c_2 e^{\left(-\frac{3}{2} - \frac{\sqrt{11}}{4}\right)x} + \frac{16}{265} \cos x + \frac{256}{795} \sin x
# $$
# 
# **4. Applying initial conditions**
# 
# Let $\lambda_1 = -\frac{3}{2} + \frac{\sqrt{11}}{4}$ and $\lambda_2 = -\frac{3}{2} - \frac{\sqrt{11}}{4}$.
# 
# Using $y(0) = 5$:
# $$
# c_1 + c_2 + \frac{16}{265} = 5 \implies c_1 + c_2 = \frac{1309}{265} = S
# $$
# 
# Now differentiate $y(x)$ to use $y'(0) = 3$:
# $$
# y'(x) = c_1 \lambda_1 e^{\lambda_1 x} + c_2 \lambda_2 e^{\lambda_2 x} - \frac{16}{265} \sin x + \frac{256}{795} \cos x
# $$
# $$
# y'(0) = c_1 \lambda_1 + c_2 \lambda_2 + \frac{256}{795} = 3 \implies c_1 \lambda_1 + c_2 \lambda_2 = \frac{2129}{795} = T
# $$
# 
# Now we have the system:
# $$
# c_1 + c_2 = S \\
# c_1 \lambda_1 + c_2 \lambda_2 = T
# $$
# Using $c_2 = S - c_1$ in the second equation:
# $$
# c_1 \lambda_1 + (S - c_1) \lambda_2 = T \implies c_1(\lambda_1 - \lambda_2) + S \lambda_2 = T
# $$
# $$
# c_1 = \frac{T - S\lambda_2}{\lambda_1 - \lambda_2}
# $$
# $$
# c_2 = S - c_1
# $$
# 
# This completes the analytical solution.

# ## The Lane-Emden Equation (Euler & RK4)
# 
# I wanted a cool-looking equation even if I don't quite grasp the meaning because... because cool.
# 
# The Lane-Emden equation describes the structure of a self-gravitating sphere of polytropic fluid:
# $$
# \frac{1}{\xi^2} \frac{d}{d\xi} \left( \xi^2 \frac{d\theta}{d\xi} \right) = -\theta^n
# $$
# Expanding the derivative yields:
# $$
# \theta'' + \frac{2}{\xi} \theta' + \theta^n = 0
# $$
# 
# We can write this as a system of first-order ODEs:
# $$
# y_0 = \theta \\
# y_1 = \theta'
# $$
# $$
# y_0' = y_1 \\
# y_1' = -\theta^n - \frac{2}{\xi} y_1
# $$
# 
# Standard initial conditions at the center ($\xi = 0$) are $\theta(0) = 1$ and $\theta'(0) = 0$. To avoid the singularity at $\xi=0$, we start integration at a small $\xi_0$ using a Taylor expansion:
# $$\theta(\xi_0) \approx 1 - \frac{\xi_0^2}{6}, \quad \theta'(\xi_0) \approx -\frac{\xi_0}{3}$$
# 
# For $n=1$, the analytical solution is $\theta(\xi) = \frac{\sin(\xi)}{\xi}$.

# In[1]:


import sys
sys.path.append('../../')
import src.integrators as ode
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


# Parameters
n = 1
xi_0 = 1e-4
xi_f = 4.0
h = 0.05

# Initial conditions at xi_0
theta_0 = 1.0 - (xi_0**2) / 6.0
theta_prime_0 = -xi_0 / 3.0
y0 = np.array([theta_0, theta_prime_0])

def lane_emden_ode(xi: float, y: np.ndarray) -> np.ndarray:
    theta = y[0]
    theta_prime = y[1]
    # To handle negative theta for non-integer n, we can use np.abs or np.nan
    # For n=1, theta^n is just theta
    d_theta = theta_prime
    d_theta_prime = -(theta**n) - (2.0 / xi) * theta_prime
    return np.array([d_theta, d_theta_prime])

# Analytical solution for n=1
def le_exact_n1(xi): 
    return np.sin(xi) / xi

# Integrate using Euler
xi_euler, y_euler = ode.integrate(lane_emden_ode, (xi_0, xi_f), y0, h, ode.euler_step)
theta_euler = y_euler[:, 0]

# Integrate using RK4 (Runge)
xi_rk4, y_rk4 = ode.integrate(lane_emden_ode, (xi_0, xi_f), y0, h, ode.rk4_runge_step)
theta_rk4 = y_rk4[:, 0]


# In[3]:


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Main plot
xi_exact = xi_rk4
theta_exact = le_exact_n1(xi_exact)
ax1.plot(xi_exact, theta_exact, 'k--', label='Analytical (n=1)', lw=2)
ax1.plot(xi_euler, theta_euler, label='Euler')
ax1.plot(xi_rk4, theta_rk4, label='RK4 (Runge)')
ax1.set_xlabel(r'$\xi$')
ax1.set_ylabel(r'$\theta(\xi)$')
ax1.set_title('Lane-Emden Equation (n=1)')
ax1.legend()
ax1.grid(True)

# Error plot
error_euler = np.abs(theta_euler - le_exact_n1(xi_euler))
error_rk4 = np.abs(theta_rk4 - le_exact_n1(xi_rk4))

ax2.plot(xi_euler[1:], error_euler[1:], label='Euler Error')
ax2.plot(xi_rk4[1:], error_rk4[1:], label='RK4 Error')
ax2.set_yscale('log')
ax2.set_xlabel(r'$\xi$')
ax2.set_ylabel('Absolute Error')
ax2.set_title('Absolute Error Comparison')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()

