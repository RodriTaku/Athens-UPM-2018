# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 10:58:40 2018

@author: RodrigoTakuro
"""

import matplotlib.pyplot as plt
my_sol = brachsim( 0., 0., 2., 1., 0.5, 1 )
plt.plot(my_sol[0,:], -my_sol[1,:])

A = 0.5171999217
theta_1 = 3.508368769

x = lambda A,theta : A*(theta - sin(theta))
y = lambda A,theta : A*(1 - cos(theta))
theta_vec = linspace(0, theta_1, 21)
t_1 = sqrt(A*1)*theta_1

plt.plot(x(A,theta_vec), -y(A,theta_vec))