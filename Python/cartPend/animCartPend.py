# -*- coding: utf-8 -*-
"""
Created on March, 2018

@author: RodrigoTakuro
"""

from numpy import array, sin, cos, sqrt, pi, append, max, min
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Circle
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#If your simulation comes from cartPendOdeInt uncomment this:
X = odeintSol.transpose()
F = array([f(t) for t in T])
U = array([u(T[i],X[0,i],X[1,i],X[2,i],X[3,i]) for i in range(T.size)])

##If your simulation comes from cartPendSim uncomment this:
#X = Sol['X']
#T = Sol['T']
#F = append(append(Sol['F'][0,1],(Sol['F'][0,2:]+Sol['F'][1,1:-1])/2),Sol['F'][1,-1])
#U = append(append(Sol['U'][0,1],(Sol['U'][0,2:]+Sol['U'][1,1:-1])/2),Sol['U'][1,-1])

f_scaling = 2.
u_scaling = 2.
cartLength = 1./2
cartHeight = 1./4

codes = [Path.MOVETO] + [Path.LINETO]*3 + [Path.CLOSEPOLY]
sides = array([[-cartLength/2, -cartHeight/2], [cartLength/2, -cartHeight/2], [cartLength/2, cartHeight/2], [-cartLength/2, cartHeight/2], [-cartLength/2, -cartHeight/2]])
verts = sides.copy()

bobCenter = array([X[0,0] + param.L*sin(X[1,0]), -param.L*cos(X[1,0])])

def animate(i):
    global verts, bobCenter, X, F, U, param
    verts[:,0] = sides[:,0] + X[0,i]
    bobCenter = array([X[0,i] + param.L*sin(X[1,i]), -param.L*cos(X[1,i])])
    circle.center = bobCenter
    line.set_xdata([X[0,i],bobCenter[0]])
    line.set_ydata([0,bobCenter[1]])
    f_arrows.set_offsets(bobCenter)
    f_arrows.set_UVC(F[i]*f_scaling,0)
    u_arrows.set_offsets([X[0,i],0])
    u_arrows.set_UVC(U[i]*u_scaling,0)
    
    return [pathpatch, circle, line, f_arrows, u_arrows, ]

path = Path(verts, codes)
pathpatch = PathPatch(path, facecolor='b', edgecolor='b', animated=False)
circle = Circle(bobCenter, sqrt(cartLength*cartHeight*param.m/(pi*param.M)), facecolor='g', animated=False)

fig, ax = plt.subplots()
ax.add_patch(pathpatch)
ax.add_patch(circle)
line, = ax.plot([X[0,0],bobCenter[0]], [0,bobCenter[1]], linewidth=2, color='k')
f_arrows = ax.quiver(bobCenter[0], bobCenter[1], F[0]*f_scaling, 0, color='m', width=0.025, scale=1, units='xy')
u_arrows = ax.quiver(X[0,0], 0, U[0]*u_scaling, 0, color='r', width=0.025, scale=1, units='xy')

ax.set_xlim([min(X[0,:])-2*max([param.L,cartLength/2]),max(X[0,:])+2*max([param.L,cartLength/2])])
ax.set_ylim([2*min(append(-param.L*cos(X[1,:]),[-cartHeight/2,-param.L/2])),2*max(append(-param.L*cos(X[1,:]),[cartHeight/2,param.L/2]))])
ax.set_aspect('equal')

ani = animation.FuncAnimation(fig, animate, frames=X.shape[1],
                              interval=25, repeat=True, blit=True)
plt.show()