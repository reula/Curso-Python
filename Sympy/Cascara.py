# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 12:52:43 2018

@author: reula
"""

from sympy import *
from sympy import legendre
from numpy import *
#from sympy.plotting import plot
import matplotlib.pyplot as mpl
#import matplotlib as mpl
#import matplotlib as mpl
#import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import animation, rc
from IPython.display import HTML
#init_printing(use_unicode=True)
x, y, z = symbols('x y z')
k, m, n = symbols('k m n', integer=True)
f, step, potential = symbols('f step potential', cls=Function)
#var('n m x')
N=200
A=lambda m: (2*(2*m+1)+1)*Integral(legendre(2*m+1,x),(x,0,1))
step=lambdify(x,Sum(A(m).doit()*legendre(2*m+1,x),(m,0,20)).doit(),"numpy")
x_vals = linspace(-1, 1, N)
z_vals = step(x_vals)
mpl.plot(x_vals, z_vals)
mpl.ylabel("step")
mpl.show()
potential=lambdify((x,y),Sum(A(m).doit()*legendre(2*m+1,x)*y**(-2*(m+1)),(m,0,20)).doit(),"numpy")
y_vals = linspace(1,10,N)
X,Y = meshgrid(x_vals,y_vals)
R,T = sqrt(X*X+Y*Y), X/sqrt(X*X+Y*Y)
#z_vals = potential(T,R)
z_vals = potential(X,Y)
fig = mpl.figure()
ax = fig.gca(projection='3d')
ax.set_zlim(-1.2, 1.2)
surf = ax.plot_surface(X,Y,z_vals, rstride=2, alpha=0.9, cstride=2, cmap=cm.coolwarm,linewidth=0, antialiased=False)
levels = arange(-1, 1, 0.2)
cset = ax.contour(X, Y, z_vals, levels, zdir='z', offset=1.5, cmap=cm.coolwarm)
cset = ax.contour(X, Y, z_vals, levels, zdir='x', offset=1.4, cmap=cm.coolwarm)
cset = ax.contour(X, Y, z_vals, levels, zdir='y', offset=-2, cmap=cm.coolwarm)
fig.colorbar(surf, shrink=0.5, aspect=5)
mpl.show()