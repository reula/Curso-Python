# -*- coding: utf-8 -*-
#Caja con lados a diferentes potenciales.

#Queremos resolver el potencial de una caja cuyas paredes son conductoras (pero aisladas entre si) y pueden admitir potenciales arbitrarios (lógicamente constantes en cada una de ellas). Por la posibilidad de superponer soluciones es suficiente considerar que todas las caras, excepto una están a potencial cero, la restante a potencial 1. La solución general será la combinación lineal arbitraria de soluciones con distintas caras a potencial unidad.
#Consideraremos el caso donde el potencial distinto de cero es en la cara z
# de la derecha. Consideraremos que la caja tiene dimensiones (a,b,c).

#Primero cargamos las librerías necesarias de Python


from sympy import *
from sympy import sin, cos, pi, sinh
import numpy as np
import matplotlib.pyplot as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from sympy.abc import x, y, z
a, b, c= symbols('a b c')
k, m, n = symbols('k m n', integer=True)
potential = symbols('potential', cls=Function)

#Definimos las dimensiones de la caja:

a=1
b=1
c=1

#definimos el autovalor

gamma_z=lambda n,m : pi*sqrt(n**2/a**2 + m**2/b**2)

#calculamos los coeficientes

A=lambda n,m: 4*Integral(sin(pi*x*n/a),(x,0,a)) * Integral(sin(pi*y*m/b),(y,0,b)) / sinh(gamma_z(n,m)*c)

# hacemos L sumandos en la aproximación al potencial

L=50

potential = lambdify((x,y,z),Sum(
    Sum(A(m,k).doit()*sin(pi*m*x/a)*sin(pi*k*y/b)*sinh(gamma_z(m,k)*z),(m,1,L)),(k,1,L))/a/b,"numpy")

# Constatamos que la condición de contorno esté bien:
    
x_vals = np.linspace(0, a, 200)
#y_vals = linspace(0, b, 200)
#X,Y = meshgrid(x_vals,y_vals)
z_vals = potential(x_vals,b/2,c)
mpl.plot(x_vals, z_vals)
mpl.ylabel("at boundary")
mpl.show()

# Ahora graficamos el potencial completo (en y=b/2)

z_vals = np.linspace(0, c, 200)
X,Y = np.meshgrid(x_vals,z_vals)
p_vals = potential(X,b/2,Y)
fig = mpl.figure()
ax = fig.gca(projection='3d')
ax.set_zlim(-0.2, 1.2)
surf = ax.plot_surface(X,Y,p_vals, rstride=2, alpha=0.9, cstride=2, cmap=cm.coolwarm,linewidth=0, antialiased=False)
levels = np.arange(0., 1., 0.2)
levels_y = np.arange(0.9, 2., 0.1)
cset = ax.contour(X, Y, p_vals, levels, zdir='z', offset=1.5, cmap=cm.coolwarm)
cset = ax.contour(X, Y, p_vals, levels, zdir='x', offset=1.4, cmap=cm.coolwarm)
cset = ax.contour(X, Y, p_vals, levels_y, zdir='y', offset=-2, cmap=cm.coolwarm)
fig.colorbar(surf, shrink=0.5, aspect=5)
mpl.show()
