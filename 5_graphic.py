#!/usr/bin/python3
from numpy import *
##################################################################################
from matplotlib.pyplot import *
from matplotlib.animation import FuncAnimation
from pylab import *
from math import *              # Libreria para usar exp()
from mpl_toolkits.mplot3d import Axes3D
##################################################################################
import warnings
warnings.filterwarnings("ignore")
##################################################################################
##################################################################################

datos1 = loadtxt('file1.dat',float)
datos2 = loadtxt('file2.dat',float)
datos3 = loadtxt('file3.dat',float)

aux1 = datos1[:]
aux2 = datos2[:]
aux3 = datos3[:,:]

radio=[]
for i in range(0,len(aux1),1):
    radio.append(aux1[i])
radio.append(1.0)

theta=[]
for i in range(0,len(aux2),1):
    theta.append(aux2[i])

temp=zeros((len(radio),len(theta)))
for i in range(0,len(radio)-1,1):
    temp[i,:]=aux3[i,:]

'''
print(radio)
print(theta)
print(temp)
'''

fig = figure()
ax = Axes3D(fig)
subplot(projection="polar")

th, r = meshgrid(theta,radio)

pcolormesh(th,r,temp,cmap='jet')#,shading='gouraud')
cbar=colorbar(orientation='horizontal')
cbar.ax.set_title(r'$Temperatura$')#\ (^{\circ}C)$')
axis('off')
savefig('solucion.eps')
savefig('solucion.png')
show()
