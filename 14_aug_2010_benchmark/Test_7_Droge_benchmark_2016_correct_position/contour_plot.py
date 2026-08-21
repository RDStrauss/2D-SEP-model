import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import random
import matplotlib.tri as tri
import matplotlib
from matplotlib.colors import LogNorm

#------------------------------------------------------------------------
# Load grids and data

R = np.loadtxt("r_grid.txt")
PHI = np.loadtxt("phi_grid.txt")

N = int(len(R))
M = int(len(PHI))

z = np.zeros(N*M,float)
x = np.zeros(N*M,float)
y = np.zeros(N*M,float)

k = 0
for i in range(0,N):
  for j in range(0,M):
    y[k] = R[i]*np.sin(PHI[j])
    x[k] = R[i]*np.cos(PHI[j])
    k = k + 1

z = np.loadtxt("output.txt")

print('Max intensity: ', np.amax(z))
z = z/np.amax(z)*100. #Normalize to maximim 100%

#------------------------------------------------------------------------
# Load position of spacecraft
data = np.loadtxt("Earth_omni_time.txt")
Earth_r = np.array([data[0,0],data[0,0]])
Earth_phi = np.array([data[0,1],data[0,1]])

data = np.loadtxt("STA_omni_time.txt")
STA_r = np.array([data[0,0],data[0,0]])
STA_phi = np.array([data[0,1],data[0,1]])

data = np.loadtxt("STB_omni_time.txt")
STB_r = np.array([data[0,0],data[0,0]])
STB_phi = np.array([data[0,1],data[0,1]])

data = np.loadtxt("Boundary_omni_time.txt")
Bounday_phi = data[0,1]

#------------------------------------------------------------------------
# Convert to Cartesian grid

ngridx = 500
ngridy = 500
xi = np.linspace(-3.0, 3.0, ngridx)
yi = np.linspace(-3.0, 3.0, ngridy)

triang = tri.Triangulation(x, y)
interpolator = tri.LinearTriInterpolator(triang, z)
xi, yi = np.meshgrid(xi, yi)
zi = interpolator(xi, yi)

#------------------------------------------------------------------------
# Variables for magnetic fieldlines later on

V_sw = 371./1.5e8*60.*60 #units of vsw in AU/hr
Omega = 2.*np.pi/25.38/24. #units of /hr

phi_plot = np.linspace(0, 2*np.pi, 100)

r = np.linspace(0.05, 3, 200)
phi = -(r - 0.05)/V_sw*Omega

#------------------------------------------------------------------------
# Do the plotting

fig = plt.figure(figsize = (10,10))
plt.title("t = 10 hrs", fontsize = 16)
plt.xlim(-3.1, 3.1)
plt.ylim(-3.1, 3.1)
plt.axes().set_aspect('equal')
plt.xlabel('x (AU)', fontsize = 16)
plt.ylabel('y (AU)', fontsize = 16)
plt.tick_params(labelsize = 14)

aap = plt.pcolormesh(xi, yi, zi, norm=LogNorm(vmax=10.0, vmin=0.001), cmap=plt.cm.jet)

plt.colorbar()

#------------------------------------------------------------------------
# Add some circles

plt.plot(3.*np.cos(phi_plot),3*np.sin(phi_plot), color = 'black',linewidth=15.0)
plt.plot(0.05*np.cos(phi_plot),0.05*np.sin(phi_plot), color = 'white',linewidth=10.0)

plt.plot(1.*np.cos(phi_plot),1.*np.sin(phi_plot), color = 'white',linewidth=2.0, linestyle = ':')
#plt.plot(0.5*np.cos(phi_plot),0.5*np.sin(phi_plot), color = 'white',linewidth=2.0, linestyle = ':')
#plt.plot(1.5*np.cos(phi_plot),1.5*np.sin(phi_plot), color = 'white',linewidth=2.0, linestyle = ':')
#plt.plot(2.0*np.cos(phi_plot),2.0*np.sin(phi_plot), color = 'white',linewidth=2.0)
#plt.plot(2.5*np.cos(phi_plot),2.5*np.sin(phi_plot), color = 'white',linewidth=2.0)

#------------------------------------------------------------------------
# Add the spacecraft

# Earth observer
Earth_x = Earth_r*np.cos(Earth_phi)
Earth_y = Earth_r*np.sin(Earth_phi)
plt.scatter(Earth_x, Earth_y, facecolor = 'green', marker = 's', s = 50, edgecolor = 'white')
phi_HMF = -(Earth_r[0] - 0.05)/V_sw*Omega
dphi = Earth_phi[0] - phi_HMF
plt.plot(r*np.cos(phi + dphi), r*np.sin(phi + dphi), color = 'green', linestyle = '-', linewidth = 2, label = 'Earth')

# STA observer
STA_x = STA_r*np.cos(STA_phi)
STA_y = STA_r*np.sin(STA_phi)
plt.scatter(STA_x, STA_y, facecolor = 'red', marker = 's', s = 50, edgecolor = 'white')
phi_HMF = -(STA_r[0] - 0.05)/V_sw*Omega
dphi = STA_phi[0] - phi_HMF
plt.plot(r*np.cos(phi + dphi), r*np.sin(phi + dphi), color = 'red', linestyle = '-', linewidth = 2, label = 'STA')

# STB observer
STB_x = STB_r*np.cos(STB_phi)
STB_y = STB_r*np.sin(STB_phi)
plt.scatter(STB_x, STB_y, facecolor = 'blue', marker = 's', s = 50, edgecolor = 'white')
phi_HMF = -(STB_r[0] - 0.05)/V_sw*Omega
dphi = STB_phi[0] - phi_HMF
plt.plot(r*np.cos(phi + dphi), r*np.sin(phi + dphi), color = 'blue', linestyle = '-', linewidth = 2, label = 'STB')

r = np.linspace(0.05, 3, 200)
phi = -(r - 0.05)/V_sw*Omega + Bounday_phi
plt.plot(r*np.cos(phi), r*np.sin(phi), color = 'black', linestyle = '--', linewidth = 2, label = 'Source')

plt.legend(ncol = 4, fontsize = 10)

fig.savefig('Intensity_contour_electrons.png', dpi=300, bbox_inches='tight')
#------------------------------------------------------------------------

plt.show()
  
