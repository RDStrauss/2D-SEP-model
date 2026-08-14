import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import random
import matplotlib.tri as tri
import matplotlib
from matplotlib.colors import LogNorm

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


print(k)

z = np.loadtxt("output.txt")

print('Max intensity')
print(np.amax(z))

z = z/np.amax(z)*100.
#z = np.log(z)

ngridx = 500
ngridy = 500
xi = np.linspace(-3.0, 3.0, ngridx)
yi = np.linspace(-3.0, 3.0, ngridy)

triang = tri.Triangulation(x, y)
interpolator = tri.LinearTriInterpolator(triang, z)
xi, yi = np.meshgrid(xi, yi)
zi = interpolator(xi, yi)

#zi = interpolate.griddata(values, z, xi,yi, method='linear', fill_value=np.nan, rescale=False)

#f = interpolate.interp2d(x, y, z, kind='linear')
#zi = f(xi, yi)
#zi = mlab.griddata(x, y, z, xi, yi, interp='linear')


V_sw = 361./1.5e8*60.*60 #units of vsw in AU/hr
Omega = 2.*np.pi/25.4/24. #units of /hr

phi_plot = np.linspace(0, 2*np.pi, 100)

r = np.linspace(0, 3, 200)
phi = -r/V_sw*Omega



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
#plt.contourf(xi, yi, np.log(zi), 255, cmap=plt.cm.jet)

#plt.tricontour(x, y, z, levels=14, linewidths=0.5, colors='k')
#cntr2 = plt.tricontourf(x, y, z, levels=255, cmap="gist_heat", vmin = 90., vmax = 100.)

#plt.plot(1.*np.cos(phi_plot),1.*np.sin(phi_plot), color = 'white',linewidth=3.0, linestyle = '-')

plt.plot(3.*np.cos(phi_plot),3*np.sin(phi_plot), color = 'black',linewidth=15.0)
plt.plot(0.05*np.cos(phi_plot),0.05*np.sin(phi_plot), color = 'white',linewidth=10.0)

plt.plot(1.*np.cos(phi_plot),1.*np.sin(phi_plot), color = 'white',linewidth=2.0, linestyle = ':')
plt.plot(0.5*np.cos(phi_plot),0.5*np.sin(phi_plot), color = 'white',linewidth=2.0, linestyle = ':')
plt.plot(1.5*np.cos(phi_plot),1.5*np.sin(phi_plot), color = 'white',linewidth=2.0, linestyle = ':')
#plt.plot(2.0*np.cos(phi_plot),2.0*np.sin(phi_plot), color = 'white',linewidth=2.0)
#plt.plot(2.5*np.cos(phi_plot),2.5*np.sin(phi_plot), color = 'white',linewidth=2.0)

# Earth observer
Earth_r = np.array([0.993490032,0.993490032])
Earth_phi = np.array([2.35E-06,2.35E-06])
Earth_phi = Earth_phi/180.*np.pi
Earth_x = Earth_r*np.cos(Earth_phi)
Earth_y = Earth_r*np.sin(Earth_phi)
plt.scatter(Earth_x, Earth_y, facecolor = 'green', marker = 's', s = 50, edgecolor = 'white')
dphi = 82.36286947
plt.plot(r*np.cos(phi + dphi), r*np.sin(phi + dphi), color = 'green', linestyle = '-', linewidth = 2)

# Mars observer
Mars_r = np.array([1.610973075,1.610973075])
Mars_phi = np.array([169.0638746,169.0638746])
Mars_phi = Mars_phi/180.*np.pi
Mars_x = Mars_r*np.cos(Mars_phi)
Mars_y = Mars_r*np.sin(Mars_phi)
plt.scatter(Mars_x, Mars_y, facecolor = 'maroon', marker = 's', s = 50, edgecolor = 'white')
dphi = 271.9469161/180.*np.pi
plt.plot(r*np.cos(phi + dphi), r*np.sin(phi + dphi), color = 'maroon', linestyle = '-', linewidth = 2)

# Bepi observer
Bepi_r = np.array([0.411289483,0.411289483])
Bepi_phi = np.array([89.94839563,89.94839563])
Bepi_phi = Bepi_phi/180.*np.pi
Bepi_x = Bepi_r*np.cos(Bepi_phi)
Bepi_y = Bepi_r*np.sin(Bepi_phi)
plt.scatter(Bepi_x, Bepi_y, facecolor = 'orange', marker = 's', s = 50, edgecolor = 'white')
dphi = 112.6313914/180.*np.pi
plt.plot(r*np.cos(phi + dphi), r*np.sin(phi + dphi), color = 'orange', linestyle = '-', linewidth = 2)

# SolO observer
SolO_r = np.array([0.803862906,0.803862906])
SolO_phi = np.array([-2.819722862,-2.819722862])
SolO_phi = SolO_phi/180.*np.pi
SolO_x = SolO_r*np.cos(SolO_phi)
SolO_y = SolO_r*np.sin(SolO_phi)
plt.scatter(SolO_x, SolO_y, facecolor = 'blue', marker = 's', s = 50, edgecolor = 'white')
dphi = 55.79048569/180.*np.pi
plt.plot(r*np.cos(phi + dphi), r*np.sin(phi + dphi), color = 'blue', linestyle = '-', linewidth = 2)

# STA observer
STA_r = np.array([0.958293441,0.958293441])
STA_phi = np.array([-37.59619789,-37.59619789])
STA_phi = STA_phi/180.*np.pi 
STA_x = STA_r*np.cos(STA_phi)
STA_y = STA_r*np.sin(STA_phi)
plt.scatter(STA_x, STA_y, facecolor = 'red', marker = 's', s = 50, edgecolor = 'white')
dphi = 28.0044515/180.*np.pi
plt.plot(r*np.cos(phi + dphi), r*np.sin(phi + dphi), color = 'red', linestyle = '-', linewidth = 2)

# PSP observer
PSP_r = np.array([0.623189194,0.623189194])
PSP_phi = np.array([-53.64955116,-53.64955116])
PSP_phi = PSP_phi/180.*np.pi 
PSP_x = PSP_r*np.cos(PSP_phi)
PSP_y = PSP_r*np.sin(PSP_phi)
plt.scatter(PSP_x, PSP_y, facecolor = 'purple', marker = 's', s = 50, edgecolor = 'white')
dphi = -13.95688608/180.*np.pi
plt.plot(r*np.cos(phi + dphi), r*np.sin(phi + dphi), color = 'purple', linestyle = '-', linewidth = 2)


r = np.linspace(0, 3, 200)
phi = -r/V_sw*Omega 
plt.plot(r*np.cos(phi), r*np.sin(phi), color = 'black', linestyle = '--', linewidth = 2, label = 'Source')

#plt.legend(ncol = 4, fontsize = 10)

fig.savefig('Intensity_contour_electrons.png', dpi=300, bbox_inches='tight')
#------------------------------------------------------------------------

#------------------------------------------------------------------------

plt.show()
  
