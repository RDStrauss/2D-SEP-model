import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------------------------------
data = np.loadtxt("./1D_SEP_model/output.txt")

time = data[:,2]
omni_directional_intensity = data[:,3]
anisotropy = data[:,4]
z = data[0,0]
r = data[0,1]

#------------------------------------------------------------------------
#Load data
data = np.loadtxt("./1D_SEP_model/Ani.csv", skiprows = 1, delimiter = ',')
time_ani = data[:,0]
ani = data[:,1]


data = np.loadtxt("./1D_SEP_model/DI.csv", skiprows = 1, delimiter = ',')
time_di = data[:,0]
di = data[:,1]

data = np.loadtxt("Earth_omni_time.txt")

time_Earth = data[:,2]
f_omni_Earth = data[:,3]
anis_Earth = data[:,4]

#normalize model di to observations
omni_directional_intensity = omni_directional_intensity/np.amax(omni_directional_intensity)*np.amax(di)
f_omni_Earth = f_omni_Earth/np.amax(f_omni_Earth)*np.amax(di)

#Plot the figures
fig = plt.figure(figsize=(5,6))

subplot2 = fig.add_subplot(211)
#subplot2.set_xlabel('Time (hr)')
#subplot2.set_title('aap')
subplot2.set_ylabel('Intensity (arb. units)')
subplot2.set_yscale('log')
subplot2.set_xscale('linear')
subplot2.set_ylim([0.005,3.])
#subplot2.ticklabel_format(useOffset=False)
subplot2.set_xlim([-0.5,12])


subplot2.plot(time,omni_directional_intensity, color = 'green', linestyle = '--', label = '1D finite difference')
subplot2.plot(time_Earth, f_omni_Earth, color = 'red', label = '2D finite difference')
subplot2.scatter(time_di, di, color = 'blue', facecolor = 'none', label = 'SDE results (Droge et al., 2010)')
subplot2.set_xticklabels([])

subplot2.legend()

subplot3 = fig.add_subplot(413)  
subplot3.set_xlabel('Time')
subplot3.set_ylabel('Anisotropy')
subplot3.set_yscale('linear')
subplot3.set_xscale('linear')
subplot3.set_ylim([-0.2,3.0])
#subplot3.ticklabel_format(useOffset=False)
subplot3.set_xlim([-0.5,12])
  
subplot3.plot(time, anisotropy, color = 'green', linestyle = '--')
subplot3.plot(time_Earth, anis_Earth, color = 'red')
subplot3.scatter(time_ani, ani, color = 'blue', facecolor = 'none')
subplot3.axhline(y = 0, color = 'black', linestyle = '--')

plt.subplots_adjust(hspace=0.15)
fig.align_ylabels([subplot2,subplot3])
fig.tight_layout

fig.savefig('Output.png', dpi=300, bbox_inches='tight')

#--------------------------------------------------------------------

plt.show()
