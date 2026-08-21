import matplotlib.pyplot as plt
import numpy as np

#----------------------------------------------------------------------

fig = plt.figure(figsize=(10,15))

subplot1 = fig.add_subplot(211)

#subplot1.set_title('7 February 2010')
subplot1.set_xlabel('Time (hrs)')
subplot1.set_ylabel('Omni-directional intensity')
subplot1.set_yscale('log')
subplot1.set_xscale('linear')
#subplot1.set_ylim([0.001,2.])
#subplot1.ticklabel_format(useOffset=False)
subplot1.set_xlim([-2.5,20])

data = np.loadtxt("STB_omni_time.txt")
subplot1.plot(data[:,2],data[:,3], color = 'blue', linestyle = '-', linewidth = 1, label = 'STB')

data = np.loadtxt("Earth_omni_time.txt")
subplot1.plot(data[:,2],data[:,3], color = 'green', linestyle = '-', linewidth = 1, label = 'ACE')

data = np.loadtxt("STA_omni_time.txt")
subplot1.plot(data[:,2],data[:,3], color = 'red', linestyle = '-', linewidth = 1, label = 'STA')

data = np.loadtxt("Mars_omni_time.txt")
subplot1.plot(data[:,2],data[:,3], color = 'black', linestyle = '-', linewidth = 1, label = 'Mars')

data = np.loadtxt("Bepi_omni_time.txt")
subplot1.plot(data[:,2],data[:,3], color = 'orange', linestyle = '-', linewidth = 1, label = 'Bepi')

data = np.loadtxt("PSP_omni_time.txt")
subplot1.plot(data[:,2],data[:,3], color = 'purple', linestyle = '-', linewidth = 1, label = 'PSP')

data = np.loadtxt("SolO_omni_time.txt")
subplot1.plot(data[:,2],data[:,3], color = 'gray', linestyle = '-', linewidth = 1, label = 'SolO')

plt.legend(loc = 4)

subplot1 = fig.add_subplot(614)

#subplot1.set_title('Anisotropy vs time')
subplot1.set_xlabel('Time (hrs)')
subplot1.set_ylabel('Anisotropy')
subplot1.set_yscale('linear')
subplot1.set_xscale('linear')
subplot1.set_ylim([-0.25,2])
subplot1.set_xlim([-2.5,20])

data = np.loadtxt("STA_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'red', linestyle = '-', label = 'STEREO A', linewidth = 1)

data = np.loadtxt("STB_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'blue', linestyle = '-', label = 'STEREO B', linewidth = 1)

data = np.loadtxt("Earth_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'green', linestyle = '-', label = 'Earth', linewidth = 1)

subplot1.axhline(y = 0, color = 'black')

plt.legend(loc = 1)
#data = np.loadtxt("omni_time_best_connection.txt")
#subplot1.plot(data[:,2],data[:,4], color = 'black', label = 'Best connection')

subplot1 = fig.add_subplot(615)

#subplot1.set_title('Anisotropy vs time')
subplot1.set_xlabel('Time (hrs)')
subplot1.set_ylabel('Anisotropy')
subplot1.set_yscale('linear')
subplot1.set_xscale('linear')
subplot1.set_ylim([-0.25,2])
subplot1.set_xlim([-2.5,20])

data = np.loadtxt("Mars_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'black', linestyle = '-', label = 'Mars', linewidth = 1)

data = np.loadtxt("Bepi_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'orange', linestyle = '-', label = 'Bepi', linewidth = 1)


subplot1.axhline(y = 0, color = 'black')

plt.legend(loc = 1)

subplot1 = fig.add_subplot(616)

#subplot1.set_title('Anisotropy vs time')
subplot1.set_xlabel('Time (hrs)')
subplot1.set_ylabel('Anisotropy')
subplot1.set_yscale('linear')
subplot1.set_xscale('linear')
subplot1.set_ylim([-0.25,2])
subplot1.set_xlim([-2.5,20])

data = np.loadtxt("PSP_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'purple', linestyle = '-', label = 'PSP', linewidth = 1)

data = np.loadtxt("SolO_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'grey', linestyle = '-', label = 'SolO', linewidth = 1)


subplot1.axhline(y = 0, color = 'black')

plt.legend(loc = 1)

#plt.legend(loc = 1)

fig.savefig('output_2.png', dpi=300)

#------------------------------------------------------------------------





plt.show()
