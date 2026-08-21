import matplotlib.pyplot as plt
import numpy as np


background = 0. # 0.00115 # data background
#----------------------------------------------------------------------

fig = plt.figure(figsize=(10,15))

subplot1 = fig.add_subplot(211)

subplot1.set_title('7 February 2010')
subplot1.set_xlabel('Time (hrs)')
subplot1.set_ylabel('Omni-directional intensity')
subplot1.set_yscale('log')
subplot1.set_xscale('linear')
subplot1.set_ylim([0.001,2.])
#subplot1.ticklabel_format(useOffset=False)
subplot1.set_xlim([-2.5,20])

data = np.loadtxt("./Test_4_Droge_benchmark_correct_position/STB_omni_time.txt")
aap = np.amax(data[:,3])
data[:,3] = data[:,3]/aap + background
subplot1.plot(data[:,2],data[:,3], color = 'blue', linestyle = '-', linewidth = 1)

data = np.loadtxt("./Test_4_Droge_benchmark_correct_position/Earth_omni_time.txt")
data[:,3] = data[:,3]/aap + background
subplot1.plot(data[:,2],data[:,3], color = 'green', linestyle = '-', linewidth = 1)

data = np.loadtxt("./Test_4_Droge_benchmark_correct_position/STA_omni_time.txt")
data[:,3] = data[:,3]/aap + background
subplot1.plot(data[:,2],data[:,3], color = 'red', linestyle = '-', linewidth = 1, label = '2D model solutions; correct position')

data = np.loadtxt("./Test_5_Droge_benchmark_footpoint_position/STB_omni_time.txt")
aap = np.amax(data[:,3])
data[:,3] = data[:,3]/aap + background
subplot1.plot(data[:,2],data[:,3], color = 'blue', linestyle = '--', linewidth = 1)

data = np.loadtxt("./Test_5_Droge_benchmark_footpoint_position/Earth_omni_time.txt")
data[:,3] = data[:,3]/aap + background
subplot1.plot(data[:,2],data[:,3], color = 'green', linestyle = '--', linewidth = 1)

data = np.loadtxt("./Test_5_Droge_benchmark_footpoint_position/STA_omni_time.txt")
data[:,3] = data[:,3]/aap + background
subplot1.plot(data[:,2],data[:,3], color = 'red', linestyle = '--', linewidth = 1, label = '2D model solutions; footpoint position')


data = np.loadtxt("./data/intensity_stereo_b.csv", delimiter = ',')
aap = np.amax(data[:,1])
data[:,1] = data[:,1]/aap
subplot1.plot(data[:,0] - 2.5,data[:,1], color = 'blue', linestyle = '-', alpha = 0.5, linewidth = 5)


data = np.loadtxt("./data/intensity_earth.csv", delimiter = ',')
data[:,1] = data[:,1]/aap
subplot1.plot(data[:,0] - 2.5,data[:,1], color = 'green', linestyle = '-', alpha = 0.5, linewidth = 5)


data = np.loadtxt("./data/intensity_stereo_a.csv", delimiter = ',')
data[:,1] = data[:,1]/aap
subplot1.plot(data[:,0] - 2.5,data[:,1], color = 'red', linestyle = '-', alpha = 0.5, linewidth = 5, label = '3D, SDE model solutions')

data = np.loadtxt("./data/data_stereo_a.csv", delimiter = ',')
data[:,1] = data[:,1]/aap
subplot1.scatter(data[:,0] - 2.5,data[:,1], edgecolor = 'red', facecolor = 'none', alpha = 0.9, label = 'Observations')

data = np.loadtxt("./data/data_earth.csv", delimiter = ',')
data[:,1] = data[:,1]/aap
subplot1.scatter(data[:,0] - 2.5,data[:,1], edgecolor = 'green', facecolor = 'none', alpha = 0.9)

data = np.loadtxt("./data/data_stereo_b.csv", delimiter = ',')
data[:,1] = data[:,1]/aap
subplot1.scatter(data[:,0] - 2.5,data[:,1], edgecolor = 'blue', facecolor = 'none', alpha = 0.9)

plt.legend(loc = 4)

subplot1 = fig.add_subplot(614)

#subplot1.set_title('Anisotropy vs time')
subplot1.set_xlabel('Time (hrs)')
subplot1.set_ylabel('Anisotropy')
subplot1.set_yscale('linear')
subplot1.set_xscale('linear')
subplot1.set_ylim([-0.25,2])
subplot1.set_xlim([-2.5,20])

data = np.loadtxt("./Test_4_Droge_benchmark_correct_position/STA_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'red', linestyle = '-', label = 'STEREO A', linewidth = 1)

data = np.loadtxt("./Test_5_Droge_benchmark_footpoint_position/STA_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'red', linestyle = '--', label = 'STEREO A', linewidth = 1)

data = np.loadtxt("./data/anisotropy_stereo_a.csv", delimiter = ',')
subplot1.plot(data[:,0] - 2.5,data[:,1], color = 'red', linestyle = '-', alpha = 0.5, linewidth = 5)

data = np.loadtxt("./data/data_anisotropy_stereo_a.csv", delimiter = ',')
subplot1.scatter(data[:,0] - 2.5,data[:,1], edgecolor = 'red', alpha = 0.9, facecolor = 'none')

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

data = np.loadtxt("./Test_4_Droge_benchmark_correct_position/STB_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'blue', linestyle = '-', label = 'STEREO B', linewidth = 1)

data = np.loadtxt("./Test_5_Droge_benchmark_footpoint_position/STB_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'blue', linestyle = '--', label = 'STEREO A', linewidth = 1)

data = np.loadtxt("./data/anisotropy_stereo_b.csv", delimiter = ',')
subplot1.plot(data[:,0] - 2.5,data[:,1], color = 'blue', linestyle = '-', alpha = 0.5, linewidth = 5)

data = np.loadtxt("./data/data_anisotropy_stereo_b.csv", delimiter = ',')
subplot1.scatter(data[:,0] - 2.5,data[:,1], edgecolor = 'blue', alpha = 0.9, facecolor = 'none')

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

data = np.loadtxt("./Test_4_Droge_benchmark_correct_position/Earth_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'green', linestyle = '-', label = 'Earth', linewidth = 1)

data = np.loadtxt("./Test_5_Droge_benchmark_footpoint_position/Earth_omni_time.txt")
subplot1.plot(data[:,2],data[:,4], color = 'green', linestyle = '--', label = 'STEREO A', linewidth = 1)

data = np.loadtxt("./data/anisotropy_earth.csv", delimiter = ',')
subplot1.plot(data[:,0] - 2.5,data[:,1], color = 'green', linestyle = '-', alpha = 0.5, linewidth = 5)

data = np.loadtxt("./data/data_anisotropy_earth.csv", delimiter = ',')
subplot1.scatter(data[:,0] - 2.5,data[:,1], edgecolor = 'green', alpha = 0.9, facecolor = 'none')

subplot1.axhline(y = 0, color = 'black')

plt.legend(loc = 1)

#plt.legend(loc = 1)

fig.savefig('output_2.png', dpi=300, bbox_inches='tight')

#------------------------------------------------------------------------





plt.show()
