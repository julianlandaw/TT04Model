import numpy as np
import matplotlib.pyplot as plt
import sys

pcl = sys.argv[5]
var1 = sys.argv[1]
var2 = sys.argv[2]
var3 = sys.argv[3]
var4 = sys.argv[4]
txtend = str(pcl) + '_VAR1_' + str(var1) + '_VAR2_' + str(var2) + '_VAR3_' + str(var3) + '_VAR4_' + str(var4)

A = np.loadtxt('TTapPCL_' + txtend + '.txt')
apds = np.loadtxt('TTbifsPCL_' + txtend + '.txt')
nais = np.loadtxt('TTnaiPCL_' + txtend + '.txt')
cais = np.loadtxt('TTcaiPCL_' + txtend + '.txt')

plt.figure(figsize=(20,10))

plt.subplot(2,2,1)
plt.plot(A[1:25000,0],A[1:25000,1])
plt.xlabel('Time (ms)');
plt.ylabel('Voltage (mV)')

plt.subplot(2,2,2);
plt.plot(apds[5:],'.')
plt.xlabel('Beat Num');
plt.ylabel('APD (ms)');

inds = np.arange(5,len(cais),2)

plt.subplot(2,2,3);
plt.plot(1000*cais[inds],'.')
plt.xlabel('Beat Num');
plt.ylabel('[Ca]_i (uM)');

plt.subplot(2,2,4);
plt.plot(nais[inds],'.')
plt.xlabel('Beat Num');
plt.ylabel('[Na]_i (mM)');

plt.savefig('TTplot_' + txtend + '.png');