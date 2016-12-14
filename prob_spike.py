# Routine to determine the probability of emiting a spike, from neuronal data (Membrane potential curve)

import pandas as pd
from numpy import *
import peakutils
from matplotlib.pyplot import *

# Receives signal and threshold
# Returns peaks index for each signal
def findNpeaks(signal, th):
	Npeaks = 0	
	index = []
	index_aux = peakutils.indexes(signal)
	for idx in index_aux:
		if signal[idx] > th:
			index.append(idx)
			Npeaks += 1
	return Npeaks, index

data_matrix = pd.read_csv("dataMatrix.dat", delimiter=",", header=None).as_matrix()
r,c = np.shape(data_matrix)

# To plot data uncomment
# figure(1)
# for i in range(0,r):
# 	plot(data_matrix[i][:])
# 	hold(True)
# show()

# Matrix number of rows and columns
r,c = np.shape(data_matrix)

# Store number of peaks in each experiment
Npeaks = []

# Store peak index
index = []

# Store threshold values
thr_values = []

# Finding peaks for each signal
for i in range(0,r):
	np,idx = findNpeaks(data_matrix[i][:], -10)
	Npeaks.append(np)
	index.append(idx)

# Finding thresholds
for i in range(0,r):
	if index[i][:] != []:
		for j in range(0, len(index[i][:])):
			P = data_matrix[i][index[i][j]-60:index[i][j]]
			# dP/dt
			P1 = diff(P,1)
			P1 = P1[:-1]
			#(d2P/dt2)
			P2 = diff(P,2)
			# Method VII
			Kp = P2 * (1 + P1**2)**1.5
			aux = argmax(Kp)
			thr_values.append(P[aux])

# Saving all potencial values in a vector
v = []
for i in range(0,r):
	for j in range(0,c):
		v.append(data_matrix[i][j])

nbins = 50
counts, center = histogram(v,nbins)
# Bin Length
bin = diff(center)
bin = bin[0]
# Creating threshold potential values
counts_1 = []
for i in range(0, len(center)):
	aux = 0
	for j in range(0, len(thr_values)):
		if thr_values[j] >= center[i]-bin/2 and thr_values[j] < center[i]+bin/2:
			aux += 1
	counts_1.append( aux )
counts_1 = counts_1[:-1]
# Probability density
gamma = 1 / sum(map(float,counts_1)/counts) 
ro = map(float,counts_1)/counts
ro = gamma * ro
# phi_v
phi_v = cumsum(ro)
plot(center[:-1],phi_v)
ylabel('Probability')
xlabel('Membrane Potential [mV]')
title('Spike probability')
show()
#savefig('prob_disp2.png', dpi=600)



