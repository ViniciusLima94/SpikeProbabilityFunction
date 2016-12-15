# Routine to determine the probability of emiting a spike, from neuronal data (Membrane potential curve)

import pandas as pd
from numpy import *
import peakutils
from operator import truediv
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

data_matrix2 = []

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
			if P[aux] < -30:
				thr_values.append(P[aux])
				data_matrix2.append( P[P <= P[aux]] )

bl  = 1.0            # Length of bins
v_m = []             # Potencial value in center of each bar
bins_pot = []        # Distribution of potential values
bins_potdisp = []    # Distribution of threshold values

for v in arange(-70,1,bl):
	v_m.append( v + bl/2 )
	bins_pot.append( 0 )
	bins_potdisp.append( 0 )
	for i in range(0, len(data_matrix2)):
		for j in range(0, len(data_matrix2[i])):
			if data_matrix2[i][j] >= v and data_matrix2[i][j] < v+bl:
				bins_pot[-1] += 1
	for i in range(0, len(thr_values)):
		if thr_values[i] >= v and thr_values[i] < v+bl:
			bins_potdisp[-1] += 1

v_c = []
phi_v = []

for i in range(0, len(bins_pot)):
	if bins_pot[i] != 0:
		phi_v.append( float(bins_potdisp[i]) /  float(bins_pot[i]))
		v_c.append(v_m[i])

plot(v_c,phi_v)
ylabel('Probability')
xlabel('Membrane Potencial [mV]')
savefig('phi_v.png', dpi=600)
show()
