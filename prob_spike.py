# Routine to determine the probability of emiting a spike, from neuronal data (Membrane potential curve)

import pandas as pd
from numpy import *
import peakutils
from operator import truediv
from matplotlib.pyplot import *

# Receives signal and threshold
# Returns peaks index for each signal
def findNpeaks(signal, th):
	Npeaks = 0 # Just count the number of peaks in a given experiment	
	index = [] # Will store peaks indexes
	# Find peaks, and store all peaks found
	index_aux = peakutils.indexes(signal)
	# For each peak found
	for idx in index_aux:
		# If this peak is greater than a given threshold
		if signal[idx] > th:
			# Store it on index list
			index.append(idx)
			# Count a peak
			Npeaks += 1
	# Return Npeaks and index
	return Npeaks, index

# Reading data
data_matrix = pd.read_csv("dataMatrix.dat", delimiter=",", header=None).as_matrix()
# Matrix number of rows and columns
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
	# For each experiment, find peaks and its indexes
	np,idx = findNpeaks(data_matrix[i][:], -10)
	# Append those values to its respective lists
	Npeaks.append(np)
	index.append(idx)
	
# data_matrix2 is a cell each line of the cell will store the potencial
# values below the found threhold for a particular expriment
data_matrix2 = []

# Finding thresholds
for i in range(0,r):
	# If it have any peaks
	if index[i][:] != []:
		for j in range(0, len(index[i][:])):
			# Separate each peak, to find the threshold
           	        # P is one peak of a given experimente
            		# P is choosed from the peak to 3ms ago
			P = data_matrix[i][index[i][j]-60:index[i][j]]
			# dP/dt
			P1 = diff(P,1)
			P1 = P1[:-1]
			#(d2P/dt2)
			P2 = diff(P,2)
			# Method VII
			Kp = P2 * (1 + P1**2)**1.5
			# Find the max Kp index
			aux = argmax(Kp)
			# Append threshold found in thr_values
			thr_values.append(P[aux])
			# Append peak in data_matrix2, excluding all values
			# greater than the found threshold
			data_matrix2.append( P[P <= P[aux]] )

bl  = 1.0            # Length of bins
v_m = []             # Potencial value in center of each bar
bins_pot = []        # Distribution of potential values
bins_potdisp = []    # Distribution of threshold values

# Dividing into bins
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
	# In order to avoid for zero division
	if bins_pot[i] != 0:
		phi_v.append( float(bins_potdisp[i]) /  float(bins_pot[i]))
		v_c.append(v_m[i])

plot(v_c,phi_v)
ylabel('Probability')
xlabel('Membrane Potencial [mV]')
savefig('phi_v.png', dpi=600)
show()
