import pandas as pd
import numpy as np
import time
from HS2SWE import HS2SWE

# Load the data from the file
# for single station testfile comparable to a reference simulation and measurement data
data_wfj = pd.read_csv("../data/data_weissfluhjoch.txt", delimiter=',')  # Use the correct delimiter if it's tab-separated
Inputdata=np.transpose(np.array(data_wfj['hs_obs'].values,ndmin=2))   #add a dimension of length 1 for point data use in python

#run the model 
swe_sim = HS2SWE(Inputdata);


#plotting routine 
import matplotlib.pyplot as plt

# Assuming data_wfj is a pandas DataFrame and swe_sim is a numpy array or pandas Series

# Create a figure with specific size (in centimeters)
plt.figure(figsize=(28/2.54, 12/2.54))  # Convert cm to inches (1 inch = 2.54 cm)

# Plot the reference SWE data (swe_ref in gray)
plt.plot(data_wfj['time'], data_wfj['swe_ref'], color='#808080', label='Reference SWE')

# Plot the simulated SWE (swe_sim in red)
plt.plot(data_wfj['time'], swe_sim, '.', color='#f54209', label='Simulated SWE')

# Plot the observed SWE (swe_obs in dark red, larger markers)
plt.plot(data_wfj['time'], data_wfj['swe_obs'], '.', color='#1e0000', markersize=12, label='Observed SWE')

# Label the y-axis
plt.ylabel("Snow water equivalent [mm]")

# Add a title
plt.title("Simulations for the Weissfluhjoch field site at 2540 m.a.s.l. near Davos, Switzerland")
plt.savefig('figures/temporal_plot')

# diagnostic check before calucalting residuals between the reference and your simulation
data_wfj['swe_ref'].shape
swe_sim.shape

residual= np.transpose(np.array(data_wfj['swe_ref'].values,ndmin=2))-swe_sim

# Create a figure with specific size (in centimeters)
plt.figure(figsize=(28/2.54, 12/2.54))  # Convert cm to inches (1 inch = 2.54 cm)

# Plot the residual between reference and your simulation
plt.plot(data_wfj['time'],residual, color='#808080', label='Reference SWE')

# Label the y-axis
plt.ylabel("Snow water equivalent [mm]")

# Add a title
plt.title("Residuals between your and the reference Simulation for Weissfluhjoch field site at 2540 m.a.s.l. near Davos, Switzerland")
plt.savefig('figures/residual_plot')





import pandas as pd
import numpy as np
import time
'''
#for investigating multiple stations test files 
# Loading data from multidimensional file station file (HS data overtime for 1:nstations)
data_multidim=np.loadtxt('../data/export_stations_no_header.txt',delimiter='\t',skiprows=1,usecols=np.arange(6,450)) #skiping non-HS data columns 
Inputdata=data_multidim            
Inputdata[Inputdata<0]=0

#skipping day1 due to multiple station without proper data
swe_sim = HS2SWE(Inputdata[1:,:]); 
#numpy.savetxt('output_all_stations_python.txt, swe_sim, fmt='%.18e', delimiter=',')


import matplotlib.pyplot as plt
# plot your SWE data:
plt.figure(figsize=(28/2.54, 12/2.54))  # Convert cm to inches (1 inch = 2.54 cm)
for i in range(swe_sim.shape[1]):
    plt.plot(np.arange(0,9669),swe_sim[:,i], color='#1e0000')
# Label the y-axis
plt.ylabel("Snow water equivalent [mm]")
plt.title("SWE for stations calculated using HS2SWE from snow height measurements")
plt.savefig('figures/temporal_plot_multi_station')

# possible compare to other simulations calculating residuals between the simulations. 
# load reference dataset from a malab simulation 
matlab_swe_sim=np.loadtxt('../data/444station_output_matlab.txt',delimiter=',')
python_swe_sim=np.loadtxt('../data/output_all_stations_python.txt',delimiter=',')

# calculate residuals relative to other data  #select what you want to compare
residual=matlab_swe_sim-swe_sim
#residual=python_swe_sim-swe_sim


# Create a figure with specific size (in centimeters)
plt.figure(figsize=(28/2.54, 12/2.54))  # Convert cm to inches (1 inch = 2.54 cm)

# Plot the reference SWE data (swe_ref in gray)
plt.plot(np.arange(0,9669),residual[:,0], color='#808080')

for i in range(swe_sim.shape[1]):
    plt.plot(np.arange(0,9669),residual[:,i], color='#1e0000')
# Label the y-axis
plt.ylabel("Snow water equivalent [mm]")

# Add a title
plt.title("Residuals for all stations")
plt.savefig('figures/residuals_all_stations')
'''
