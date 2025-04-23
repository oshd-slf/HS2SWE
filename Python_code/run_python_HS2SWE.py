import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from HS2SWE import HS2SWE

def single_station_run():
    # Load the data from the file
    data_wfj = pd.read_csv("../data/data_weissfluhjoch.txt", delimiter=',')  # Adjust delimiter if needed
    Inputdata = np.transpose(np.array(data_wfj['hs_obs'].values, ndmin=2))  # Add a dimension for point data

    # Run the model
    swe_sim = HS2SWE(Inputdata)

    # Plotting routine
    plt.figure(figsize=(28 / 2.54, 12 / 2.54))  # Convert cm to inches
    plt.plot(data_wfj['time'], data_wfj['swe_ref'], color='#808080', label='Reference SWE')
    plt.plot(data_wfj['time'], swe_sim, '.', color='#f54209', label='Simulated SWE')
    plt.plot(data_wfj['time'], data_wfj['swe_obs'], '.', color='#1e0000', markersize=12, label='Observed SWE')
    plt.ylabel("Snow water equivalent [mm]")
    plt.title("Simulations for the Weissfluhjoch field site at 2540 m.a.s.l. near Davos, Switzerland")
    plt.savefig('figures/temporal_plot')

    # Calculate residuals
    residual = np.transpose(np.array(data_wfj['swe_ref'].values, ndmin=2)) - swe_sim

    # Plot residuals
    plt.figure(figsize=(28 / 2.54, 12 / 2.54))  # Convert cm to inches
    plt.plot(data_wfj['time'], residual, color='#808080', label='Residuals')
    plt.ylabel("Snow water equivalent [mm]")
    plt.title("Residuals between your and the reference Simulation for Weissfluhjoch field site")
    plt.savefig('figures/residual_plot')

def multiple_station_run():
    # Load data from multidimensional file
    data_multidim = np.loadtxt('../data/export_stations_no_header.txt', delimiter='\t', skiprows=1, usecols=np.arange(6, 450))
    Inputdata = data_multidim
    Inputdata[Inputdata < 0] = 0  # Replace negative values with 0

    # Skip day 1 due to missing data
    swe_sim = HS2SWE(Inputdata[1:, :])

    # Plot SWE data
    plt.figure(figsize=(28 / 2.54, 12 / 2.54))  # Convert cm to inches
    for i in range(swe_sim.shape[1]):
        plt.plot(np.arange(0, swe_sim.shape[0]), swe_sim[:, i], color='#1e0000')
    plt.ylabel("Snow water equivalent [mm]")
    plt.title("SWE for stations calculated using HS2SWE from snow height measurements")
    plt.savefig('figures/temporal_plot_multi_station')
    
if __name__ == "__main__":
    print("Select mode:")
    print("1. Single station run")
    print("2. Multiple station run")
    choice = input("Enter 1 or 2: ")

    if choice == "1":
        single_station_run()
    elif choice == "2":
        multiple_station_run()
    else:
        print("Invalid choice. Please run the script again and select 1 or 2.")