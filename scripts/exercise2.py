# Exercise 2: Plotting MD Data

# In this exercise, we're going to plot the output of MD simulations using matplotlib library

# First we're going to extract and MD data from the log file using pandas

import pandas as pd
import matplotlib.pyplot as plt

# Read the data from the log file into a Pandas DataFrame
data = pd.read_csv("Cu_md.log", sep='\s+', skiprows=1)
data.columns = ['Time[ps]', 'Etot[eV]', 'Epot[eV]', 'Ekin[eV]', 'T[K]']

# Plot each column against Time[ps]
plt.figure(figsize=(10, 6))
plt.plot(data['Time[ps]'], data['T[K]'])
plt.xlabel('Time / ps')
plt.ylabel('Temperature / K')
plt.title('MD Simulation Data')
plt.legend()
plt.grid(True)
plt.show()

# Now it's your turn! Have a go at plotting the energy values againsts time, try and plot them together.