from dataclasses import field
from tkinter import N
from traceback import print_stack
import numpy as np
from scipy.integrate import odeint

import matplotlib.pyplot as plt


# Define the differential equation dy/dt = -y
def model(y, t):
    return -y


# Time points
t = np.linspace(0, 5, 100)

# Initial condition
y0 = 1

# Solve ODE
solution = odeint(model, y0, t)

# Plot results
plt.figure(figsize=(8, 6))
plt.plot(t, solution, "b-", label="y(t)")
plt.xlabel("Time")
plt.ylabel("y(t)")
plt.title("Solution to dy/dt = -y")
plt.grid(True)
plt.legend()
plt.show()

import pandas as pd

df = pd.DataFrame({"t": t, "y": solution[:, 0]})
mean_values = df.aggregate("mean")
df.between_time("0:00", "0:01")
print(mean_values)
list(df.columns)

import polars as pl
