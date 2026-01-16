import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def plot_orbit(file):
    df = pd.read_csv(file)

    t = df["t"]
    x = df["x"]
    y = df["y"]
    z = df["z"]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x, y, z, c=t, cmap="viridis", s=1)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(f"3D Trajectory of {file}")
    ax.legend()

    plt.show(block=False)

def plot_orbit_planar(file):
    df = pd.read_csv(file)

    t = df["t"].to_numpy()
    x = df["x"].to_numpy()
    y = df["y"].to_numpy()

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap='viridis', norm=plt.Normalize(t.min(), t.max()))
    lc.set_array(t)
    lc.set_linewidth(2)

    fig, ax = plt.subplots()
    line = ax.add_collection(lc)
    ax.autoscale()
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title(f"Planar Trajectory of {file}")

    cbar = fig.colorbar(line, ax=ax)
    cbar.set_label("Time")

    plt.show(block=False)
    
def plot_radius(file):
    df = pd.read_csv(file)

    t = df["t"]
    x = df["x"]
    y = df["y"]
    z = df["z"]

    r = r = np.sqrt(x**2 + y**2 + z**2)

    fig = plt.figure()
    ax = fig.add_subplot()

    ax.scatter(t, r)

    ax.set_xlabel("radius")
    ax.set_ylabel("time")
    ax.set_title(f"radius of {file}")
    ax.legend()

    plt.show()

#plot_orbit_planar("../benchmark_output/Dormand-Prince-5/bw_1.csv")
#plot_orbit_planar("../benchmark_output/Dormand-Prince-5/fw_1.csv")
#plot_orbit_planar("../benchmark_output/Dormand-Prince-5/both.csv")
#
#plot_orbit_planar("../benchmark_output/Dormand-Prince-5/estimated_surrogate.csv")
#plot_orbit_planar("../benchmark_output/Dormand-Prince-5/exact_surrogate.csv")
#plot_orbit_planar("../benchmark_output/Dormand-Prince-5/exact_and_estimated_surrogate.csv")

# plot

# input("DO IT")