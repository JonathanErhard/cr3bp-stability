import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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

    t = df["t"]
    x = df["x"]
    y = df["y"]

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.scatter(x, y, c=t, cmap='viridis', label="Trajectory", linewidth=2)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title(f"planar Trajectory of {file}")
    ax.legend()

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

# plot_orbit_planar("data/CR3BP_L3_forward.csv")
# plot_orbit_planar("data/CR3BP_L3_backward.csv")
plot_orbit_planar("data/CR3BP_rnd_forward.csv")
plot_orbit_planar("data/CR3BP_rnd_backward.csv")
plot_orbit_planar("data/CR3BP.csv")
input("DO IT")