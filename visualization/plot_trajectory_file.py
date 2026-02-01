import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D


def plot_orbit(file, output_path=None):
    df = pd.read_csv(file)

    t = df["t"]
    x = df["x"]
    y = df["y"]
    z = df["z"]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    sc = ax.scatter(x, y, z, c=t, cmap="viridis", s=1)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title(f"3D Trajectory of {file}")

    fig.colorbar(sc, ax=ax, label="Time")

    if output_path is not None:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show(block=False)


def plot_two_orbits_planar(file1, file2, output_path=None):
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    def make_lc(df, color):
        x = df["x"].to_numpy()
        y = df["y"].to_numpy()

        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = LineCollection(segments, colors=color, linewidths=2)
        return lc

    lc1 = make_lc(df1, "tab:blue")
    lc2 = make_lc(df2, "tab:orange")

    fig, ax = plt.subplots()

    ax.add_collection(lc1)
    ax.add_collection(lc2)

    ax.autoscale()
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title("Planar Trajectories")

    # Legend (LineCollection does not auto-register)
    legend_handles = [
        Line2D([0], [0], color="tab:blue", lw=2, label=file1),
        Line2D([0], [0], color="tab:orange", lw=2, label=file2),
    ]
    ax.legend(handles=legend_handles)

    if output_path is not None:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show(block=False)

def plot_orbit_planar(file, output_path=None):
    df = pd.read_csv(file)

    t = df["t"].to_numpy()
    x = df["x"].to_numpy()
    y = df["y"].to_numpy()

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(
        segments,
        cmap="viridis",
        norm=plt.Normalize(t.min(), t.max())
    )
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

    if output_path is not None:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show(block=False)


def plot_radius(file, output_path=None):
    df = pd.read_csv(file)

    t = df["t"]
    x = df["x"]
    y = df["y"]
    z = df["z"]

    r = np.sqrt(x**2 + y**2 + z**2)

    fig, ax = plt.subplots()
    ax.scatter(t, r, s=5)

    ax.set_xlabel("Time")
    ax.set_ylabel("Radius")
    ax.set_title(f"Radius of {file}")

    if output_path is not None:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
