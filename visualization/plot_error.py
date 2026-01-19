import pandas as pd
import matplotlib.pyplot as plt

def plot_error(file):
    df = pd.read_csv(file)
    t = df["t"].to_numpy()
    rel_err = df["rel"].to_numpy()
    abs_err = df["abs"].to_numpy()

    plt.figure()
    plt.plot(t, abs_err, label="Absolute Error")
    plt.plot(t, rel_err, label="Relative Error")

    plt.xlabel("t [s]")
    plt.ylabel("errors")
    plt.title("Absolute and Relative errors over time")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show(block=False)
