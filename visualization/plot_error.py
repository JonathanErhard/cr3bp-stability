import pandas as pd
import matplotlib.pyplot as plt

def plot_error(file):
    df = pd.read_csv(file)
    t = df["t"]
    rel_err = df["rel"]
    abs_err = df["abs"]

    plt.figure()
    plt.plot(t.to_numpy(), abs_err.to_numpy(), label="Absolute Error")
    plt.plot(t.to_numpy(), rel_err.to_numpy(), label="Relative Error")

    plt.xlabel("t [s]")
    plt.ylabel("errors")
    plt.title("Absolute and Relative errors over time")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show(block=False)
