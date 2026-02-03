import pandas as pd
import matplotlib.pyplot as plt

def plot_error(file,output_path=None):
    df = pd.read_csv(file)
    t = df["t"].to_numpy()
    rel_err = df["rel"].to_numpy()
    abs_err = df["abs"].to_numpy()

    fig = plt.figure()
    plt.plot(t, abs_err, label="Absolute Error")
    plt.plot(t, rel_err, label="Relative Error")

    plt.xlabel("t [s]")
    plt.ylabel("errors")
    plt.title("Absolute and Relative errors over time")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()

    if output_path is not None:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show(block=False)

def plot_abs_error(file,output_path=None):
    df = pd.read_csv(file)
    t = df["t"].to_numpy()
    abs_err = df["abs"].to_numpy()

    fig = plt.figure()
    plt.plot(t, abs_err, label="Absolute Error")

    plt.xlabel("t [s]")
    plt.ylabel("abs error")
    plt.title("Absolute error over time")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()

    if output_path is not None:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show(block=False)