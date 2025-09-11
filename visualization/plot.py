import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("data.csv",header=None)

_, axes = plt.subplots(1, 1, figsize=(9, 9))

axes.plot(data[0].to_numpy(), data[1].to_numpy(), linewidth=2, color='darkblue', ls='-', label='integral')

plt.show()