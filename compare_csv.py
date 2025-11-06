import pandas as pd
import numpy as np

def compare_files(file1, file2):
    df1 = pd.read_csv(file1,header=None)
    df2 = pd.read_csv(file2,header=None)
    derivaties = []
    for index,row in df.iterrows():
        state = [np.float64(row['x']),np.float64(row['y']),np.float64(row['dx']),np.float64(row['dy'])]
        print(f"energy: {energy(state,mu)}")
        derivaties.append(prtbp(0,state,mu))
        #print(np.array2string(prtbp(0,state,1.215e-2),separator=',',precision=3)[1:-1])