import pandas as pd
import numpy as np

def prtbp(t, x, mu):
    mu1 = 1 - mu  # mass of larger primary (nearest origin on left)
    mu2 = mu      # mass of smaller primary (furthest from origin on right)


    r3 = ((x[0] + mu2) ** 2 + x[1] ** 2) ** 1.5  # r: distance to m1, LARGER MASS
    R3 = ((x[0] - mu1) ** 2 + x[1] ** 2) ** 1.5  # R: distance to m2, smaller mass
    print(f"r1sqr= {((x[0] + mu2) ** 2 + x[1] ** 2)}")
    print(f"r1: {r3}")
    print(f"r2: {R3}")
    print(f"state: {x}")

    xdot = np.zeros(4)
    xdot[0] = x[2]
    xdot[1] = x[3]
    xdot[2] = x[0] - (mu1 * (x[0] + mu2) / r3) - (mu2 * (x[0] - mu1) / R3) + 2 * x[3]
    xdot[3] = x[1] - (mu1 * x[1] / r3) - (mu2 * x[1] / R3) - 2 * x[2]

    print(f"x_dot: {xdot}")
    
    return xdot



# Energy function for the CR3BP system
def energy(X, mu):
    # mu1=1-mu;mu2=mu;
    mu1 = 1 - mu
    mu2 = mu

    # Vsqrd= X(:,3).^2+X(:,4).^2 ;
    Vsqrd = X[2]**2 + X[3]**2

    # Ubar = - mu1./sqrt((X(:,1)+mu2).^2+X(:,2).^2) ...
    #        - mu2./sqrt((X(:,1)-mu1).^2+X(:,2).^2) ...
    #        - 0.5*(X(:,1).^2+X(:,2).^2) - 0.5*mu1*mu2 ;
    Ubar = - mu1 / np.sqrt((X[0] + mu2)**2 + X[1]**2) \
           - mu2 / np.sqrt((X[0] - mu1)**2 + X[1]**2) \
           - 0.5 * (X[0]**2 + X[1]**2) - 0.5 * mu1 * mu2

    # E = 0.5*Vsqrd + Ubar ;
    E = 0.5 * Vsqrd + Ubar

    # end
    return E


def run_pertbp(file):
    df = pd.read_csv(file)
    derivaties = []
    for index,row in df.iterrows():
        state = [np.float64(row['x']),np.float64(row['y']),np.float64(row['dx']),np.float64(row['dy'])]
        print(f"energy: {energy(state,mu)}")
        derivaties.append(prtbp(0,state,mu))
        #print(np.array2string(prtbp(0,state,1.215e-2),separator=',',precision=3)[1:-1])

# main
mu = 1.215e-2
np.set_printoptions(legacy='1.25')
file = "build/CR3BP.csv"
run_pertbp(file)