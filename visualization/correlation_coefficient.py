import numpy as np

def correlation_coefficient(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.shape != y.shape:
        raise ValueError("Input lists must have the same length")

    return np.corrcoef(x, y)[0, 1]


hamiltonian = [
        7.56710, 5.15426, 7.73095, 7.37927, 5.65139, 5.93459,
        14.00290, 4.08783, 9.69049, 5.41153, 3.71193, 5.12941,
        5.06777, 3.37614, 5.49035, 5.56518, 1.81516, 4.30766,
        3.99741, 1.25056, 4.83029, 4.65137, 3.77703, 5.10692
    ]

round_trip = [
        1.73441e-03, 2.15947e-02, 1.10226e-06,
        3.39271e-06, 2.10479e-02, 2.94123e-06,
        5.78864e-10, 6.52330e-03, 5.14535e-09,
        6.98536e-04, 1.67912e+00, 2.51654e-06,
        5.10312e-02, 1.98606e+00, 1.77416e-05,
        1.14103e-04, 2.14778e+00, 6.62741e-03,
        1.23745e-01, 2.18191e+00, 2.66206e-05,
        3.57416e-04, 2.57443e+00, 1.55446e-04
    ]

print(correlation_coefficient(hamiltonian, round_trip))