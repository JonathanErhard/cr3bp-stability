import sys

import plot_trajectory
import plot_error

integrator_name = sys.argv[1]

directories_rbc = ["halo_stable","halo_unstable","gateway_L2_southern","gateway_L2_northern"]
directories_surrogate_p1 = ["orbit_right","orbit_above"]

for directory in directories_rbc:
    # plot_trajectory.plot_orbit_planar(f"../benchmark_output/{integrator_name}/{directory}/bw_1.csv")
    # plot_trajectory.plot_orbit_planar(f"../benchmark_output/{integrator_name}/{directory}/fw_1.csv")
    # plot_trajectory.plot_orbit_planar(f"../benchmark_output/{integrator_name}/{directory}/fw_bw.csv")
    plot_error.plot_error(f"../benchmark_output/{integrator_name}/{directory}/hamiltonian_error.csv")

for directory in directories_surrogate_p1:
    pass
    # plot_trajectory.plot_orbit_planar(f"../benchmark_output/{integrator_name}/{directory}/estimated_surrogate.csv")
    # plot_trajectory.plot_orbit_planar(f"../benchmark_output/{integrator_name}/{directory}/exact_surrogate.csv")
    # plot_trajectory.plot_orbit_planar(f"../benchmark_output/{integrator_name}/{directory}/exact_and_estimated_surrogate.csv")

input("press any key")