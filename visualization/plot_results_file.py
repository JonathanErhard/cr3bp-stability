import sys

import plot_trajectory_file
import plot_error

integrator_name = sys.argv[1]

directories_rbc = ["halo_stable","halo_unstable","gateway_L2_southern","gateway_L2_northern","transfer_orbit","distant_retrograde"]
directories_surrogate_p1 = ["orbit_right","orbit_above"]

dir_prefix = f"../benchmark_output/{integrator_name}"

for directory in directories_rbc:
    plot_trajectory_file.plot_orbit_planar(f"{dir_prefix}/{directory}/bw.csv", output_path=f"../benchmark_output/{integrator_name}/{directory}/bw.png")
    plot_trajectory_file.plot_orbit_planar(f"{dir_prefix}/{directory}/fw.csv", output_path=f"../benchmark_output/{integrator_name}/{directory}/fw.png")
    plot_trajectory_file.plot_two_orbits_planar(f"{dir_prefix}/{directory}/bw.csv", f"{dir_prefix}/{directory}/fw.csv", output_path=f"../benchmark_output/{integrator_name}/{directory}/both.png")
    plot_error.plot_error(f"../benchmark_output/{integrator_name}/{directory}/hamiltonian_error.csv", f"../benchmark_output/{integrator_name}/{directory}/hamiltonian_error.png")

for directory in directories_surrogate_p1:
    plot_trajectory_file.plot_orbit_planar(f"{dir_prefix}/{directory}/estimated_surrogate.csv", output_path=f"../benchmark_output/{integrator_name}/{directory}/estimated_surrogate.png")
    plot_trajectory_file.plot_orbit_planar(f"{dir_prefix}/{directory}/exact_surrogate.csv", output_path=f"../benchmark_output/{integrator_name}/{directory}/exact_surrogate.png")
    plot_trajectory_file.plot_two_orbits_planar(f"{dir_prefix}/{directory}/estimated_surrogate.csv", f"{dir_prefix}/{directory}/exact_surrogate.csv", output_path=f"../benchmark_output/{integrator_name}/{directory}/estimated_vs_exact.png")
# input("press any key")