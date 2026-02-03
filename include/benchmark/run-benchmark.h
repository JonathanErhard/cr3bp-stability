    std::stringstream ss;
    ss << " benchmark results:\n";
    ss << "round-trip closure and hamiltonian:\n";
// loop through round-back-closure starting states
    for(int index = 0;index<INTEGRATOR_PARAMETERS::rbc::starting_positions.size();index++){
        auto x0 = INTEGRATOR_PARAMETERS::rbc::starting_positions[index];

        const auto& [fwd_traj,bwd_traj] = cr3bp_benchmarks::roundtrip_closure_benchmark(
            cr3bp_model,
            subject,
            INTEGRATOR_PARAMETERS::rbc::starting_positions[index],
            0.0,
            INTEGRATOR_PARAMETERS::integration_time,
            INTEGRATOR_PARAMETERS::grid_resolution
        );
        //auto l2_distances = compare_trajectories_isochronous(fwd_traj,bwd_traj);
        //double last_L2 = l2_distances.back();
        double last_L2 = l2_error(fwd_traj.begin()->second, bwd_traj.back().second);
        
        // calculate hamiltonian error of the forward trajectory
        auto hamiltonian_error = cr3bp_benchmarks::hamiltonian_conservation_benchmark(fwd_traj, INTEGRATOR_PARAMETERS::mu);
        double max_hamiltonian_error = std::max_element(
            hamiltonian_error.begin(),
            hamiltonian_error.end(),
            [](const cr3bp_benchmarks::numeric_error& a, const cr3bp_benchmarks::numeric_error& b){
                return abs(a.rel_error) < abs(b.rel_error);
            }
        )->rel_error;

        ss << INTEGRATOR_PARAMETERS::rbc::names[index] << ": last L2 error = " << last_L2 << "\n";
        ss << INTEGRATOR_PARAMETERS::rbc::names[index] << ": max hamiltonian error (log)= " << -std::log10(max_hamiltonian_error) << "\n";

        // save files for plots
        std::string directory_path = file_path_prefix + '/' + INTEGRATOR_PARAMETERS::rbc::names[index];
        
        std::filesystem::create_directories(directory_path);
        
        save_trajectory(fwd_traj, directory_path + "/fw.csv");
        save_trajectory(bwd_traj, directory_path + "/bw.csv");
        
        save_numeric_error(hamiltonian_error, directory_path + "/hamiltonian_error.csv");
    }

    ss << "\nsurrogate_p1 model benchmark:\n";

    for(int index = 0;index<INTEGRATOR_PARAMETERS::surrogate_p1::starting_positions.size();index++){
        auto x0 = INTEGRATOR_PARAMETERS::surrogate_p1::starting_positions[index];

        // loop through surrogate_p1 starting positions
        const auto& [exact_traj,estimated_traj] = cr3bp_benchmarks::surrogate_p1_benchmark(
            p1_surrogate_model,
            subject,
            x0,
            0,
            INTEGRATOR_PARAMETERS::integration_time,
            INTEGRATOR_PARAMETERS::grid_resolution
        );
            
        auto l2_distances = compare_trajectories_isochronous(exact_traj,estimated_traj);
        double last_L2 = l2_distances.back().abs_error;

        // safe results to file
        std::string directory_path = file_path_prefix + '/' + INTEGRATOR_PARAMETERS::surrogate_p1::names[index];
        
        std::filesystem::create_directories(directory_path);
        
        ss << INTEGRATOR_PARAMETERS::surrogate_p1::names[index] << ": last error = " << last_L2 << "\n";
        // save trajectories
        save_trajectory(exact_traj, directory_path + "/exact_surrogate.csv");
        save_trajectory(estimated_traj, directory_path + "/estimated_surrogate.csv");

    }
    std::cout << "Files written to " << file_path_prefix << "/<csv-files>.\n";
    std::cout << integrator_name <<  ss.str() << "\n####################\n\n";