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
        double max_l2 = compare_trajectories_isochronous(fwd_traj,bwd_traj);
        
        // calculate hamiltonian error of the forward trajectory
        auto hamiltonian_error = cr3bp_benchmarks::hamiltonian_conservation_benchmark(fwd_traj, INTEGRATOR_PARAMETERS::mu);
        
        // save files for plots
        std::string directory_path = file_path_prefix + '/' + INTEGRATOR_PARAMETERS::rbc::names[index];
        
        std::filesystem::create_directories(directory_path);
        
        save_trajectory(fwd_traj, directory_path + "/fw.csv");
        save_trajectory(bwd_traj, directory_path + "/bw.csv");

        // trajectory_type joined_traj;
        // std::copy(fwd_traj.begin(),fwd_traj.end(),std::back_inserter(joined_traj));
        // joined_traj.insert(joined_traj.end(),bwd_traj.begin(),bwd_traj.end());
        // save_trajectory(joined_traj, directory_path + "/fw_bw.csv");
        
        save_numeric_error(hamiltonian_error, directory_path + "/hamiltonian_error.csv");
    }

    for(int index = 0;index<INTEGRATOR_PARAMETERS::surrogate_p1::starting_positions.size();index++){
        auto x0 = INTEGRATOR_PARAMETERS::surrogate_p1::starting_positions[index];

        // loop through surrogate_p1 starting positions
        const auto& [exact_traj,estimated_traj] = cr3bp_benchmarks::surrogate_p1_benchmark(
            p1_surrogate_model,
            subject,
            x0,
            0,
            100,
            INTEGRATOR_PARAMETERS::grid_resolution
        );
            
        // safe results to file
        std::string directory_path = file_path_prefix + '/' + INTEGRATOR_PARAMETERS::surrogate_p1::names[index];
        
        std::filesystem::create_directories(directory_path);
        
        // save trajectories
        save_trajectory(exact_traj, directory_path + "/exact_surrogate.csv");
        save_trajectory(estimated_traj, directory_path + "/estimated_surrogate.csv");
        
        // save both trajectories one after the other for plotting purposes
        // trajectory_type t1;
        // std::copy(exact_traj.begin(),exact_traj.end(),std::back_inserter(t1));
        // t1.insert(t1.end(),estimated_traj.begin(),estimated_traj.end());
        // save_trajectory(t1, directory_path + "/exact_and_estimated_surrogate.csv");
    }
    std::cout << "Files written to " << file_path_prefix << "/<csv-files>.\n";