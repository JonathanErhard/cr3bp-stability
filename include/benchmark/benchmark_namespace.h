#pragma once

// Benchmarking namespace for the CR3BP benchmarks
namespace cr3bp_benchmarks{
    using state_type = std::array<double,6>;
    using trajectory_type = std::vector<std::pair<double,state_type>>;
};