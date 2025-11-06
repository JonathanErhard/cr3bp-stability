
# N Body Propagation

This repository contains implementations for basic numeric integrators and a model for the CR3BP.

## examples

The examples folder contains 3 files. CR3CP_euler.cpp and CR3CP_rk.cpp calculate a trajectory using the model defined in include/model/CR3BP. The initial state, step_size and mu are hardcoded in the cpp file

RK-test.cpp contains a basic application of a RK3 integrator.

## HOWTO compile

To compile any of the cpp files create and navigate into the build directory

```bash
mkdir build && cd build
```

to compile CR3BP_euler.cpp, use

```bash
g++ ../examples/CR3BP_euler.cpp -I ../include/ -o CR3BP_euler.out
```

to compile CR3BP_rk.cpp, use

```bash
g++ ../examples/CR3BP_rk.cpp -I ../include/ -o CR3BP_rk.out
```

to compile RK-test.cpp, use

```bash
g++ ../examples/RK-test.cpp -I ../include/ -o RK-test.out
```

## HOWTO run

Following the compilation process, the executables are created in the build folder. If you run CR3BP_euler.out or CR3BP_rk.out, an outputfile is created in the build folder. It contains a timestamp followed by the state in rotating frame (bary-centered).

## HOWTO plot

To plot a trajectory use

```bash
cd visualization
mkdir data
mv ../build/CR3BP.csv data
python3 plot_trajectory.py
```

By default, a projection onto the x-y plane of the trajectory is plotted in the rotating frame. The filename and plotting method can be configured in visualization/plot_trajectory.py