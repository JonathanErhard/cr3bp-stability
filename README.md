
# cr3bp stability

This repository implements a stability benchmarking tool as well as implementations of cr3bp dynamics and explicit and implicit RK-integrators

## HowTo compile

required libraries:
eigen, heyoka, boost, g++15, c++23

Compile benchmarks using

```bash
mkdir build && cd build
cmake ..
make
```

if CMake struggles to find a library, add

```bash
-DLIBRARY_NAME_ROOT=/path/to/folder
```

to your cmake command

example:

```bash
cmake -DBOOST_ROOT=/usr/include/boost ..
```

I have created a file TROUBLESHOOT.md in case make struggles to link to heyokas dependencies.

## HowTo run benchmarks

To run the benchmarks just run the executables with

```bash
./benchmark-heyoka
./benchmark-fehlberg78
./benchmark-DP5
./benchmark-bulirsch_stoer
```

the data files will appear in benchmark_output in csv format

## HowTo plot

There are two files for plotting the results. 

### plot_results_file.py

Saves the plots to a file next to the corresponding csv files

run using

```bash
cd util
python3 plot_results_file.py integrator_name
```

### plot_results.py

Shows the plots instead of saving them to a file
run using

```bash
cd util
python3 plot_results.py integrator_name
```


## files

### src

The src folder contains benchmark files in ~benchmark/*~ , test-files to ensure the the correctness of the implementation as well as files to propagate an orbit using different integrators.

### Util

The util folder contains files to plot results and a script to calculate the correlation coefficient of two lists. It also contains a data folder with .csv files of reference orbits. Note that some of the scripts were written with AI assistence, since reading matlibplot documentation is not my passion

### include

the include folder contains header files that implement the numerical integrators as well as header files for the mathematical model and files for the benchmarking suite

### benchmark_output

Contains csv and plots of the orbits and their associated errors