# nsgaii

In this repository there are two implementations of the 
Nondominated Sorting Genetic Algorithm as it is described in the paper of [Deb et al., 2000](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=996017)

Out intention is not to write an optimization algorithm that includes all possible options a genetic algorithm may have, but rather we include only the parts that are work in the problems we are facing on our labs. Therefore, from the dozens of available selection, crossover nad mutation operators we have selected only a few that work well in our cases.

At the moment we provide two implementations. The Matlab and the Cpp.

## Matlab 
The workflow for the matlab implementation is the following:
1. First define the options:
```
nsga_opt = nsgaiiOptions(n);
```
where ``` n ``` is the number of desicion variables.

The ```nsga_opt``` variable containts all the options of the Genetic Algorithm. 

2. After configuring the options create an anonymoys inline function as
```
f = @(x)testObjectiveFunctions(x, 'KUR');
```
3. Then solve the multiobjective optimization problem by running
```
ParetoSolutions = nsgaii(f, nsga_opt);
```

The repository contains the script [runTestFunctions](https://github.com/UCD-GW-Nitrate/nsgaii/blob/master/runTestFunctions.m) which show how to run several test functions that are included in the above paper.

## C++
The matlab implementation works well when the objective function is not very time consuming. However, in hydrology this is never the case, so the c++ implementation comes to the rescue. The only advantage of C++ is that distributes the objective function load over a number requested processors.
Other than that all the work is executed by only one processor, the one with 0 rank. 

1. The first step to run the c++ is to prepare a config file. To get a list of available options run
```
nsgaii.exe -h
```

2. After the configuration file is setup run:
```
C:\"Program Files\Microsoft MPI"\Bin\mpiexec.exe -n 4 .\nsgaii.exe -c D:\giorgk\Documents\GitHub\nsgaii\CPP\config.dat
```
This will distribute the objective function evaluations to 4 processors.
