# Langevin Simple

This program simulates particles immersed on a fluid under significant brownian agitation. 
The moviment of the particles is integrated from a Generalized Langevin Equation that is a form of the Newton's secound law. The novelty is the implementation of a integro-differential equation that discribes a non-Markovian GLE. 

## Compilation

The Fortran code can be compiled using mpiifort (Intel Fortran Compiler) or mpifort (GCC Fortran). You may change the make file to choose the compiler that you wish. 

## Requirements
GNU Fortran 
```
sudo apt-get install gfortran # Ubuntu Based distros
sudo zypper in gcc-fortran # OpenSUSE
```
Or [Intel Fortran Compiler](https://software.intel.com/en-us/fortran-compilers)

If using GNU Fortran, install also the OpenMPI package
```
sudo apt-get install openmpi-bin libopenmpi-dev # Ubuntu 
sudo zypper in openmpi openmpi-devel # OpenSUSE
```
For the Python code it is necessary the pyevtk package to convert the csv files to vtk. 
The package progressbar2 is optional. 

For the visualization, you can use Paraview or some program that reads VTK files. 


## Input and Output

The input needs a configuration file named settings.ini , which contains the physical properties and numerial parameters for the simulation, and files containing the positions of the particles. Optionally a file containing the initial velocity of the particles can also be provided. 

The output contains the position and velocities of the particles. 

## Postprocessing

The program automatically converts the .csv results to VTK, so one may use this result to postprocess on paraview. 
Some python scripts are provided to calculate autocorrelations of velocities and concentration distributions. The theory behind these program is summarized at the Master Thesis that will be published soon. 

## References

Some of the code on this repository was based on code published by 
Jason R. Blevins ACM Fortran Forum 28(3), 2–7, 2009.
https://github.com/jannisteunissen/config_fortran

## Acknowledgments 

The author to acknowledge the support of this work by FAPESP-Brazil under the grant number 2017/05643-8 and FAEPEX Unicamp. 


