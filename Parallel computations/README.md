During my studying I've been given a task to make solver for three dimension Cauchy problem for Laplace's equation with different types of boundary conditions.

Main requirement of this task were to make three solvers: straight (just a program), MPI realization (using MPI and 3d separation of a volume), MPI + OpenMP realization (just like second one but also using OpenMP to increase speed).

Every node of comm_world group takes it's own part of simulation volume, information exchange by using SendRecv. 

All computations were made using University supercomputer Bluegene.

Below you can see 2d slice of finite result animated in matlab.

![Inf image](https://github.com/DMuhayev/Portfolio/blob/master/Parallel%20computations/result.gif?raw=true)
