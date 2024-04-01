# EPFL ME467: Turbulence
### Time-marching the KSE
The following program (`/programs/example1.m`) is an example of how to utilize the function `KSE_integrate` to advance the KSE forward in time.
```
close all; clear; clc

addpath('../functions/')

L = 39;
N = 64;
symm = true;
T = 500;
dt = 0.1;

[x,~] = domain(L,N);
u0 = sin(2*pi*x/L);

v0 = field2vector(u0,N,symm);
[vT,~] = KSE_integrate(v0,T,dt,0,L,N,symm);

uT = vector2field(vT,N,symm);

figure
  plot(x,uT,'LineWidth',2)
  hold on; grid on
  plot(x,u0,'LineWidth',2)
  xlabel('x'); ylabel('u')
  legend('u(T)','u(0)')
```
This program can be broken down into the following steps:
- We first tell MATLAB to look for functions in a folder of the relative path `../functions/`.
- Then, we define the domain length `L = 39`, number of grid points `N = 64`, center symmetry `symm = true`, integration time `T = 500` and time step size `dt = 0.1`.
- To define the initial condition, we first construct the grid points using `[x,~] = domain(L,N)`. Then, the initial condition is defined as the sine function `u0 = sin(2*pi*x/L)`. Note that this initial condition is consistent with the imposed center symmetry.
- We construct the state vector associated with `u0` by calling `v0 = field2vector(u0,N,symm)`. Finally, we call `[vT,~] = KSE_integrate(v0,T,dt,0,L,N,symm)` where `vT` is the state vector of the KSE obtained by advancing `v0` for time `T`.
- Since `vT` is a state vector and not the vector of grid values, we call `uT = vector2field(vT,N,symm)` to transform `vT` to its corresponding physical field `uT` and, then, plot `u0` and `v0` in one figure.

 **Note 1**: By setting the 4th argument of `KSE_integrate` to `0`, no intermediate snapsahot is stored and `vT` is a column vector. In order to store intermediate snapshots, you can call `[vt,t] = KSE_integrate(v0,T,dt,dt_store,L,N,symm)`. As a result, every `dt_store` the snapshot is appended to `vt` as a new column. The vector `t` contains the time instance assocaited with each column of `vt`.

### Search for a periodic orbit
The following program (`/programs/example2.m`) is an example of how to utilize the function `search4PO` to converge a periodic orbit from a guess. Here, the same sine function as the previous program and a period of `T_guess=50` are chosen for demonstration.
```
close all; clear; clc

addpath('../functions/')

L = 39;
N = 64;
symm = true;
dt = 0.1;

[x,~] = domain(L,N);
u0 = sin(2*pi*x/L);

u_guess = field2vector(u0,N,symm);
T_guess = 50;

[u_best,T_best] = search4PO(u_guess,T_guess,dt,L,N,symm); 
```
This crude guess does not converge, and the following will be displayed
```
                                         Norm of      First-order   Trust-region
 Iteration  Func-count     f(x)          step         optimality    radius
     0         23         5319.86                      3.37e+03               1
     1         46         4366.73              1        1.6e+03               1
   ...        ...             ...            ...            ...             ...
    75       1550         2335.65      0.0145519           22.4          0.0146

Solver stopped prematurely.

fsolve stopped because it exceeded the iteration limit,
options.MaxIterations = 7.500000e+01.

Error using search4PO (line 36)
Search for unstable periodic orbit failed...
```
In case of a successful search, `u_best` is the converged state vector and `T_best` is the converged period of the periodic orbit. In that case, if you advance `u_best` as the initial condition for the time interval of `T_best` (see the first example above), the trajectory closes on itself. Note that `u_guess` and `u_best` are both state vectors.
