
The project consists of five source files

constants   :- various constant declarations, including the reading in of the opacity-file (on the last line)
functions   :- functions for the inner workings of the integration loop, including the RK4 solver, the temperature                gradient (and convective stability) function, opacacity, plus various other parameter functions
sanitycheck :- opacity and upper convection zone santity checks
plotting    :- functions for plotting all the various parameters
main        :- outer part of the integration loop, including initialization, final result output, and the production of                   models with varying initial parameters

In the models-folder you will find .jld-files with various stars I have produced, which can be used for plotting without having to run the integration loop each time, as described below.

The program is designed to be run interactively in the Julia terminal by calling the various functions. 
To start, include the main-file, which will import all the files and functions. I.e., run

julia> include("main.jl")

A model can be computed by calling the energyTransport-function, which takes the input (in this order)
- maximum number of integration iterations (default 5000)
- array with [R0, rho0, T0, (P0)], where P0 is optional, default is solar values
- name of the model, if you want to save the result as a JLD-file (similar to .npy files)

Running this will print the progress of the model every 100 steps, and output the final results. 
Example of running:

-------------------------------------------------------------------------
julia> energyTransport(5000, [1.2*Rs, 1.0*ρs, 0.9*Ts], "testing")

Initiating Stellar Model with Initial Parameters
R0 = 1.20e+00 Rsol | T0 = 9.00e-01 Tsol | ρ0 = 1.00e+00 ρsol | P0 = 9.00e-01 Psol

Input parameter out of bounds, extrapolating opacity
Input parameter out of bounds, extrapolating opacity
Input parameter out of bounds, extrapolating opacity
Step: 100 | m: 100.0000 | L: 100.0000 | r: 99.9476 %
Step: 200 | m: 100.0000 | L: 100.0000 | r: 99.8461 %
Step: 300 | m: 100.0000 | L: 100.0000 | r: 99.6735 %
....
Step: 2200 | m: 93.7348 | L: 100.0000 | r: 42.4375 %
Step: 2300 | m: 89.0317 | L: 99.9996 | r: 37.0198 %
Step: 2400 | m: 81.3757 | L: 99.9946 | r: 31.7544 %
Step: 2500 | m: 69.2824 | L: 99.9389 | r: 26.5410 %
Step: 2600 | m: 50.9122 | L: 99.3397 | r: 21.1249 %
Step: 2700 | m: 24.9764 | L: 93.1570 | r: 14.7744 %
Step: 2800 | m: 2.2259 | L: 60.4905 | r: 6.2480 %

Final parameter values:
Mass: 0.013 %
Radius: 3.008 %
Luminosity: 52.615 %
Outer Convection Zone Size: 1.709 %
Core Size: 21.745 %
-------------------------------------------------------------------------

To plot something (cross-section, gradients, etc), you can either load a pre-produced model, or produce a new one directly by giving the plot function initial parameters. I.e. you can do

julia> cross_section([Rs, rhos, Ts], nothing)

or

julia> cross_section(nothing, "sol")


where "sol" is a model which has already been produced with the same initial values. This can be done with all the plotting functions.

-------------------------------------------------------------------------

If you want to plot everything for one model, you can run plot_everything, and either give it initial values and a name, or give it just a model name if the model already exists for in a .jld file. For example

julia> plot_everything(solar_inits, "sol")

or

julia> plot_everything(nothing, "sol")

since "sol.jld" already exists in the models-folder. 

-------------------------------------------------------------------------

To plot everything from the best model I found, you can run

julia> plot_everything(nothing, "best_model")

or 

julia> plot_everything(best_inits, "best_model")

where the latter is if you want to compute everything anew. If you want to just perform the integration loop and print the final results of the best model, you can do

julia> energyTransport(5000, best_inits, nothing);

-------------------------------------------------------------------------

If you want to perform any of the sanity checks, just do 

julia> convection_zone_test()

or

julia> read_opacity_check()