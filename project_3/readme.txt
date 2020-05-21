########################
Short summary of this readme:
------------------------------
To run the whole program:

In the REPL terminal (in the src directory) run 

julia> include("main.jl")
julia> savedata(300, true)

To animate:

run animate.py

######################################################

The src folder contains the code, and includes the files

- solver.jl: the hydro-solver, including initialization, difference schemes, time step and hydro_solver
- main.jk: the main program, that you can call to save the data using fvis3
- animate.py: a python program for animating the saved data. 
- fvis3.py: the visualization module

The data folder contains the results from the experiments, and this is where experiments are saved. The folder
with the main simulation (300 seconds with perturbation) is called "disturbed_300s", and also contains videos
and snapshots of all parameters. The "60s" and "100s" contains animations of the temperature without any added perturbation for 60 and 100 seconds.

To run the module, open the Julia REPL terminal in the src folder, and type

julia> include("main.jl")

You will be prompted whether to run a small simulation in order to precompile the program.
This is recommended as it will speed up the main simulations.

You can now call the function "savedata" with arguments simulation time in seconds, and true/false whether to add
a perturbation to the initial temperature. For example

julia> savedata(100, true)

will run a simulation for 100 seconds with added perturbation. The data will be saved in data/disturbed_100s.



You can also run the module interatively after including the main file. By running

julia> solver = initialize(disturbe=true/false);

you can initialize the solver and access its variables as you would with a normal class in python. You can
then compute the hydrodynamic equations, the time step, boundary conditions, etc. You can also check if the system
is in hydrostatic equilibrium by running

julia> solver = equilibrium(solver);

If you want to perform one time step, you simply do

julia> hydro_solver(solver)

And the updated variables can be accessed with (for example)

julia> solver.w

To reset the solver, simply re-initialize it

julia> solver = initialize(disturbe=true/false);



To animate data, run the python file animate.py. You will be asked from what folder you would
like to animate, which parameter, and whether to save the animation.