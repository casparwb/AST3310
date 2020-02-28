The whole program consists of four files:

- constants.jl: this file contains all constans and parameters which are used in the calculations. This includes all the isotope masses, various physical constants, and parameters for the reactions themselves, such as the change in mass and the energy of the produced neutrinos.

- reactionrates.jl: this file contains individual functions for computing the various reaction rates (lambdas), and a function for computing the reaction rates per unit mass, r, given a lambda.

- functions.jl: here you will find all the functions used in the program. This includes a function for computing the energy production, a function for calculating the neutrino loss, performing the sanity check, producing the plot, and computing number densities and reaction rates given a temperature and density.

- main.jl: this is the program which can be run produce all the data sequentially. 

To run the program, you can either do (in a terminal with julia installed)

-> julia main.jl

note that this will be quite slow, as Julia requires precompilation of modules. The better way is to either open a REPL terminal, and run

julia> include("functions.jl")
julia> include("main.jl")

this will first import and compile everything, and then the main program will be run. The advantage is that is you want to re-run the main program, this will now happen very fast instead of having to wait to compilation every time. 
Additionally, you can now run every function interactively in the terminal, for example if you just want to see the energy produced by the different branches given the solar temperature, you can run

julia> chain_E = energyprod(T, \rho)[end]

or if you want to produce the plot

julia> plot_T(1000)

or however many grid points you want.