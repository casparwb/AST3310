using PyCall
include("solver.jl")

function savedata(t, disturbe=false)
    """
    Function for saving the data using the fvis3 module.
    t: total simulation time [s]
    disturbe: whether to add gaussian perturbation to the temperature
    """

    ### initialise solver ###
    solver = initialize(disturbe=disturbe)

    # wrapper for step-function (fvis3 needs one that doesn't take any arguments)
    step() = hydro_solver(solver)

    # necessary in order to import python modules from current folder #
    pushfirst!(PyVector(pyimport("sys")."path"), "")

    if disturbe
        savepath = "../data/disturbed_" * "$(t)" * "s"
    else
        savepath = "../data/" * "$(t)" * "s"
    end

    ### save data using the visualization module by calling python ###
    py"""
    import fvis3
    vis = fvis3.FluidVisualiser()

    def savedata(t, solver, step, savepath):

        vis.save_data(t, step,\
                      rho=solver.ρ, u=solver.u, w=solver.w, \
                      e=solver.e, P=solver.P, T=solver.T, sim_fps=1.0, folder=savepath)

    """

    py"savedata"(t, solver, step, savepath)
end


function run_first()
    """ Function to run when importing for the first time, to precompile """
    solver = initialize()
    for i = 1:10 hydro_solver(solver) end
end

println("Would you like to run once to precompile? (y/n)")
ans = lowercase(readline(stdin))
if ans == "y" || ans == "yes"
    run_first()
    println("Done!")
else
    nothing
end
