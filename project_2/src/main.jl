include("functions.jl")
include("constants.jl")
include("energyprod/core.jl")
include("plotting.jl")
include("sanitycheck.jl")

function energyTransport(N=5000, initial_parameters=solar_inits, modelname=nothing)
    """
    Function for modeling the stellar interior

    Input:
    ------------------------
    N:                  int, maximum number of integration steps
    initial_parameters: array of floats, array with [R0, ρ0, T0, (P0)],
                        where P0 is optional
    modelname:          string, if you want to save the model


    Output:
    ------------------------
    params: Dict, dictionairy with arrays of various parameter values at each step
    ∇s: array of floats, array with the temperature gradient at each step
    stable: array of ints, array with 1/0, whether the star was
            convectively stable or not at given mass shell/step (yes/no)
    """

    ### check if pressure is given as initial parameter ###
    ### if not, compute using the given density ###
    if length(initial_parameters) == 4
        R0, ρ, T0, P0 = initial_parameters
    elseif length(initial_parameters) == 3
        R0, ρ, T0 = initial_parameters
        P0 = P_func(ρ, T0)
    else
        throw("Please give three or four initial parameters")
    end

    ### initialization ###
    m0 = M0 = Ms # solar mass
    L0 = Ls      # solar luminosity
    u0 = [R0, P0, T0, L0]

    println("\nInitiating Stellar Model with Initial Parameters")
    @printf("R0 = %.2e Rsol | T0 = %.2e Tsol | ρ0 = %.2e ρsol | P0 = %.2e Psol\n \n", R0/Rs, T0/Ts, ρ/ρs, P0/Ps)

    ∇_ad = P0/(T0*ρ*Cp) # adiabatic temperature gradient

    ### arrays for storing various parameters of interest ###
    stable = zeros(N) # for storing when we have stability
    ∇s = zeros(N)     # for the gradients

    params = Dict("m" => zeros(N),
                  "T" => zeros(N),
                  "L" => zeros(N),
                  "P" => zeros(N),
                  "R" => zeros(N),
                  "ρ" => zeros(N),
                  "∇ad" => ∇_ad)

    p = 0.01     # varible step size parameter

    ### begin integration loop ###

    i = 1    # iteration counter
    stop = false
    while ! stop

        ### compute step length ###
        gradients = f(u0, m0, ρ, ∇_ad)
        Δms = abs.(p*(u0 ./ gradients))     # pV/f
        Δm = -minimum(Δms)

        ### check convective stability at current mass shell
        ### and compute corresponding temperature gradient ###
        ∇T, stability = get_∇(u0, m0, ρ, ∇_ad)

        ### save relevant values to arrays ###
        stable[i] = stability
        ∇s[i] = ∇T
        params["m"][i] = m0
        params["R"][i] = u0[1]
        params["P"][i] = u0[2]
        params["T"][i] = u0[3]
        params["L"][i] = u0[4]
        params["ρ"][i] = ρ


        ### Print progress ###
        if (i%100) == 0

            mprogress = (m0/M0)*100
            Lprogress = (u0[4]/L0)*100
            rprogress = (u0[1]/R0)*100

            @printf("Step: %d | m: %.4f | L: %.4f | r: %.4f %s\n", i, mprogress, Lprogress, rprogress, "%")
        end

        ### Integrate one step and update parameters ###
        u, m = step(f, u0, m0, Δm, ρ, ∇_ad)

        u0 = u
        m0 = m
        ρ = ρ_func(u0[2], u0[3])

        ### check stopping criterion ###
        if ((i+1) >= N || u0[1] < 0 || u0[4] < 0 || m0 < 0)
            stop = true
            break
        end

        i += 1

    end # end while

    if i >= N
        println("Maximum number of iterations reached!")
    else
        ### extract the values which were calculated (remove the zeros) ###
        stable      = stable[1:i-1]
        ∇s          = ∇s[1:i-1]
        params["m"] = params["m"][1:i-1]
        params["R"] = params["R"][1:i-1]
        params["P"] = params["P"][1:i-1]
        params["T"] = params["T"][1:i-1]
        params["L"] = params["L"][1:i-1]
        params["ρ"] = params["ρ"][1:i-1]
    end # end if

    ### print final results ###
    print_final_results(params, stable)

    ### save model ###
    if modelname != nothing
        savepath = "../models/"*modelname*".jld"
        save(savepath, "params", params, "stability", stable, "∇s", ∇s)
        println("\nModel saved as  $savepath !\n")
        return params, ∇s, stable
    else
        return params, ∇s, stable
    end # end if

end # end energyTransport func

function print_final_results(params, stable)
    """
    Function for printing the final output of the model
    after the integration loop.
        - params["X"][1] is the initial value.
    """
    ### find final progress ###
    mProgress = params["m"][end]/params["m"][1]*100
    RProgress = params["R"][end]/params["R"][1]*100
    LProgress = params["L"][end]/params["L"][1]*100

    ### find size of core ###
    core_edge = sum((params["L"]/params["L"][1] .>= 0.995) .== 0)
    core_size = 100*reverse(params["R"])[core_edge]/params["R"][1]

    ### find size of outer convection zone ###
    ctr = 1
    for i = 1:length(stable)-1
        s0 = stable[i]
        s = stable[i+1]
        if s0 == 1 && s == 1     # edge radiation zone, skip
            continue
        elseif s == 0            # outer convection zone, what we want
            ctr += 1
        elseif s == 1 && s0 == 0 # outer radiation zone, stop
            break
        end
    end

    conv_zone_size = 100*(1 - params["R"][ctr]/params["R"][1])

    println("\nFinal parameter values:")
    println("Mass: $(round(mProgress, digits=3)) %")
    println("Radius: $(round(RProgress, digits=3)) %")
    println("Luminosity: $(round(LProgress, digits=3)) %")
    println("Outer Convection Zone Size: $(round(conv_zone_size, digits=3)) %")
    println("Core Size: $(round(core_size, digits=3)) %")

end


function some_models()
    """ Produces and plots stars with varying initial parameter values """
    init_params = [[Rs, ρs, Ts],
                   [0.8*Rs, ρs, Ts],
                   [2*  Rs, ρs, Ts],
                   [Rs, ρs, 0.8*Ts],
                   [Rs, ρs, 2*Ts],
                   [Rs, 2*ρs, Ts],
                   [Rs, 0.8*ρs, Ts],
                   [Rs, ρs, Ts, 2*Ps],
                   [Rs, ρs, Ts, 0.8*Ps]]

    modelnames = ["sol", "low_R0", "high_R0",
                  "low_T0", "high_T0",
                  "high_rho", "low_rho",
                  "high_P0", "low_P0"]

    for (inits, modelname) in zip(init_params, modelnames)
        energyTransport(5000, inits, modelname)
    end

    for modelname in modelnames
        cross_section(nothing, modelname)
    end

end
