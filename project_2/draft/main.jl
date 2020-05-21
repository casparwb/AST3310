include("functions.jl")
include("constants.jl")
include("energyprod/core.jl")
include("plotting.jl")

function energyTransport(N=5000, initial_parameters=solar_inits, modelname="nothing")

    """
    Idea:
    u: array with [r(m), P(m), T(m), L(m)]
    initial_parameters: array with [R0, ρ0, T0, (P0)] (P0 is optional)
    """

    if length(initial_parameters) == 4
        R0, ρ, T0, P0 = initial_parameters
    elseif length(initial_parameters) == 3
        R0, ρ, T0 = initial_parameters
        P0 = P_func(ρ, T0)
    else
        throw("Please give three or four initial parameters")
    end

    m0 = M0 = Ms
    L0 = Ls
    u0 = [R0, P0, T0, L0]
    println("\nInitiating Stellar Model with Initial Parameters")
    @printf("R0 = %.2e Rsol | T0 = %.2e Tsol | ρ0 = %.2e ρsol | P0 = %.2e Psol\n \n", R0/Rs, T0/Ts, ρ/ρs, P0/Ps)

    ∇_ad = P0/(T0*ρ*Cp)


    p = 0.01 # variable step length parameter
    i = 1    # iteration counter

    ### arrays for storing various parameters of interest ###
    stable = zeros(N)
    ∇s = zeros(N)

    params = Dict("m" => zeros(N),
                  "T" => zeros(N),
                  "L" => zeros(N),
                  "P" => zeros(N),
                  "R" => zeros(N),
                  "ρ" => zeros(N),
                  "∇ad" => ∇_ad)

    p = 0.01
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
        save_values(i, params, u0, m0, ρ)


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


        if ((i+1) >= N || u0[1] < 0 || u0[4] < 0 || m0 < 0)
            stop = true

            mprogress = (m0/M0)*100
            Lprogress = (u[4]/L0)*100
            rprogress = (u[1]/R0)*100

            @printf("Step: %d | m: %.4f | L: %.4f | r: %.4f %s\n", i, mprogress, Lprogress, rprogress, "%")
            break
        end

        i += 1

    end # end while

    if i >= N
        println("Maximum number of iterations reached!")
    else
        ### extract the values which were caulculated  ###
        stable = stable[1:i-1]
        ∇s = ∇s[1:i-1]
        extract_values(i-1, params)
    end # end if

    ### print final results ###
    print_final_results(params, stable)

    ### save model ###
    if modelname != nothing
        savepath = "models/"*modelname*".jld"
        save(savepath, "params", params, "stability", stable, "∇s", ∇s)
        return params, ∇s, stable
    else
        return params, ∇s, stable
    end # end if

end # end energyTransport func

function print_final_results(params, stable)
    mProgress = params["m"][end]/params["m"][1]*100
    RProgress = params["R"][end]/params["R"][1]*100
    LProgress = params["L"][end]/params["L"][1]*100

    core_edge = sum((params["L"]/params["L"][1] .>= 0.995) .== 0)
    core_size = 100*reverse(params["R"])[core_edge]/params["R"][1]

    ctr = 1
    for i = 1:length(stable)-1
        s0 = stable[i]
        s = stable[i+1]
        if s0 == 1 && s == 1 # edge radiation zone, skip
            continue
        elseif s == 0        # outer convection zone, what we want
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

function save_values(i, params, u, m, ρ)

    params["m"][i] = m
    params["R"][i] = u[1]
    params["P"][i] = u[2]
    params["T"][i] = u[3]
    params["L"][i] = u[4]
    params["ρ"][i] = ρ

end

function extract_values(imax, params)
    params["m"] = params["m"][1:imax]
    params["R"] = params["R"][1:imax]
    params["P"] = params["P"][1:imax]
    params["T"] = params["T"][1:imax]
    params["L"] = params["L"][1:imax]
    params["ρ"] = params["ρ"][1:imax]
end

function some_models()
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

    # # R0 = 1.392e+09 | P0 = 1.551e+04 | T0 = 5.770e+03 | L0 = 3.846e+26
end
