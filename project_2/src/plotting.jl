
function cross_section(inits=nothing, modelname=nothing)
    """ Plots a cross section of the star model """

    if inits == nothing && modelname != nothing
        loadpath = "../models/" * modelname * ".jld"
        params, stable = load(loadpath, "params", "stability")
    elseif inits != nothing
        params, ∇s, stable = energyTransport(5000, inits, modelname)
    else
        throw("Please give model name or initial parameters")
    end

    ### initial values in solar units ###
    inits = Dict([key => params[key][1] for key in collect(keys(params))])
    R = inits["R"]/Rs
    T = inits["T"]/Ts
    ρ = inits["ρ"]/ρs
    P = inits["P"]/Ps

    core_L = 0.995*params["L"][1]
    r = params["R"] ./ params["R"][1]

    rmax = 1.2*r

    fig = figure(figsize=(12, 12))
    ax = fig.gca()

    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect(:equal)
    ax.set_xlabel(L"R/R_0", fontsize=25)
    ax.set_ylabel(L"R/R_0", fontsize=25)


    init_string = "R = $(round(R, digits=2)) R⊙, ρ = $(round(ρ, digits=2)) ρ⊙, T = $(round(T, digits=2)) T⊙, P = $(round(P, digits=2)) P⊙"
    ax.set_title(init_string, fontsize=30)

    ### plot every 20th step ###
    idxs = collect(1:20:length(stable))
    ss = stable[idxs]

    for (i, s) in zip(idxs, ss)
        if params["L"][i] >= core_L # outside core
            if s == 1 # radiation
                radius = r[i]
                circle_yellow = matplotlib.pyplot.Circle((0, 0), radius,
                             color=:yellow, fill=false)
                ax.add_artist(circle_yellow)
            elseif s == 0 # convection
                radius = r[i]
                circle_red = matplotlib.pyplot.Circle((0, 0), radius,
                             color=:red, fill=false)
                ax.add_artist(circle_red)
            end
        else # inside core
            if s == 1 # radiation
                radius = r[i]
                circle_cyan = matplotlib.pyplot.Circle((0, 0), radius,
                             color=:cyan, fill=false)
                ax.add_artist(circle_cyan)
            elseif s == 0 # convection
                radius = r[i]
                circle_navy = matplotlib.pyplot.Circle((0, 0), radius,
                             color=:navy, fill=false)
                ax.add_artist(circle_navy)
            end
        end

    end

    # create legends
    circle_red    = matplotlib.pyplot.Circle((2*rmax, 2*rmax), 0.1*rmax, color=:red, fill=true)
    circle_yellow = matplotlib.pyplot.Circle((2*rmax, 2*rmax), 0.1*rmax, color=:yellow, fill=true)
    circle_blue   = matplotlib.pyplot.Circle((2*rmax, 2*rmax), 0.1*rmax, color=:blue, fill=true)
    circle_cyan   = matplotlib.pyplot.Circle((2*rmax, 2*rmax), 0.1*rmax, color=:cyan, fill=true)

    ax.legend([circle_red,circle_yellow,circle_cyan,circle_blue],
              ["Convection outside core","Radiation outside core","Radiation inside core","Convection inside core"]
              , fontsize=30)

    if modelname != nothing
      savepath = "../figures/cross_section_" * modelname *".png"
      savefig(savepath)
  end
    return fig

end


function plot_nablas(inits=nothing, modelname=nothing)
    """ Plot the temperature gradients """

    if inits == nothing && modelname != nothing
        loadpath = "../models/" * modelname * ".jld"
        params, stable, ∇s = load(loadpath, "params", "stability", "∇s")
    elseif inits != nothing
        params, ∇s, stable = energyTransport(5000, inits, modelname)
    else
        throw("Please give model name or initial parameters")
    end

    ∇_ad = params["∇ad"]

    r = params["R"] / params["R"][1]

    ∇_star = zeros(length(r))
    ∇_stable = zeros(length(r))

    for i = 1:length(r)
        T = params["T"][i]
        R = params["R"][i]
        P = params["P"][i]
        L = params["L"][i]
        ρ = params["ρ"][i]
        m = params["m"][i]

        u = R, P, T, L

        κ = get_opacity(T, ρ)
        g = gravity_acc(m, R)
        Hp = scale_height(T, g)

        ∇_stable[i] = get_∇_stable(L, R, T, κ, Hp, ρ)
        ∇_star[i] = ∇T_con(u, m, κ, ρ, Hp, g, ∇_ad, ∇_stable[i])[1]

    end


    fig = figure(figsize=(12, 8))
    plot(r, log10.(∇_stable), label=L"\nabla_{{stable}}")
    plot(r, log10.(∇_star), label=L"\nabla^*")
    plot([r[end], r[1]], [log10(∇_ad), log10(∇_ad)], label=L"\nabla_{{ad}}")

    xlabel(L"R / R_0", fontsize=20)
    ylabel(L"\log_{{10}}{{\nabla T}}", fontsize=20)
    title("Temperature gradients (for p = 0.01)", fontsize=25)
    legend(fontsize=30,loc="lower center")

    if modelname != nothing
        savepath = "../figures/nablas_" * modelname * ".png"
        savefig(savepath)
    end
    return fig

end

function plot_params(inits=nothing, modelname=nothing)
    """ Plot the parameters as function of radius """

    if inits == nothing && modelname != nothing
        loadpath = "../models/" * modelname * ".jld"
        params, stable, ∇s = load(loadpath, "params", "stability", "∇s")
    elseif inits != nothing
        params, ∇s, stable = energyTransport(5000, inits, modelname)
    else
        throw("Please give model name or initial parameters")
    end


    R = params["R"] ./ params["R"][1]

    fig = figure(figsize=(16, 12))

    subplot(321)
    plot(R, params["T"]/params["T"][1])
    ylabel(L"T/T_0", fontsize=20)

    subplot(322)
    plot(R, params["m"]/params["m"][1])
    ylabel(L"M/M_0", fontsize=20)

    subplot(323)
    plot(R, params["L"]/params["L"][1])
    ylabel(L"L/L_0", fontsize=20)

    subplot(324)
    plot(R, log10.(params["ρ"]))
    ylabel(L"\log_{{10}}{ρ} \quad [kg/m^3]", fontsize=20)
    xlabel(L"R/R_0", fontsize=20)

    subplot(325)
    plot(R, log10.(params["P"]))
    ylabel(L"log_{{10}}{P} \quad [Pa]", fontsize=20)
    xlabel(L"R/R_0", fontsize=20)

    if modelname != nothing
        savepath = "../figures/params_" * modelname * ".png"
        savefig(savepath)
    end

    return fig
end


function plot_flux(inits=nothing, modelname=nothing)
    """ Plot the relative energy flux from radiation and convection """

    if inits == nothing && modelname != nothing
        loadpath = "../models/" * modelname * ".jld"
        params, stable, ∇s = load(loadpath, "params", "stability", "∇s")
    elseif inits != nothing
        params, ∇s, stable = energyTransport(5000, inits, modelname)
    else
        throw("Please give model name or initial parameters")
    end

    fig = figure(figsize=(12, 8))

    R = params["R"]
    L = params["L"]
    T = params["T"]
    ρ = params["ρ"]
    m = params["m"]
    P = params["P"]
    ∇_ad = params["∇ad"]

    Fcon = zeros(length(R))
    Frad = zeros(length(R))

    for i = 1:length(R)
        κ = get_opacity(T[i], ρ[i])
        g = gravity_acc(m[i], R[i])
        Hp = scale_height(T[i], g)
        lm = Hp

        u = [R[i], P[i], T[i], L[i]]

        ∇_stable = get_∇_stable(L[i], R[i], T[i], κ, Hp, ρ[i])
        ∇_star, ξ = ∇T_con(u, m[i], κ, ρ[i], Hp, g, ∇_ad, ∇_stable)

        if ∇_stable < ∇_ad
            Fcon[i] = 0
            Frad[i] = (16*σ*T[i]^4)/(3*κ*ρ[i]*Hp)*∇_stable
        else
            Fcon[i] = ρ[i]*Cp*T[i]*sqrt(g)*Hp^(-3/2)*(lm/2)^2*ξ^3
            Frad[i] = (16*σ*T[i]^4)/(3*κ*ρ[i]*Hp)*∇_star
        end

    end

    Ftot = Fcon .+ Frad
    Frad_rel = Frad ./ Ftot
    Fcon_rel = Fcon ./ Ftot

    R /= params["R"][1]
    plot(R, Frad_rel, label="Radiation")
    plot(R, Fcon_rel, label="Convection")
    legend(loc=:center, fontsize=30)

    xlabel(L"R/R_0", fontsize=25)
    ylabel(L"F/F_{{tot}}", fontsize=25)

    if modelname != nothing
        savepath = "../figures/flux_" * modelname * ".png"
        savefig(savepath)
    end

    return fig

end


function plot_energyproduction(inits=nothing, modelname=nothing)
    """ Plot the energy production as function of radius """

    if inits == nothing && modelname != nothing
        loadpath = "../models/" * modelname * ".jld"
        params, stable, ∇s = load(loadpath, "params", "stability", "∇s")
    elseif inits != nothing
        params, ∇s, stable = energyTransport(5000, inits, modelname)
    else
        throw("Please give model name or initial parameters")
    end

    R = params["R"]
    T = params["T"]
    ρs = params["ρ"]

    chains = ["PPI" "PPII" "PPIII" "CNO"]        # names of all chains/branches
    energies = zeros(length(R), length(chains)) # array for storing resulting energies
    energies_rel = zeros(length(R), length(chains)) # for ϵ(r)/ϵmax
    eps = zeros(length(R))

    """ Produce Data """
    for i = 1:length(R)
        result = energyprod(T[i], ρs[i])
        ### result is a dictionairy we need to ###
        ### unpack in the correct order ###
        energies[i,:] = [result[3][chain] for chain in chains]
        energies_rel[i,:] = energies[i,:]
        energies[i,:] = energies[i,:]/sum(energies[i,:])
        eps[i] = sum(collect(values(result[2])))
    end # end for

    ϵmax = maximum(eps)
    energies_rel /= ϵmax

    """ Plotting """
    fig = figure(figsize=(12, 6))
    ax = fig.gca()
    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

    colors = [:blue, :orange, :green, :red]

    R /= params["R"][1]
    for i = 1:4
        ax.plot(R, energies[:,i], label=chains[i], color=colors[i])
        ax2.plot(R, energies_rel[:,i], color=colors[i], linestyle="--")
    end

    """ Plot ϵ/ϵmax with separate y-axis """
    color = :purple
    ax2.set_ylabel(L"\epsilon/\epsilon_{{max}}", fontsize=20)
    #ax2.plot(R, eps, linestyle="--")
    #ax2.tick_params(axis="y")

    fig.legend(fontsize=25, loc=:center)

    ax.set_ylabel(L"\epsilon/\epsilon_{{tot}}", fontsize=25)
    ax.set_xlabel(L"R/R_0", fontsize=25)

    if modelname != nothing
        savepath = "../figures/energyprod_" * modelname * ".png"
        savefig(savepath)
    end

    return fig
end

function plot_everything(inits, modelname)
    """ Plot everything for one model """

    cross_section(inits, modelname)
    plot_nablas(nothing, modelname)
    plot_flux(nothing, modelname)
    plot_energyproduction(nothing, modelname)
    plot_params(nothing, modelname)

end
