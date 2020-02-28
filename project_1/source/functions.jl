include("constants.jl")
include("reactionrates.jl")
using Printf, PyPlot

function energyprod(T, ρ)
    """
    Function for calculating the energy produced
        by the PP- and CNO-fusion chains.

        Input:
        --------
        T9: temperature of the stellar core [T^9 K]
        ρ: density of the stellar core [kg/m^3]
        reactions: array of string with reaction names
        Eνs: array with neutrino energies from each reaction in J
        δm: mass change in each reaction (in the same order as reactions)
        r: reactions rates, unscaled


        Output:
        ---------
        ϵ: total energy produced
        Q: energy produced by each reaction
        E: total energy produced by each reaction
        chain_E:  total energy produced by each reaction chain
        """

        """ initialization """

        T9 = T*1e-9
        r = get_reactionrates(T9, ρ) # compute reaction rates with the given temperature and density
        r = rescale_reactionrates(r) # rescale reaction rates

        c2 = c*c        # to save CPU work

        # energy output (including neutrino energy)
        Q       = Dict([reaction => reactions[reaction].δm*c2
                                    for reaction in reactionNames])

        # energy output (excluding neutrino energy)
        Q_prime = Dict([reaction => Q[reaction] - reactions[reaction].Eν
                                    for reaction in reactionNames])

        # energy production
        E       = Dict([reaction => Q_prime[reaction]*r[reaction]
                                    for reaction in reactionNames])


        # compute and save energy produced by each reaction chain
        chain_E = Dict()
        chain_E["PPI"]   = (Q_prime["pp"] + Q_prime["33"])*r["33"]
        chain_E["PPII"]  = (Q_prime["pp"] + Q_prime["34"] + Q_prime["e7"] + Q_prime["17_prime"])*r["17_prime"]
        chain_E["PPIII"] = (Q_prime["pp"] + Q_prime["34"] + Q_prime["17"])*r["17"]
        chain_E["CNO"]   = Q_prime["p14"]*r["p14"]


        return Q, E, chain_E

end # end energyprod

function get_number_density(ρ)
    """
        Compute the number densities of elements
        given a core mass density

        Input
        --------------
        ρ: float, density of stellar core [kg/m^3]

        Output
        --------------
        Dict with keys equal to the different elements, and values
        with their number densities [m^-3]
    """

    return  Dict( "p"    => ρ*X/m_u,
                  "3He" => ρ*Y_3He/(3*m_u),#m_3He,
                  "4He" => ρ*Y_4He/(4*m_u),#m_4He,
                  "7Li" => ρ*Z_7Li/(7*m_u),#m_7Li,
                  "7Be" => ρ*Z_7Be/(7*m_u),#m_7Be,
                  "14N" => ρ*Z_14N/(14*m_u),#m_14N,
                  "e"   => ρ/(2*m_u)*(1 + X)
                  )
end

function get_reactionrates(T9, ρ)
    """
        Function for computing the reaction rates rik for
        all reactions ik.

        Input
        --------------------
        T9: float, temperature of stellar core [T^9 K]
        ρ:  float, stellar core mass density [kg/m^3]

        Output
        --------------------
        reactionrates: dict of floats, reaction rates for all reactions [kg^-1 s^-1]
    """
    n = get_number_density(ρ) # compute istope number densities

    reactions = ["pp", "33", "34", "e7", "17", "17_prime", "p14"]

    # array with all reaction rate functions
    λs = [λpp, λ33, λ34, λe7, λ17, λ17_prime, λp14]

    # reaction indices, for the Kronecker delta
    iks = [("p", "p"),
           ("3He", "3He"), ("3He", "4He"),
           ("7Be", "e"), ("7Be", "p"),
           ("7Li", "p"), ("14N", "p")]

    reactionrates = Dict() # for storing the rates for each reaction

    # loop through all reactions and their reaction rates and indices
    for (r, λ, (i, k)) in zip(reactions, λs, iks)

        i == k ? δ = 1 : δ = 0 # check if i=k, set δ=1 if true, else 0

        if r == "e7" # e7 needs ρ for the electron number density (electron capture limit)
            reactionrates[r] = rik(n[i], n[k], ρ, λ(T9, ρ), δ)
        else
            reactionrates[r] = rik(n[i], n[k], ρ, λ(T9), δ)
        end # end if
    end # end for

    return reactionrates
end # end get_reactionrates


function rescale_reactionrates(r::Dict)
    """
    Rescales reaction rates such that
    no reaction uses more mass than is produced
    by the preceding reaction

    Input
    ------
    r: dictionairy of reaction rates

    Returns
    --------
    r with rescaled values
    """

    if r["pp"] < (2*r["33"] + r["34"])
        scaler = r["pp"]/(2*r["33"] + r["34"])
        r["33"] *= scaler
        r["34"] *= scaler
    end # end if

    if r["34"] < (r["e7"] + r["17"])
        scaler = r["34"]/(r["e7"] + r["17"])
        r["e7"] *= scaler
        r["17"] *= scaler
    end # end if

    r["17_prime"] = r["e7"]

    return r
end # end rescale_reactionrates


function plot_T(N, saveplot=false)
    """ Function for plotting the energy production of the
    different reaction chains for various stellar core temperatures """


    Ts = range(4, stop=9, length=N)           # logspace temperatures to compute for
    chains = ["PPI" "PPII" "PPIII" "CNO"]        # names of all chains/branches
    energies = zeros(length(Ts), length(chains)) # array for storing resulting energies

    """ Produce Data """
    for i = 1:N
        T = 10^Ts[i]
        result = energyprod(T, ρ)[end]
        # result is a dictionairy we need to
        # unpack in the correct order
        energies[i,:] = log10.([result[chain] for chain in chains])
    end # end for

    """ Plotting """
    fig, (ax1, ax2) = subplots(nrows=1, ncols=2, figsize=(12, 6))

    ax1.plot(Ts, energies)

    ax1.axvline(x=log10(T), label="Solar Temperature")

    # plot 60% of T-values for making a close-up plot:
    n = Int(N/100*60)
    for (i, chain) in enumerate(chains)
        ax2.plot(Ts[n:end], energies[n:end,i], label=chain)
    end # end if

    ax2.axvline(x=log10(T))

    fig.legend()

    ax1.set_ylabel("log10 E [J kg^-1 s^-1]")
    ax1.set_xlabel("log10 T [K]")
    ax2.set_xlabel("log10 T [K]")

    suptitle("Stellar energy production as a function of temperature")

    if saveplot == true
        savefig("energy_temperature.png", dpi=500)
        println("Plot saved to folder.")
    end # end if

    gcf() # get current figure: for displaying plot

end # end plot_T

function sanitycheck(energies)
    """ Perform sanity check with the produced energies """

    println("Performing Sanity Check\n")

    errors = Dict([reaction => abs(sanities[reaction] - energies[reaction]*ρ)/sanities[reaction]
                               for reaction in reactionNames])

    println("\n  Sanity Check Benchmark Error\n-------------------------------")
    for key in collect(keys(errors))
        @printf("   %-5s   :  %-5.2f %s\n", key, round(errors[key]*100, digits=2), "%")
    end # end for

end # end sanitycheck

function neutrinoLoss(Q)
    """
    Compute neutrino loss from each reaction chain in percentage

    Input:
    ------------------
    Q: array of floats, energy output from each reaction [kg m^2 s^-2]

    Output:
    ------------------
    prints the neutrino loss from each branch in percentage
    """

    loss = Dict()

    loss["PPI"]   = 2*Eνs[1]/(Q["pp"] + Q["33"])
    loss["PPII"]  = (Eνs[1] + Eνs[4])/(Q["pp"] + Q["34"] + Q["e7"] + Q["17_prime"])
    loss["PPIII"] = (Eνs[1] + Eνs[6])/(Q["pp"] + Q["34"] + Q["17"])
    loss["CNO"]   = (Eνs[1] + Eνs[7])/Q["p14"]

    println("\n  Branch Neutrino Energy Loss\n-------------------------------")
    for chain in collect(keys(loss))
        Loss = round(loss[chain]*100, digits=3)
        @printf("     %5s      | %10.3f %s \n", chain, Loss, "%")
    end # end for
end # end neutrinoLoss
