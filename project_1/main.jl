include("constants.jl")
include("reactionrates.jl")
include("functions.jl")



function energyprod(T, ρ, sanitycheck=false)
    """
    Function for calculating the energy produced
        by the PP- and CNO-fusion chains.

        Input:
        --------
        T: temperature of the stellar core [K]
        ρ: density of the stellar core [kg/m^3]
        sanitycheck: true/false, whether to perform the sanity checks and print errors

        Output:
        ---------
        prints the total energy production and returns
        the energy produced by each of the reaction chains
        """

        """ initialization """
        reactions = ["pp", "33", "34", "e7", "17_prime", "17", "p14"]
        Eνs = [0.265, 0, 0, 0.815, 0, 6.711, 0.707 + 0.997] * MeV # neutrino energies

        # misc
        T9 = T*1e-9     # Convert temperature in K to K^9
        c2 = c*c        # to save CPU work
        ϵ = 0           # initialize total energy gain

        # precalculate reaction rates and mass differences
        r = get_reactionrates(T9, ρ)
        δm = get_δms()

        # rescale reaction rates
        r = rescale_reactionrates(r)

        # array for storing energy production from each reaction
        energies = []
        Qs = []
        neutrino_loss = zeros(4)
        chain_energies = zeros(4)

        println("\nStarting reaction chains with T = $(round(T, digits=2)) K.")
        for (i, reaction) in enumerate(reactions)
            Q = δm[reaction]*c2 - Eνs[i]

            δE = Q*r[reaction]
            ϵ += δE

            push!(energies, δE)
            push!(Qs, Q + Eνs[i])

        end # end for
        println("Reaction chains complete.")
        println("Total energy produced: $(round(ϵ, digits=4)) J s^-1 kg^-1.")

        chain_energies[1] = energies[1] + energies[2]
        chain_energies[2] = energies[1] + sum(energies[3:5])
        chain_energies[3] = energies[1] + energies[6]
        chain_energies[4] = energies[7]

        neutrino_loss[1] = Eνs[1]/(Qs[1] + Qs[2])
        neutrino_loss[2] = (Eνs[1] + Eνs[4])/(Qs[1] + sum(Qs[3:5]))
        neutrino_loss[3] = Eνs[6]/(Qs[1] + Qs[6])
        neutrino_loss[4] = Eνs[7]/Qs[7]

        neutrino_loss *= 100

        """ Sanity check """
        if sanitycheck
            println("Performing sanity check.")

            errors = [abs(sanities[reaction] - energies[i]*ρ)/sanities[reaction]
                      for (i, reaction) in enumerate(reactions)]

            for (i, error) in enumerate(errors)
                println("Sanity check $i error: $(round(error*100, digits=2)) %")
            end # end for
        end # end if

        return chain_energies

    end # end energyprod
