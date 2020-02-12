include("main.jl")

using PyPlot

function plot_T(N::Int=100, saveplot=false)
    """ Function for plotting the energy production of the
    different reaction chains for various stellar core temperatures """

    Ts = range(4, stop=9, length=N)

    labels = ["PPI", "PPII", "PPIII", "CNO"]

    energies = zeros(length(Ts), length(labels))
    for i = 1:length(Ts)
        energies[i,:] = energyprod(10^Ts[i], ρ)
    end


    fig, (ax1, ax2) = subplots(nrows=1, ncols=2, figsize=(12, 6))

    n = Int(N/100*70)
    for i = 1:length(labels)
        ax1.plot(Ts, log10.(energies[:,i]))
        ax2.plot(Ts[n:end], log10.(energies[n:end,i]), label=labels[i])
    end
    fig.legend()

    ax1.set_ylabel("log10 ϵ [J kg^-1 s^-1]")
    ax1.set_xlabel("log10 T [K]")
    ax2.set_xlabel("log10 T [K]")

    suptitle("Stellar energy production as a function of temperature")

    if saveplot savefig("energy_temperature.png", dpi=500) end
    gcf()


end

plot_T(100, true)
