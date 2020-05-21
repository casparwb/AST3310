include("functions.jl")
include("constants.jl")
include("energyprod/core.jl")

function energyTransport(N=Inf, verbose=false)

    """
    Idea:
    u: array with [r(m), P(m), T(m), L(m)]
    """

    m0 = M0
    ρ = ρ0
    P0 = P_func(ρ0, T0)

    u0 = [R0, P0, T0, L0]
    println(u0)
    println()

    ∇_ad = P0/(T0*ρ*Cp)

    # println(∇_ad) # 0.40000720723418803

    mStop = 0.01*m0
    LStop = 0.01*L0
    rStop = 0.01*R0

    p = 0.01 # variable step length parameter
    i = 1    # iteration counter

    ### arrays for storing various parameters of interest ###
    stable = zeros(10000)
    ∇s = zeros(10000)
    r = zeros(10000)

    # println("               ∂r∂m           |           ∂P∂m           |           ∂T∂m            |            ∂L∂m            |            Δm")

    stop = false
    while ! stop

        ### compute step length ###
        gradients = f(u0, m0, ρ, ∇_ad)
        Δms = abs.(p*(u0 ./ gradients))     # pV/f
        Δm = -minimum(Δms)

        ### check convective stability at current mass shell ###
        ∇T, stability = get_∇(u0, m0, ρ, ∇_ad)

        ### save relevant values to arrays ###
        stable[i] = stability
        r[i]      = u0[1]
        ∇s[i]     = ∇T

        # println("$(gradients[1])    |   $(gradients[2])    |    $(gradients[3])    |    $(gradients[4])    |    $(Δm)")

        """ Integrate one step and update parameters """
        u, m = step(f, u0, m0, Δm, ρ, ∇_ad)
        u0 = u
        m0 = m
        ρ = ρ_func(u0[2], u0[3])


        """ Print progress """
        if (i%100) == 0

            mprogress = (m0/M0)*100
            Lprogress = (u[4]/L0)*100
            rprogress = (u[1]/R0)*100

            @printf("Step: %d | m: %.4f | L: %.4f | r: %.4f %s\n", i, mprogress, Lprogress, rprogress, "%")
        end


        if (i > N || u0[1] <= rStop || u0[4] <= LStop || m0 <= mStop)
            stop = true
        end

        i += 1


    end # end while


    ### extract the values which were caulculated ###
    stable = stable[1:i]
    r      = r[1:i]
    ∇s     = ∇s[1:i]

    """ plotting """
    unstable = stable .== 0
    stable   = stable .== 1

    ∇_star   = ∇s[unstable]
    ∇_stable = ∇s[stable]

    r_rel = r ./ R0

    r_rel_star   = r_rel[unstable]
    r_rel_stable = r_rel[stable]

    plot(r_rel, log10.(∇s), label=false)
    plot!(r_rel_star, log10.(∇_star), label=L"\nabla*")
    plot!(r_rel_stable, log10.(∇_stable), label=L"\nabla_{{stable}}")
    plot!(r_rel, [log10(∇_ad), log10(∇_ad)], label=L"\nabla_{{ad}}")
    xlabel!(L"R / R_0")
    ylabel!(L"∇T")
    title!("Temperature gradients (for p = $p)")
    # savefig("figures/gradients.png")


end
