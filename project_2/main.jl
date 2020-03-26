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

    mStop = 0.05*m0
    LStop = 0.05*L0
    rStop = 0.05*R0

    p = 0.01 # variable step length parameter
    i = 1

    stable = []
    ∇s = (stable = [], star = [], ad = [])
    r = []

    push!(r, R0)
    push!(∇s.ad, ∇_ad)

    println("               ∂r∂m           |           ∂P∂m           |           ∂T∂m            |            ∂L∂m            |            Δm")

    stop = false
    while ! stop

        """ check convective stability """
        ∇_stable, ∇_star = get_∇s(u0, m0, ρ, ∇_ad)


        if ∇_star > ∇_ad
            ∇T = ∇_stable
            stability = 1
        else
            ∇T = ∇_star
            stability = 0
        end

        # println(u0[3])

        """ Compute step length """
        gradients = f(u0, m0, ρ, ∇_ad)
        Δms = abs.(p*u0 ./ gradients)     # pV/f
        Δm = -minimum(Δms)

        println("$(gradients[1])    |   $(gradients[2])    |    $(gradients[3])    |    $(gradients[4])    |    $(Δm)")


        """ Integrate one step and update parameters """
        u, m = step(f, u0, m0, Δm, ρ, ∇_ad)
        u0 = u
        m0 = m
        ρ = ρ_func(u0[2], u0[3])

        """ Save quantities for plotting """
        ∇_convert = u[2]/u[2]*gradients[2]
        push!(stable, stability)
        push!(∇s.stable, ∇_stable*∇_convert)
        push!(∇s.star, ∇_star*∇_convert)
        push!(r, u[1])


        if (i%100) == 0
            # mprogress = (1.0 - m0/M0)*100
            # Lprogress = (1.0 - u[4]/L0)*100
            # rprogress = (1.0 - u[1]/R0)*100

            mprogress = (m0/M0)*100
            Lprogress = (u[4]/L0)*100
            rprogress = (u[1]/R0)*100

            @printf("Step: %d | m: %.4f | L: %.4f | r: %.4f %s\n", i, mprogress, Lprogress, rprogress, "%")
            # println(u0)
            # println()
        end
        i += 1

        if (u0[1] <= rStop || u0[4] <= LStop || m0 <= mStop || i > N)
            stop = true
        end

    end # end while

    plot(r ./ R0, ∇s.stable)
    plot!(r ./ R0, ∇s.star)

end
