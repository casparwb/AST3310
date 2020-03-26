include("functions.jl")
include("constants.jl")
include("energyprod/core.jl")

function study_type_stability()
    """
    Idea:
    u: array with [r(m), P(m), T(m), L(m)]
    """

    m0 = M0
    ρ = ρ0
    P0 = P_func(ρ0, T0)
    ϵ = sum(collect(values(energyprod(Tsol, ρsol)[end])))

    u0 = [R0, P0, T0, L0]

    println(u0)

    mStop = 0.05*m0
    LStop = 0.05*L0
    rStop = 0.05*R0

    p = 0.01 # variable step length parameter
    i = 1

    stable = []
    ∇s = (stable = [], star = [], ad = [])



    """ check convective stability """
    ∇_stable, ∇_ad, ∇_star = get_∇s(u0, m0, ρ, ϵ)
    # println(∇_stable)
    # println(∇_ad)
    # println(∇_star)

    if ∇_star > ∇_ad
        ∇T = ∇_stable
        stability = 1
    else
        ∇T = ∇_star
        stability = 0
    end

    """ Save gradients (for plotting) """
    push!(stable, stability)
    push!(∇s.stable, ∇_stable)
    push!(∇s.star, ∇_star)
    push!(∇s.ad, ∇_ad)


    """ Compute step length """
    gradients = f(u0, m0, ρ, ϵ, ∇T)
    Δms = abs.(p*u0 ./ gradients)     # pV/f
    Δm = minimum(Δms)


    """ Integrate one step and update parameters """
    @code_warntype step(f, u0, m0, Δm, ϵ, ρ, ∇T)
    # u, m = step(f, u0, m0, Δm, ϵ, ρ, ∇T)
    # u0 = u
    # m0 = m
    # ρ = ρ_func(u0[2], u0[3])


end
