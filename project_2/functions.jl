using Interpolations, Printf, Polynomials, Plots, PlotThemes
pyplot()
theme(:juno)

function f(u, m, ρ, ∇_ad)
    """
    Function for computing gradients for r(m), P(m), T(m), L(m).
     """

    r, P, T, L = u

    ∇_stable, ∇_star = get_∇s(u, m, ρ, ∇_ad)
    if ∇_star > ∇_ad
        ∇T = ∇_stable
    else
        ∇T = ∇_star
    end

    # ∇_star > ∇_ad ? ∇T = ∇_stable : ∇T = ∇_star

    ϵ = sum(collect(values(energyprod(u[3], ρ)[2])))

    """ gradients """
    ∂r∂m = 1.0/(4*π*r^2*ρ)
    ∂P∂m = -G*m/(4*π*r^4)
    ∂L∂m = ϵ
    ∂T∂m = ∇T*T/P*∂P∂m

    return [∂r∂m, ∂P∂m, ∂T∂m, ∂L∂m]

end

pressure_scale_height(T, g) = k_B*T/(μ*m_u*g)
geometric_factor(lm) = 4.0/lm^2
gravity_acc(m, r) = G*m/r^2

function get_∇s(u, m, ρ, ∇_ad)

    r, P, T, L = u
    κ = get_opacity(T, ρ)
    κ = 10^(κ - 1) # convert from log(κ) in cgs units to κ in SI unit
    g = gravity_acc(m, r)
    Hp = pressure_scale_height(T, g)
    # ∇_stable = -3*κ*L/(256*π^2*σ*r^4*T^3)
    ∇_stable = 3*κ*ρ*Hp*L0/(64*π*r^2*σ*T^4)
    ∇_star = ∇T_con(u, m, κ, ρ, Hp, g, ∇_ad, ∇_stable)

    return ∇_stable, ∇_star

end

function step(f, u0, m0, Δm, ρ, ∇_ad)
    """ Perform one RK4 step with step length Δm"""

    K1 = Δm*f(u0,        m0,        ρ, ∇_ad)
    K2 = Δm*f(u0+0.5*K1, m0+0.5*Δm, ρ, ∇_ad)
    K3 = Δm*f(u0+0.5*K2, m0+0.5*Δm, ρ, ∇_ad)
    K4 = Δm*f(u0+K3,     m0+Δm,     ρ, ∇_ad)

    u = u0 - (1/6)*(K1 + 2*K2 + 2*K3 + K4)
    m = m0 - Δm
    return u, m
end

function ρ_func(P, T)
    return μ*m_u/(k_B*T)*(P - 4*σ*T^4/(3*c))
end

function P_func(ρ, T)
    """ P = Pgas + Prad (p)"""

    Pgas = ρ*k_B*T/(μ*m_u)
    Prad = a*T^4/3

    return Pgas + Prad
end

function ∇T_con(u, m, κ, ρ, Hp, g, ∇_ad, ∇_stable)

    r, P, T, L = u

    """ ξ equation parameters """
    δ = 1#T/V*N*k_b/P
    α = 1#V*N*k_B*T/P
    U = 64*σ*T^3/(3*κ*ρ^2*Cp)*sqrt(Hp/(g*δ))

    lm = 1.0*Hp

    Ω = geometric_factor(lm)

    a = 1
    b = U/lm^2
    c = U^2*Ω/lm^2
    d = U/lm^2*(∇_ad - ∇_stable)

    ξ = roots(Poly([d, c, b, a]))

    if any(typeof.(ξ) .== Complex{Float64})
        ξ = real(ξ[imag(ξ) .== 0.0])[1] # only get the real root, where imaginary part is 0
    else
        ξ = ξ[1] # if it does not have any complex roots
    end


    ∇_star = ∇_stable - lm^2/U*ξ^3 # ∇*


    return ∇_star

end

function read_opacity(filename="opacity.txt")

    f = open(filename)
    lines = readlines(f)                 # read all lines into array
    nrows = length(lines) - 2
    ncols = length(split(lines[1])) - 1

    logκ = zeros(Float64, nrows, ncols)
    logT = zeros(Float64, nrows)

    logR = [parse(Float16, w) for w in split(lines[1])[2:end]]

    for (row, ln) in enumerate(lines[3:end])
        words = split(ln)

        logT[row] = parse(Float64, words[1])
        logκ[row,:] = parse.(Float64, words[2:end])

    end # end for

    close(f)

    return logR, logT, logκ
end

function get_opacity(T_in, ρ_in, verbose=false, round=false)
    logR, logT, logκ = κ_table

    logR_in = log10(ρ_in/(T_in*1e-6)^3) # convert from ρ to log R
    logT_in = log10(T_in)

    # if round
    #     logR_in = round(logR_in, digits=2)
    #     logT_in = round(logT_in, digits=2)
    # end

    # check if exact value is present
    T_eq = logT_in .== logT
    R_eq = logR_in .== logR

    if any(T_eq) && any(R_eq) # if exact value is present
        rowidx = argmax(T_eq)
        colidx = argmax(R_eq)

        return logκ[rowidx, colidx]

    else                     # else interpolate

        if verbose println("Input parameters not in data file, interpolating") end
        knots = (logT, logR)
        itp = interpolate(knots, logκ, Gridded(Linear()))

        if logT_in < logT[1] || logT_in > logT[end] || logR_in < logR[1] || logR_in > logR[end]
            if verbose println("Input parameter out of bounds, extrapolating") end

            etp = extrapolate(itp, Flat())

            return etp(logT_in, logR_in)

        else
            return itp(logT_in, logR_in)
        end # end inner if
    end # end outer if

end # end function
