using Interpolations, Printf, Polynomials, PyPlot, JLD


function step(f, u0, m0, Δm, ρ, ∇_ad)
    """ Perform one RK4 step with step length Δm"""

    K1 = Δm*f(u0,        m0,        ρ, ∇_ad)
    K2 = Δm*f(u0+0.5*K1, m0+0.5*Δm, ρ, ∇_ad)
    K3 = Δm*f(u0+0.5*K2, m0+0.5*Δm, ρ, ∇_ad)
    K4 = Δm*f(u0+K3,     m0+Δm,     ρ, ∇_ad)

    u = u0 + (1/6)*(K1 + 2*K2 + 2*K3 + K4)
    m = m0 + Δm
    return u, m
end

""" some simple parameter functions (like lambda functions in python)"""
scale_height(T, g) = k_B*T/(μ*m_u*g)              # pressure scale height
geometric_factor(lm) = 4.0/lm                     # gas parcel geometric factor S/(Q*d)
gravity_acc(m, r) = G*m/r^2                       # gravitational acceleration at current mass shell
ρ_func(P, T) = μ*m_u/(k_B*T)*(P - 4*σ*T^4/(3*c))  # density
P_func(ρ, T) = ρ*k_B*T/(μ*m_u) + 4*σ*T^4/(3*c)    # pressure

""" stable temperature gradient function """
get_∇_stable(L, r, T, κ, Hp, ρ) = 3*κ*ρ*Hp*L/(64*π*r^2*σ*T^4)

function f(u, m, ρ, ∇_ad)
    """
    Function for computing gradients for r(m), P(m), T(m), L(m) at
    current mass shell m.

    Input
    ------------------
    u: array of floats with [r, P, T, L]
    m: float, mass of current mass shell
    ρ: float, density at current mass shell
    ∇_ad: float, adiabatic temperature gradient (constant)

    Output
    ------------------
    - array with [∂r∂m, ∂P∂m, ∂T∂m, ∂L∂m]
    """

    r, P, T, L = u

    ### compute temperature gradient ###
    ∇T = get_∇(u, m, ρ, ∇_ad)[1]

    ### energy production at current mass shell ###
    ϵ = sum(collect(values(energyprod(T, ρ)[2])))

    ### compute and return gradients ###
    ∂r∂m = 1.0/(4*π*r^2*ρ)
    ∂P∂m = -G*m/(4*π*r^4)
    ∂L∂m = ϵ
    ∂T∂m = ∇T*(T/P)*∂P∂m

    return [∂r∂m, ∂P∂m, ∂T∂m, ∂L∂m]

end


function get_∇(u, m, ρ, ∇_ad)
    """
    Function for checking convective stability at current
    mass shell m, and return the corresponding temperature
    gradient ∇T.

    Input:
    -----------------
    u: array of floats with [r, P, T, L]
    m: float, mass of current mass shell
    ρ: float, density at current mass shell
    ∇_ad: float, adiabatic temperature gradient (constant)

    Output
    ------------------
    - tuple with
        - ∇T: float, (either ∇* or ∇stable)
        - stability: bool, whether the system is convectively stable or not
    """

    r, P, T, L = u

    ### compute opacity ###
    κ = get_opacity(T, ρ)

    ### compute gravitational acceleration and pressure scale height ###
    g = gravity_acc(m, r)
    Hp = scale_height(T, g)

    ### compute temperature gradients and check convective stability ###
    ∇_stable = get_∇_stable(L, r, T, κ, Hp, ρ)
    # ∇_stable = get_∇_stable(κ, L, P, T, m)
    if ∇_stable < ∇_ad
        # stable
        ∇T = ∇_stable
        stability = 1
    else
        # unstable
        ∇T = ∇T_con(u, m, κ, ρ, Hp, g, ∇_ad, ∇_stable)[1]
        stability = 0
    end

    return ∇T, stability
end


function ∇T_con(u, m, κ, ρ, Hp, g, ∇_ad, ∇_stable)
    """ Function for computing the temperature gradient for
    convective energy transport (∇*)

    Output
    --------------------
    - tuple with
        - ∇*
        - ξ
    """


    r, P, T, L = u

    ### ξ equation parameters ###
    U = 64*σ*T^3/(3*κ*ρ^2*Cp)*sqrt(Hp/(g))
    lm = Hp
    Ω = geometric_factor(lm)

    ### aξ^3 + bξ^2 + cξ + d = 0 ###
    a = 1
    b = U/lm^2
    c = U^2*Ω/lm^3
    d = U/lm^2*(∇_ad - ∇_stable)

    ### compute roots ###
    ξ = roots(Poly([d, c, b, a]))

    ### extract real root ###
    if any(typeof.(ξ) .== Complex{Float64})
        ξ = real(ξ[imag(ξ) .== 0.0])[1]    # i.e. where imaginary part is 0
    else
        ξ = ξ[1] # if it does not have any complex roots
    end

    ### compute ∇* ###
    ∇_star = ξ^2 + U*Ω/lm*ξ + ∇_ad

    return ∇_star, ξ

end # end ∇T_con function

function read_opacity(filepath="opacity.txt")
    """ Function for reading the opacity table file and saving to arrays """

    f = open(filepath)
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
end # end read_opacity function

function get_opacity(T_in, ρ_in)
    """
    Function for read opacity from file, or interpolate/extrapolate
    if necessary.

    Input
    ---------------------
    T_in: float, temperature in [K]
    ρ_in: float, density in [kg/m^3]

    Output
    ---------------------
    - κ: opacity in SI units
    """

    logR, logT, logκ = κ_table # read values from file

    ρ_in *= 1e-3 # convert to cgs units

    logR_in = log10(ρ_in/((T_in*1e-6)^3)) # convert from ρ to log R
    logT_in = log10(T_in)


    # check if exact value is present
    T_eq = logT_in .== logT
    R_eq = logR_in .== logR


    if any(T_eq) && any(R_eq) # if exact value is present
        colidx = argmax(T_eq)
        rowidx = argmax(R_eq)
        logκ_out = logκ[colidx, rowidx]
        κ = 10 ^ (logκ_out - 1)
        return κ

    else     # else interpolate

        knots = (logT, logR)
        itp = interpolate(knots, logκ, Gridded(Linear()))

        if logT_in < logT[1] || logT_in > logT[end] || logR_in < logR[1] || logR_in > logR[end]
            println("Input parameter out of bounds, extrapolating opacity")

            etp = extrapolate(itp, Line())

            logκ_out = etp(logT_in, logR_in)
            κ = 10 ^ (logκ_out - 1)
            return κ

        else
            logκ_out = itp(logT_in, logR_in)
            κ = 10 ^ (logκ_out - 1)
            return κ
            
        end # end inner if
    end # end outer if

end # end get_opacity function
