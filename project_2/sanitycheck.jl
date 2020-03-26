include("functions.jl")

function read_opacity_check()

    logT = [3.750, 3.755, 3.755, 3.755, 3.755,
            3.770, 3.780, 3.795, 3.770, 3.775,
            3.780, 3.795, 3.800]


    logR = [-6.00, -5.95, -5.80, -5.70, -5.55,
            -5.95, -5.95, -5.95, -5.80, -5.75,
            -5.70, -5.55, -5.50]

    R = 10 .^ logR
    T = 10 .^ logT
    ρ = (R .* ((T*1e-6) .^ 3))

    n = length(logT)
    trueValues = (logκ_cgs = [-1.55, -1.51, -1.57, -1.61, -1.67,
                              -1.33, -1.20, -1.02, -1.39, -1.35,
                              -1.31, -1.16, -1.11],

                  κ_SI     = [2.84, 3.11, 2.68, 2.46, 2.12, 4.70,
                              6.25, 9.45, 4.05, 4.43, 4.94, 6.89,
                              7.69]*1e-3)


    approxκ = get_opacity.(T, ρ)
    approxValues = (logκ_cgs = approxκ,
                    κ_SI = ((10 .^ approxκ)/10))



    errors_logκ_cgs = abs.((approxValues.logκ_cgs .- trueValues.logκ_cgs) ./ trueValues.logκ_cgs)

    errors_κ_SI = abs.((approxValues.κ_SI .- trueValues.κ_SI) ./ trueValues.κ_SI)

    println("Opacity Sanity Check")
    println("\n    No. | log κ [cgs] | κ [SI]")
    println("-----------------------------------")
    for i = 1:n
            @printf("    %-5d %.2f %s        %.2f %s\n", i, errors_logκ_cgs[i]*100, "%", errors_κ_SI[i]*100, "%")
    end

end



function convection_zone_test()
    """ Example 5.1 """

    T = 0.9e6
    ρ = 55.9
    R = 0.84*R0
    M = 0.99*M0
    κ = 3.98
    ϵ = sum(collect(values(energyprod(Tsol, ρsol)[2])))
    α = 1
    P = P_func(ρ, T)
    g = gravity_acc(M, R)
    Hp = pressure_scale_height(T, g)
    ∇_ad = P/(T*ρ*Cp)
    # ∇_stable = -3*κ*L0/(256*π^2*σ*R^4*T^3)
    ∇_stable = 3*κ*ρ*Hp*L0/(64*π*R^2*σ*T^4)

    # L = L0 - ϵ*M
    #
    # println(-3*κ*L0)
    # println(256*π^2*σ*R^4*T^3)

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

    v = sqrt(g*lm^2/(4*Hp))*ξ



    trueValues = Dict("Hp" => 32.5e6,
                      "U"  => 5.94e5,
                      "ξ"  => 1.175e-3,
                      "∇*" => 0.400,
                      "v" => 65.62,
                      "∇_stable" => 3.26,
                      "∇_ad" => 2/5)

    approxValues = Dict("Hp" => Hp,
                      "U"  => U,
                      "ξ"  => ξ,
                      "∇*" => ∇_star,
                      "v" => v,
                      "∇_stable" => ∇_stable,
                      "∇_ad" => ∇_ad)

    println("  Quantity   | True Value | Numerical Value | Error\n---------------------------------------------------")
    for key in collect(keys(trueValues))
        error = abs((trueValues[key] - approxValues[key])/trueValues[key])
        @printf("  %-10s | %-10.2e | %-10.2e      | %.2f %s\n", key, trueValues[key], approxValues[key], error*100, "%")
    end


end
