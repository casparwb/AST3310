include("functions.jl")

function read_opacity_check()
    """ Sanity check for computing opacity """

    """ parameters to check for """
    logT = [3.750, 3.755, 3.755, 3.755, 3.755,
            3.770, 3.780, 3.795, 3.770, 3.775,
            3.780, 3.795, 3.800]


    logR = [-6.00, -5.95, -5.80, -5.70, -5.55,
            -5.95, -5.95, -5.95, -5.80, -5.75,
            -5.70, -5.55, -5.50]

    R = 10 .^ logR
    R *= 1e3 # convert to SI units
    T = 10 .^ logT
    ρ = (R .* ((T*1e-6) .^ 3))



    trueValues = (κ_SI     = [2.84, 3.11, 2.68, 2.46, 2.12, 4.70,
                              6.25, 9.45, 4.05, 4.43, 4.94, 6.89,
                              7.69]*1e-3,
                  logκ_cgs = [-1.55, -1.51, -1.57, -1.61, -1.67,
                              -1.33, -1.20, -1.02, -1.39, -1.35,
                              -1.31, -1.16, -1.11])



    approxκ = get_opacity.(T, ρ)
    approxValues = (κ_SI = approxκ,
                    logκ_cgs = log10.(approxκ*10))



    errors_logκ_cgs = abs.((approxValues.logκ_cgs .- trueValues.logκ_cgs) ./ trueValues.logκ_cgs)

    errors_κ_SI = abs.((approxValues.κ_SI .- trueValues.κ_SI) ./ trueValues.κ_SI)

    println("Opacity Sanity Check")
    println("\n    No. | log κ [cgs] | κ [SI]")
    println("-----------------------------------")
    for i = 1:length(logT)
            @printf("    %-5d %.2f %s        %.2f %s\n", i, errors_logκ_cgs[i]*100, "%", errors_κ_SI[i]*100, "%")
    end

end



function convection_zone_test()
    """ Example 5.1 """

    T = 0.9e6
    ρ = 55.9
    R = 0.84*Rs
    M = 0.99*Ms
    L = Ls
    κ = 3.98
    ϵ = sum(collect(values(energyprod(Tsol, ρsol)[2])))
    α = 1
    P = P_func(ρ, T)
    g = gravity_acc(M, R)
    Hp = scale_height(T, g)
    lm = Hp
    ∇_ad = P/(T*ρ*Cp)
    ∇_stable = get_∇_stable(L, R, T, κ, Hp, ρ)


    """ ξ equation parameters """
    u = R, P, T, L
    U = 64*σ*T^3/(3*κ*ρ^2*Cp)*sqrt(Hp/(g))
    ∇_star, ξ = ∇T_con(u, M, κ, ρ, Hp, g, ∇_ad, ∇_stable)

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
