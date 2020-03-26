
k = 1e-6/NA # constant to convert from cm^3 to m^3 and divide by NA (Avogadros constant)

function λpp(T9)

    returnvalue = 4.01e-15*T9^(-2/3)*exp(-3.380*T9^(-1/3))
                 *(1 + 0.123*T9^(1/3) + 1.09*T9^(2/3) + 0.938*T9)

    return returnvalue*k

end

function λ33(T9)
    returnvalue = 6.05*1e10*T9^(-2/3)*exp(-12.276*T9^(-1/3))
                 *(1 + 0.034*T9^(1/3) - 0.522*T9^(2/3) - 0.124*T9
                 + 0.353*T9^(4/3) + 0.213*T9^(5/3))

    return returnvalue*k
end

function λ34(T9)

    T9_ = T9/(1 + 4.95*1e-2*T9)

    returnvalue = 5.61e6*T9_^(5/6)*T9^(-3/2)*exp(-12.826*T9_^(-1/3))

    return returnvalue*k
end

function λe7(T9, ρ)
    if T9 < 1e-3 # include upper limit of 7Be electron capture
        ne = get_number_density(ρ)["e"]
        returnvalue = 1.57e-7/ne
    else
        returnvalue =  1.34e-10*T9^(-1/2)*(1 - 0.537*T9^(1/3) + 3.86*T9^(2/3) +
                       0.0027*T9^(-1)*exp(2.515*10^-3*T9^-1))
   end

    return returnvalue*k

end

function λ17_prime(T9)
    T9_ = T9/(1 + 0.759*T9)
    returnvalue = 1.096e9*T9^(-2/3)*exp(-8.472*T9^(-1/3))
                - 4.830e8*T9_^(5/6)*T9^(-3/2)*exp(-8.472*T9_^(-1/3))
                + 1.06e10*T9^(-3/2)*exp(-30.442*T9^(-1))

    return returnvalue*k
end

function λ17(T9)

    returnvalue =  3.11e5*T9^(-2/3)*exp(-10.262*T9^(-1/3)) +
                   2.53e3*T9^(-3/2)*exp(-7.306*T9^(-1))

    return returnvalue*k
end

function λp14(T9)
    returnvalue = 4.90e7*T9^(-2/3)*exp(-15.228*T9^(-1/3) - 0.092*T9^2)
               *(1 + 0.027*T9^(1/3) - 0.778*T9^(2/3) - 0.149*T9
               + 0.261*T9^(4/3) + 0.127*T9^(5/3))
               + 2.37e3*T9^(-3/2)*exp(-3.011*T9^(-1))
               + 2.19e4*exp(-12.53*T9^(-1))

    return returnvalue*k
end


function rik(ni, nk, ρ, λ, δik=0)
    """
    Compute reaction rate (not per unit mass)

    Input
    ----------------
    ni, nk: number density of element i and k
    ρ: density of stellar core
    λ: reaction rate per unit mass for reaction between elements i and k
    δik: dirac delta (1 if i = k (ni = nk))

    Output
    -------------
    Reaction rate (not per unit mass)
    """
    r = ni*nk*λ/(ρ*(1 + δik))
    return r
end
