function get_number_density(ρ)
    return  Dict("p" => ρ*X/m_p,
                  "3He" => ρ*Y_3He/m_3He,
                  "4He" => ρ*Y_4He/m_4He,
                  "7Li" => ρ*Z_7Li/m_7Li,
                  "7Be" => ρ*Z_7Be/m_7Be,
                  "14N" => ρ*Z_14N/m_14N,
                  "e"   => ρ/(2*m_u)*(1 + X)
                  )
end

function get_reactionrates(T9, ρ)
    n = get_number_density(ρ)

    reactions = ["pp", "33", "34", "e7", "17", "17_prime", "p14"]

    λs = [λpp, λ33, λ34, λe7, λ17, λ17_prime, λp14]

    iks = [("p", "p"),
           ("3He", "3He"), ("3He", "4He"),
           ("7Be", "e"), ("7Be", "p"),
           ("7Li", "p"), ("14N", "p")]

    reactionrates = Dict()
    for (r, λ, (i, k)) in zip(reactions, λs, iks)
        i == k ? δ = 1 : δ = 0
        reactionrates[r] = rik(n[i], n[k], ρ, λ(T9), δ)
    end

    return reactionrates
end



function get_δms()
    """ Produce and return a dictionairy
       of mass differences for each reaction """

    δm = Dict("pp"       => (m_p + m_p + m_p) - (m_3He),
              "33"       => ((m_3He + m_3He) - (m_4He + 2*m_p)),
              "34"       => (m_3He + m_4He) - m_7Be,
              "e7"       => (m_7Be) - m_7Li,
              "17_prime" => (m_7Li + m_p) - 2*m_4He,
              "17"       => (m_7Be + m_p) - 2*m_4He,
              "p14"      => (4*m_p) - m_4He)

    return δm
end

function rescale_reactionrates(r::Dict)
    """
    Rescales reaction rates such that
    no reaction uses more mass than is produced
    by the preceding reaction

    Input
    ------
    r: dictionairy of reaction rates

    Returns
    --------
    r with rescaled values
    """

    if r["pp"] < (2*r["33"] + r["34"])
        r["pp"] = 2*r["33"] + r["34"]
    end
    if r["34"] < (r["e7"] + r["17"])
        scaler = r["34"]/(r["e7"] + r["17"])
        r["e7"] *= scaler
        r["17"] *= scaler
    end

    r["17_prime"] = r["e7"]

    return r

end
