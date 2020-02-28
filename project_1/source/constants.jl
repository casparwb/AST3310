
""" General physics constant """

m_u   = 1.660539066e-27         # atomic mass unit [kg]
m_p   = 1.007825032241*m_u      # Hydrogen-1 mass [kg]
m_3He = 3.016029322650*m_u      # Helium-3 mass [kg]
m_4He = 4.002603254130*m_u      # Helium-4 mass [kg]
m_7Be = 7.016928720000*m_u      # Beryllium-7 mass [kg]
m_7Li = 7.016003437000*m_u      # Lithium-7 mass [kg]
m_14N = 14.00307400446*m_u      # Nitrogen-14 mass [kg]
m_e   = 9.1094e-31              # electron mass [kg]
k_B   = 1.3806e-23              # Boltzmann constant [m^2 kg s^-2 K^-1]
σ     = 5.6704e-8               # Stefan-Boltzmann constant [W m^-2 K^-4]
c     = 2.9979e8                # speed of light [m s^-1]
a     = 4*σ/c                   # radiation density constant
h     = 6.626e-34               # Planck constant [J s]
ϵ0    = 8.854e-12               # Vacuum permativity [s^4 A^2 kg^-1 m^-3]
MeV   = 1.602176e-13            # Joule per MeV
NA    = 6.02214e23              # Avogadro's number

""" Reaction specific constants """
Eνs = [0.265, 0, 0, 0.815, 0, 6.711, 0.707 + 0.997]*MeV # neutrino energies

# dictionairy with the neutrino energy and mass difference for each reaction
# data can be extracted as (example) reactions["pp"].Eν and reactions["pp"].δm
reactions = Dict("pp"       => (Eν = 0.265*MeV, δm = (3*m_p) - (m_3He)),
                 "33"       => (Eν = 0*MeV,     δm = ((m_3He + m_3He) - (m_4He + 2*m_p))),
                 "34"       => (Eν = 0*MeV,     δm = (m_3He + m_4He) - m_7Be),
                 "e7"       => (Eν = 0.815*MeV, δm = (m_7Be - m_7Li)),
                 "17_prime" => (Eν = 0*MeV,     δm = (m_7Li + m_p) - 2*m_4He),
                 "17"       => (Eν = 6.711*MeV, δm = (m_7Be + m_p) - 2*m_4He),
                 "p14"      => (Eν = 1.704*MeV, δm = 4*m_p - m_4He))

reactionNames = collect(keys(reactions)) # array with all reaction names

""" Mass fraction of each atomic species """

X = 0.7
Y_4He = 0.29
Y_3He = 1e-10
Z_7Li = 1e-7
Z_7Be = 1e-7
Z_14N = 1e-11


""" Solar parameters (core)"""

ρ = 1.62e5 # density [kg/m^3]
T = 1.57e7 # temperature [K]

""" Sanity test values """

sanities = Dict("pp" => 4.04e2,
                "33" => 8.68e-9,
                "34" => 4.86e-5,
                "e7" => 1.49e-6,
                "17_prime" => 5.29e-4,
                "17" => 1.63e-6,
                "p14" => 9.18e-8)
