
""" General physics constant """

m_u   = 1.660539066e-27            # atomic mass unit [kg]
m_p = 1.007825032241*m_u
m_3He = 3.0160293*m_u
m_4He = 4.002602*m_u
m_8Be = 8.00530510*m_u      # Beryllium-8 mass
m_7Be = 7.01692872*m_u      # Beryllium-8 mass
m_7Li = 7.016003437*m_u     # Lithium-7 mass
m_8B  = 8.0246073*m_u        # Boron-8 mass
m_15O = 15.0030656*m_u      # Oxygen-15 mass
m_12C = 12.0000000*m_u              # Carbon-12 mass
m_13C = 13.003354*m_u       # Carbon-13 mass
m_13N = 13.00573861*m_u     # Nitrogen-13 mass
m_14N = 14.00307400446*m_u  # Nitrogen-14 mass
m_15N = 15.0001088989*m_u   # Nitrogen-15 mass
m_e   = 9.1094e-31            # electron mass [kg]
k_B = 1.3806e-23            # Boltzmann constant [m^2 kg s^-2 K^-1]
σ = 5.6704e-8               # Stefan-Boltzmann constant [W m^-2 K^-4]
c = 2.9979e8                # speed of light [m s^-1]
a = 4*σ/c                   # radiation density constant
h = 6.626e-34               # Planck constant [J s]
ϵ0 = 8.854e-12              # Vacuum permativity [s^4 A^2 kg^-1 m^-3]
MeV = 1.602176e-13          # Joule per MeV
NA = 6.02214e23             # Avogadro's number


""" Mass fraction of each atomic species """

X = 0.7
Y_4He = 0.29
Y_3He = 1e-10
Z_7Li = 1e-7
Z_7Be = 1e-7
Z_14N = 1e-11

"""X and Y are the fractional abundances by weight of hydrogen and
helium respectively """

""" Solar parameters """

# Solar core

ρ = 1.62e5 # density [kg/m^3]
T = 1.57e7 # temperature [K]
P = 3.45e16 # pressure [Pa]


""" Sanity test values """

sanities = Dict("pp" => 4.04e2,
                "33" => 8.68e-9,
                "34" => 4.86e-5,
                "e7" => 1.49e-6,
                "17_prime" => 5.29e-4,
                "17" => 1.63e-6,
                "p14" => 9.18e-8)
