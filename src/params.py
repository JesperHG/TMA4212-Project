"""
Physical and model parameters for the heat equation with vulcanization reaction.
"""

# Thermal properties
rho = 1100      # density [kg/m³]
cp  = 1900      # heat capacity [J/(kg K)]
k   = 0.2       # thermal conductivity [W/(m K)]
Hr  = 2.5e5     # total reaction heat [J/kg]
h   = 200       # heat transfer coefficient [W/(m² K)]

# Chemical kinetics parameters
A = 5e5         # pre-exponential factor [s⁻¹]
E = 8.0e4       # activation energy [J/mol]
R = 8.314       # gas constant [J/(mol K)]

# Reaction exponents
m = 1.2
n = 1.8

# Temperatures
Ts = 393        # threshold temperature [K]
T0 = 300        # initial temperature [K]
Tm = 443        # mold temperature [K]

# Geometry
L = 0.02        # domain half-size [m]

# Derived constants
c     = Hr / (cp * (Tm - T0))
tdiff = L**2 * rho * cp / k
gamma = h * L / k
