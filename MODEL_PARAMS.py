# Model parameters for heat equation with reaction

# Thermal properties
rho = 1100  # density [kg/m³]
cp = 1900  # heat capacity [J/(kg K)]
k = 0.2  # thermal conductivity [W/(m K)]
Hr = 2.5e5  # total reaction heat [J/kg]
h = 200  # heat transfer coefficients [W/(m² K)]

# Chemical kinetics parameters
A = 5e5  # pre-exponential factor [s⁻¹]
E = 8.0e4  # activation energy [J/mol]
R = 8.314  # gas constant [J/(mol K)]

# Reaction parameters
m = 1.2  # reaction exponent
n = 1.8  # reaction exponent
Ts = 393  # threshold temperature [K]

# Task parameters
L = 0.02 # [m] Domain size
T0 = 300 # [K] Initial temperature
Tm = 443 # [K] Mold temperature

# Constants
c = Hr/(cp*(Tm-T0))
tdiff = L**2 * rho * cp / k 
gamma = h*L/k