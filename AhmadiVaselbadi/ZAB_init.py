import numpy as np


i_ext = 100  # External current (A/m2)
t_final = 100

" Chage transfer inputs "
F = 96485    # Faraday's constant, C/mol equivalent charge
R = 8.3145     # Gas constant, J/mol-K
T = 298.15 # Temperature, K
P_ca = 101325 #Ambient pressure, Pa

"Cathode"
k_f_eq_ca = 2.37e-11 # Chemical forward rate constant of cathode reaction, mol/s
k_r_eq_ca = 2.07e-11 # Chemical reverse rate constant of cathode reaction, mol/s

alpha_ca = 0.5
n_ca = 2 # Number of electron transferred
C_0_ca = P_ca/R/T 

# Composition including CO2
X_O2_ca = 0.20
X_H2O_ca = 0.1
X_N2_ca = 0.58
X_OH_ca = .1 
X_CO2_ca = 0.3
X_CO3_ca = 0


X_wCO2 = np.array([X_O2_ca, X_H2O_ca,X_N2_ca, X_OH_ca, X_CO2_ca, X_CO3_ca])
C_wCO2 = X_wCO2 * C_0_ca

#Composition with no CO2
X_O2_ca = 0.20
X_H2O_ca = 0.1
X_N2_ca = 0.78
X_OH_ca = 0.1
X_CO2_ca = 0.0001
X_CO3_ca = 0


X_noCO2 = np.array([X_O2_ca, X_H2O_ca,X_N2_ca, X_OH_ca, X_CO2_ca, X_CO3_ca])
C_noCO2 = X_noCO2 * C_0_ca

# Gas Diffusion coefficients
D_CO3 = 8.04e-10
D_H2O = 2.59e-5
D_CO2 = 1.65e-5
D_N2 = 2.798e-5
D_O2 = 2.438e-5
D_OH = 4.94e-9 

# Diffusion coeffiecient array
D_k = np.array([D_O2, D_H2O, D_N2, D_OH, D_CO2, D_CO3]) #m2/s

# Gibbs free energy data
g_H2O = -228.572e3 #J/mol
g_O2 = 0
g_OH = -157.2e3 #J/mol

# Standard Gibbs free energy of the cathode reaction
delg_0_ca = 2*g_OH - g_H2O

rho_h2O = 997.048 #kg/mol
mu_air = 2.08e-5 #kg/m-s


"Anode"
k_f_eq_an = 4.78e-11 # Chemical forward rate constant of anode reaction, mol/s
k_r_eq_an = 3.53e-8 # Chemical forward rate constant of anode reaction, mol/s

V_zn = 9.16e-6 #m3/mol
V_zno = 14.5e-6 #m3/mol
k_ZnO = 0.25e-3 #  Chemical rate costant of ZnO Formation, m3/s

alpha_an = 0.5
n_an = 2 # Number of electron transferred


"Electolyte"
# Considering a 32wt% KOH at the electrolye-GDL interface

P_sep = 10e3  # Pressure at electrolyte-GDL interface
rho_sep = 1301 #kg/m3

# Composition of air at electrolyte-GDL interface
X_H2O_sep = .25 
X_O2_sep = .45 #mol/m3
X_OH_sep = .0001 #mol/m3
X_CO3_sep = 1e-8 #mol/m3
X_CO2_sep = 0
X_N2_sep = 0.3

X_k_sep = np.array([X_O2_sep, X_H2O_sep,X_N2_sep, X_OH_sep, X_CO2_sep, X_CO3_sep])
C_k_sep = X_k_sep*P_sep/R/T


" Microstructure and geometry "

d_part_an = 75e-6 # Anode particle radius, m

n_Brugg = -0.5

d_part_ca = 0.5e-6 # Cathode particle radius, m
r_p_ca = 2e-6

H_an = 4.5e-3  # Anode thickness, m
H_sep = 0.1e-3 # Electrolyte separator thickness, m
H_ca = 0.3e-3  # Cathode thickness, m

eps_zn_an = .185 # Volume fraction of zinc in anode

eps_inact_an = 0  # Volume fraction of inactive material in anode
eps_inact_sep = .185 # Volume fraction of inactive material in separator
eps_inact_ca = .185 # Volume fraction of inactive material in cathode

eps_elyte_an = 0.515 # Volume fraction of elyte in anode
eps_elyte_ca = 0.515 # Volume fraction of elyte in cathode
eps_elyte_sep = 0.515 # Volume fraction of elyte in separator

eps_air_an = 0.3 # Volume fraction of gas phase in anode
eps_air_sep =.3 # Volume fraction of gas phase in separator
eps_air_ca = 0.3 # Volume fraction of gas phase in cathode

A_sep = 1.33e-3 #m2

" Charge transfer inputs "
# Initial voltages (used to calculate phi_dl guesses)
phi_an_0 = 0
phi_elyte_0 = 0.6
phi_ca_0 = 1.1

C_dl_an = 2500 # Anode Double layer capacitance, F/m2 of dl interface
C_dl_ca = 140 # Cathode Double layer capacitance, F/m2 of dl interface

dphi_dl_an_0 = phi_elyte_0 - phi_an_0 
dphi_dl_ca_0 = phi_ca_0 - phi_elyte_0
