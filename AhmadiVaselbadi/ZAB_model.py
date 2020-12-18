"""
    This file runs and executes a youyr model, calculating the cell/device/system properties of interest as a function of time. 

    The code structure at this level is meant to be very simple, and delegates most functions to lower-level modules, which can be updated with new capabilties, over time.  The code:

        1 - Reads inputs and initializes the model
        
        2 - Calls the residual function and integrates over the user-defined    
            time span.
        3 - The simulation then returns the solution vector at the end of the 
            integration, which can be processed as needed to plot or analyze any quantities of interest.
"""

#
from ZAB_init import *
from scipy.integrate import solve_ivp
import numpy as np
from matplotlib import pyplot as plt
from math import exp, log

##################################################################################################
"Calculation of Exchange current density of the Cathode"

def prop(C_k_ca):
    pi_rev_ca = C_k_ca[3]**2 
    pi_fwd_ca = C_k_ca[0]**0.5*C_k_ca[1]

    # Gibbs free energy calculation
    delg_ca = delg_0_ca + R*T*log(pi_rev_ca/pi_fwd_ca)
    dphi_eq_ca = delg_ca /(-n_ca*F)

    k_fwd_ca = k_f_eq_ca*exp(-alpha_ca*n_ca*F*dphi_dl_ca_0/(R*T))
    k_rev_ca = k_r_eq_ca*exp((1-alpha_ca)*n_ca*F*dphi_dl_ca_0/(R*T))

    i_o_ca = n_ca*F*k_fwd_ca**(1-alpha_ca)*k_rev_ca**(alpha_ca)*pi_fwd_ca**(1-alpha_ca)*pi_rev_ca**alpha_ca
    
    
    ss = np.array((dphi_eq_ca, i_o_ca))
    return ss
    
##################################################################################################
"Anode properties"
i_o_an  = 2.5 # Fixed by user
dphi_eq_an = .8 # Fixed bu user

SV_0 = np.array([dphi_dl_an_0, dphi_dl_ca_0])
time_span = np.array([0, t_final])

##################################################################################################
"Calculation of electric potentials of the cathode and the electrolyte - CO2 included"

def residual_0(t,SV):
    dSV_dt = np.zeros_like(SV)
        # Anode Interface:
    
    eta_an = SV[0] - dphi_eq_an
    i_Far_an = i_o_an*(exp(-n_an*F*alpha_an*eta_an/R/T)
                            - exp(n_an*F*(1-alpha_an)*eta_an/R/T))
    i_dl_an = i_ext - i_Far_an
    dSV_dt[0] = -i_dl_an/C_dl_an

        # Cathode Interface:
    eta_ca = SV[1] - prop(C_wCO2)[0]
    i_Far_ca = prop(C_wCO2)[1]*(exp(-n_ca*F*alpha_ca*eta_ca/R/T)
                          - exp(n_ca*F*(1-alpha_ca)*eta_ca/R/T))
    i_dl_ca = i_ext - i_Far_ca

    dSV_dt[1] = -i_dl_ca/C_dl_ca
    return dSV_dt


sol_wCO2 = solve_ivp(residual_0,time_span,SV_0,rtol=1e-4, atol=1e-6)

##################################################################################################
"Calculation of electric potentials of the cathode and the electrolyte - No CO2"

def residual_1(t,SV):
    dSV_dt = np.zeros_like(SV)
        # Anode Interface:
    
    eta_an = SV[0] - dphi_eq_an
    i_Far_an = i_o_an*(exp(-n_an*F*alpha_an*eta_an/R/T)
                            - exp(n_an*F*(1-alpha_an)*eta_an/R/T))
    i_dl_an = i_ext - i_Far_an
    dSV_dt[0] = -i_dl_an/C_dl_an

        # Cathode Interface:
    eta_ca = SV[1] - prop(C_noCO2)[0]
    i_Far_ca = prop(C_noCO2)[1]*(exp(-n_ca*F*alpha_ca*eta_ca/R/T)
                          - exp(n_ca*F*(1-alpha_ca)*eta_ca/R/T))
    i_dl_ca = i_ext - i_Far_ca


    dSV_dt[1] = -i_dl_ca/C_dl_ca
    return dSV_dt


sol_noCO2 = solve_ivp(residual_1,time_span,SV_0,rtol=1e-4, atol=1e-6)



V_elyte = sol_wCO2.y[0,:]
V_ca = V_elyte + sol_wCO2.y[1,:]
plt.plot(sol_wCO2.t,V_elyte)
plt.plot(sol_wCO2.t,V_ca)
plt.xlabel('Time (s)',fontsize=14)
plt.title ('CO2 included')
plt.ylabel('Electric Potential (V)',fontsize=14)
plt.legend(['Electrolyte','Cathode'])
plt.show()

V_elyte1 = sol_noCO2.y[0,:]
V_ca1 = V_elyte1 + sol_noCO2.y[1,:]
plt.plot(sol_noCO2.t,V_elyte1)
plt.plot(sol_noCO2.t,V_ca1)
plt.title ('No CO2')
plt.xlabel('Time (s)',fontsize=14)
plt.ylabel('Electric Potential (V)',fontsize=14)
plt.legend(['Electrolyte','Cathode'])
plt.show()


###########################################################################################################

"Cell potential calculation for varying current"
def ZA_cellvoltage(resid, i_ext, SV_0, plot_flag): 
    time_span = np.array([0,100])
    soln = solve_ivp(resid,time_span,SV_0,rtol=1e-4, atol=1e-6)

    return soln.y[:,-1]

i_array = np.linspace(0, 100, 100)
V_cell_0 = np.zeros_like(i_array)
V_cell_1 = np.zeros_like(i_array)

SV = np.zeros((SV_0.size,i_array.size))

for j,i_ext in enumerate(i_array):
    plot = 0
    SV[:,j] = ZA_cellvoltage(residual_0,i_ext, SV_0, plot)
    SV_0 = SV[:,j]


V_cell_0 = SV[0,:] + SV[1,:]
plt.plot(i_array,V_cell_0)
plt.title ('CO2 included')
plt.xlabel('Current',fontsize=14)
plt.ylabel('Cell Potential (V)',fontsize=14)
plt.show()



for j,i_ext in enumerate(i_array):
    plot = 0
    SV[:,j] = ZA_cellvoltage(residual_1,i_ext, SV_0, plot)
    SV_0 = SV[:,j]


V_cell_1 = SV[0,:] + SV[1,:]
plt.plot(i_array,V_cell_1)
plt.title ('No CO2')
plt.xlabel('Current (A)',fontsize=14)
plt.ylabel('Cell Potential (V)',fontsize=14)
plt.show()

########################################################################################
"  Fluc calculation of different species for ZAB cathode"

def transport(X_k_ca):
    
    Kg = eps_air_ca**3*d_part_ca**2/(72*eps_air_ca**-1*(1-eps_air_ca)**2)
    Vconv = -Kg*(P_sep-P_ca)/mu_air/H_ca

    xavg = np.zeros_like(X_k_ca)
    xavg = (X_k_ca + X_k_sep)/2
    
    Dkeff = eps_air_ca**1.5*D_k
    dk = (X_k_sep - X_k_ca)/H_ca
    Vdif =  -1*Dkeff*dk/xavg
    Vtot = Vconv + Vdif
    C_k = xavg*P_ca/R/T
    N_k = C_k*Vtot
    return N_k


N_wCO2 = np.array([transport(X_wCO2)[0], transport(X_wCO2)[1], transport(X_wCO2)[2]])
N_noCO2 = np.array([transport(X_noCO2)[0], transport(X_noCO2)[1], transport(X_noCO2)[2]])

fig, ax = plt.subplots()

labels = ['O2', 'H2O', 'N2']
x = np.arange(len(labels))

width = 0.35
ax.bar(x+width/2,N_wCO2,width)
ax.bar(x-width/2,N_noCO2,width)
ax.legend(['CO2 included','No CO2'],frameon=False)

ax.set_xticks(x)
ax.set_xticklabels(labels)

ax.set_ylabel('Flux (mol/m$^2$/s)',fontsize=14)
ax.set_xlabel('Species',fontsize=14)

