# Python code to make the recoil spectrum in Xe
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

# dR/dE_{R} = ( rho_0 * sigma_A / 2 * M_chi * mu_A^2 ) * F^2(q) * [integral^(v_max)_(v_min) {f(v)/v} d^3 v]
# rho_0    -> kg m^-3
# sigma_A  -> m^2
# M_chi    -> keV m^-2 s^2
# mu_A     -> kg^2
# F^2(q)   -> dimensionless
# integral -> m^-1 s

# ====> kg m^-1 keV^-1 m^2 s^-2 kg^-2 m^-1 s
# ====> kg^-1 keV^-1 s^-1


# energies
E_LOW   = 0.1  # keV
E_HIGH  = 200  # keV
keV2GeV = 1e-6 # 1 keV = 1e-6 GeV
BIN     = 0.1  # keV
NBINS   = ((E_HIGH-E_LOW)/BIN)+1
E_R     =np.linspace(E_LOW*keV2GeV,E_HIGH*keV2GeV,NBINS)    # GeV (1e-6 GeV = 1 keV, 2e-4 = 200 keV)


# error function defintion
def Nesc(errZ):
    return (math.erf(errZ)-((2*errZ*math.exp(-pow(errZ,2)))/pow(math.pi,0.5)))

#########################################

A       =131.293             # mass number of xenon
amu2gev =0.931149            # amu -> GeV/c^2
m_N     =A*amu2gev           # nucleus mass (GeV/c^2)
m_n     =amu2gev             # nucleon mass
m_X     =50.                 # WIMP mass (GeV/c^2)
mu_N    =(m_N*m_X)/(m_N+m_X) # reduced mass of the WIMP-nucleus system
mu_n    =(m_n*m_X)/(m_n+m_X) # reduced mass of the WIMP-nucleon system
#sigma_n =1.e-36             # WIMP-nucleon cross-section (cm^2)
sigma_n = 1.0e-45  # (roughly) Xe1T sensitivity at 50 GeV/c^2
rho_0   =0.3                 # local dark matter density GeV/c^2 cm^-3                

# gives us the zero momentum-transfer or 'standard' cross section
sigma_A = A*A*(mu_N/mu_n)*(mu_N/mu_n)*sigma_n # cm^2
print(sigma_A)
# compute nuclear form factor (energy dependent)
s             =1.e-15                         # m
R             =1.2*pow(A,(1./3.))*1e-15       # m
r             =pow((pow(R,2)-5*pow(s,2)),0.5) # m
h_bar_MeV_fm  =197.3                          # MeV*fm (see Lewin-Smith, p8)
h_bar_GeV_m   =197.3*1e-15*1e-3               # GeV*m
F2q = []

# astrophysical constants
v_esc   =544. # km /s
v_0     =220. # km /s
v_earth =245. # km /s
y       =v_earth/v_0
z       =v_esc/v_0
vMin = []
integral = []
contrib  =0
GeV2J        =1.6022e-10  # GeV -> J
csqrd_const  =3.e8*3.e8   # m^2 /s^2

index    =0        
Rate     = []
for e_r in E_R:
    # form factor
    q             =pow((2*m_N*e_r),0.5)           # GeV (natural units, c = 1)
    dimlessqr     =(q*r)/h_bar_GeV_m              # q*r / h_bar_GeV_m is dimensionless
    dimlessqs     =(q*s)/h_bar_GeV_m              # q*s / h_bar_GeV_m is dimensionless
    j             =(math.sin(dimlessqr)/(pow(dimlessqr,2)))-(math.cos(dimlessqr)/dimlessqr)  # bessel function
    f2q           =pow(((3*j)/(dimlessqr)),2)*math.exp(-pow(dimlessqs,2))
    F2q.append(f2q)

    # velocity integration
    m_N_kg  = (m_N*GeV2J)/csqrd_const
    e_r_J   = e_r*GeV2J
    mu_N_kg = (mu_N*GeV2J)/csqrd_const
    v_min   = (pow(((m_N_kg*e_r_J)/(2*pow((mu_N_kg),2))),0.5))/1.e3
    x       = v_min/v_0    
    vMin.append(v_min)
    if x<abs(y-z):
        contrib=(1./(2*Nesc(z)*v_0*y))*((math.erf(x+y))-(math.erf(x-y))-((4./(pow(math.pi,0.5)))*y*math.exp(-pow(z,2))))
    if x>abs(y-z) and x<(y+z):
        contrib=(1./(2*Nesc(z)*v_0*y))*((math.erf(z))-(math.erf(x-y))-((2./(pow(math.pi,0.5)))*(y+z-x)*math.exp(-pow(z,2))))    
    if x>(y+z):
        contrib=0
    integral.append(contrib)
        
    # combine
    # rate = ((f2q*rho_0*sigma_A)/(2*m_X*mu_N*mu_N))*contrib
    rate = (rho_0/m_X) * (sigma_A / (2*mu_N*mu_N))  * f2q * contrib
    Rate.append(rate)

# dR/dE_{R} = ( rho_0 * sigma_A / 2 * M_chi * mu_A^2 ) * F^2(q) * [integral^(v_max)_(v_min) {f(v)/v} d^3 v]
# rho_0    -> kg m^-3
# sigma_A  -> m^2
# M_chi    -> keV m^-2 s^2
# mu_A^2   -> kg^2
# F^2(q)   -> dimensionless
# integral -> m^-1 s 

rho_0_cnv    = GeV2J*1.e6/csqrd_const          # GeV c^-2 cm^-3 /->/ kg m^2 s^-2 m^-2 s^2 cm^-3 /->/ kg cm^-3 /->/ kg m^-3
sigma_A_cnv  = 1.e-4                           # cm^2 -> m^2 
M_chi_cnv    = 1.e6/csqrd_const                # GeV c^-2 /->/ keV m^-2 s^2
mu_A_cnv     = pow(GeV2J,2)/pow(csqrd_const,2) # GeV^2 c^-4 /->/ kg^2 m^4 s^-4 m^-4 s^-4 /->/ kg^2
integral_cnv = 1.e-3                           # km^-1 s /->/ m^-1 s
sec2day_cnv  = 86400.                          # s^-1 to day^-1
tot_cnv      = rho_0_cnv*sigma_A_cnv*integral_cnv*sec2day_cnv/(M_chi_cnv*mu_A_cnv)

print "Total conversion = ",tot_cnv

Rate_kgdkeV =[i*tot_cnv for i in Rate]
Rate_kgd    =sum(Rate_kgdkeV)*(BIN)
LZEvents    =Rate_kgd*7000.*1000.

print "Total number of events per kg per day",Rate_kgd
print "Total number of LZ events = ",LZEvents
print "Total number of LZ events per day = ",LZEvents/1000.

# plot
# plt.plot(E_R*1e6,F2q)
# plt.yscale("log")
# plt.xlabel('Energy (keV)')
# plt.ylabel('Form Factor')
# plt.savefig("FormFactor.png")
# plt.show()

# plt.plot(vMin,integral)
# plt.xlabel('v_{min} [km s^{-1}]')
# plt.ylabel('Integral of f(v)/v above v_{min} [km s^{-1}]')
# plt.savefig("Integral.png")
# plt.show()

plt.plot(E_R*1e6,Rate_kgdkeV)
plt.yscale("log")
plt.xlabel('Energy [keV]')
plt.ylabel('Rate [cts/kg/day/keV]')
plt.savefig("RecoilSpec.png")
plt.show()
