import numpy as n
import matplotlib.pyplot as pyp

## Import the modules
import DarkMatterUtilities.Constants as dmcon
import DarkMatterUtilities.Targets as dmt

## Initialize a complex target
nat_xe_target = dmt.NaturalXenonTarget()

## Create a plot of the Lindhard
## factor (fraction of recoil
## energy into signal) vs recoil
## energy.
## Result is calculated as a 
## weighted sum over isotopic
## abundances.
## ==============================
E_recoil_keV = n.linspace(start=0.0 , stop=30.0, num=1000)
L_factors = n.zeros(len(E_recoil_keV))
for i in n.arange(len(E_recoil_keV)):
	L_factors[i] = nat_xe_target.LindhardFactor(E_recoil_keV[i])

pyp.figure()
pyp.plot(E_recoil_keV, L_factors)

pyp.xlabel("Recoil energy [keV]")
pyp.ylabel("Lindhard factor")
pyp.title(nat_xe_target.Name)

## Create a plot of the angular
## dependence of recoil energy
## as a fraction of incoming
## particle energy (default 
## option is incoming neutron, 
## can change this with the 
## argument: 
## _incoming_particle_mass_GeV ).
## Result is calculated as a 
## weighted sum over isotopic
## abundances.
## ==============================
angle_vals = n.linspace(start=-n.pi , stop=n.pi, num=1000)
frac_vals = n.zeros(len(angle_vals))
for i in n.arange(len(angle_vals)):
	frac_vals[i] = nat_xe_target.RecoilEnergyAngularDist(angle_vals[i])

pyp.figure()
pyp.plot(angle_vals, frac_vals)

pyp.xlabel("Angle from Incident Particle [rad]")
pyp.ylabel("Recoil Energy Fraction")
pyp.title(nat_xe_target.Name)

## Create a plot of the maximum
## recoil energy of a nucleus
## for an incoming  particle with 
## a specific mass as a function
## of incomint particle energy. 
## Result is calculated as a 
## weighted sum over isotopic
## abundances.
## ==============================
in_energy_vals  = n.linspace(start=0.0, stop=100.0, num=1000)
out_energy_vals = n.zeros(len(in_energy_vals))
for i in n.arange(len(in_energy_vals)):
	out_energy_vals[i] = nat_xe_target.RecoilEnergyMax_AnyParticle_keV(dmcon.m_neutron_GeV , in_energy_vals[i])

pyp.figure()
pyp.plot(in_energy_vals, out_energy_vals, label="Neutron, M="+str(dmcon.m_neutron_GeV)+" GeV")
pyp.legend(loc='lower right')

pyp.xlabel("Incident Particle Energy [keV]")
pyp.ylabel("Maximum Recoil Energy [keV]")
pyp.title(nat_xe_target.Name)

## Create a plot of the maximum
## recoil energy of a nucleus
## for dark matter (uses the MW
## escape velocity for maximal
## incident energy).
## Result is calculated as a 
## weighted sum over isotopic
## abundances.
## ==============================
DM_massGeV_vals = n.linspace(start=5.0, stop=1000.0, num=1000)
out_energy_vals = n.zeros(len(DM_massGeV_vals))
for i in n.arange(len(DM_massGeV_vals)):
	out_energy_vals[i] = nat_xe_target.RecoilEnergyMax_DM_keV(DM_massGeV_vals[i])

pyp.figure()
pyp.plot(DM_massGeV_vals, out_energy_vals)

pyp.xlabel("Incident Dark Matter Mass [GeV"+r"/c$^2$"+"]")
pyp.ylabel("Maximum Recoil Energy [keV]")
pyp.title(nat_xe_target.Name)

## Create a plot of the minimum
## mass a detector of this target
## is sensitive to for a given
## energy threshold.
## Result is calculated as a 
## weighted sum over isotopic
## abundances.
## ==============================
E_thresh_keV_vals = n.logspace(start=-2, stop=2, num=10000)
min_DM_mass_GeV  = n.zeros(len(E_thresh_keV_vals))
for i in n.arange(len(E_thresh_keV_vals)):
	min_DM_mass_GeV[i] = nat_xe_target.MinimumDetectableMass_GeV(E_thresh_keV_vals[i])

pyp.figure()
pyp.plot(E_thresh_keV_vals, min_DM_mass_GeV)

# pyp.xscale('log')
# pyp.yscale('log')

pyp.xlabel("Detector threshold [keV]")
pyp.ylabel("Minimum Detectable Dark Matter Mass [GeV"+r"/c$^2$"+"]")
pyp.title(nat_xe_target.Name)

## ==============================
pyp.show()
