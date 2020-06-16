import numpy as n
import matplotlib.pyplot as pyp

## Import the modules
import DarkMatterUtilities.Constants as dmcon
import DarkMatterUtilities.Targets as dmt

## Initialize a few simple targets
avg_xe_target = dmt.SimpleTarget( 131.293 , 54.0 , "Xe")
avg_ge_target = dmt.SimpleTarget(  72.630 , 32.0 , "Ge")
avg_ar_target = dmt.SimpleTarget(  39.948 , 18.0 , "Ar")
avg_si_target = dmt.SimpleTarget(  28.085 , 14.0 , "Si")
avg_ne_target = dmt.SimpleTarget(  20.180 , 10.0 , "Ne")
avg_he_target = dmt.SimpleTarget(   4.0   ,  2.0 , "He")
avg_h_target  = dmt.SimpleTarget(   1.0   ,  1.0 , "H" )
avg_e_target  = dmt.SimpleTarget( 0.00054858,1.0 , "e-")

## Create a plot of the Lindhard
## factor (fraction of recoil
## energy into signal) vs recoil
## energy.
## ==============================
E_recoil_keV = n.linspace(start=0.0 , stop=30.0, num=1000)
L_factors_xe = n.zeros(len(E_recoil_keV))
L_factors_ge = n.zeros(len(E_recoil_keV))
L_factors_ar = n.zeros(len(E_recoil_keV))
L_factors_si = n.zeros(len(E_recoil_keV))
L_factors_ne = n.zeros(len(E_recoil_keV))
L_factors_he = n.zeros(len(E_recoil_keV))
L_factors_h  = n.zeros(len(E_recoil_keV))
L_factors_e  = n.zeros(len(E_recoil_keV))

for i in n.arange(len(E_recoil_keV)):
	L_factors_xe[i] = avg_xe_target.LindhardFactor(E_recoil_keV[i])
	L_factors_ge[i] = avg_ge_target.LindhardFactor(E_recoil_keV[i])
	L_factors_ar[i] = avg_ar_target.LindhardFactor(E_recoil_keV[i])
	L_factors_si[i] = avg_si_target.LindhardFactor(E_recoil_keV[i])
	L_factors_ne[i] = avg_ne_target.LindhardFactor(E_recoil_keV[i])
	L_factors_he[i] = avg_he_target.LindhardFactor(E_recoil_keV[i])
	L_factors_h[i]  = avg_h_target.LindhardFactor(E_recoil_keV[i])
	L_factors_e[i]  = avg_e_target.LindhardFactor(E_recoil_keV[i])

pyp.figure()

pyp.plot(E_recoil_keV, L_factors_xe, label=avg_xe_target.Name)
pyp.plot(E_recoil_keV, L_factors_ge, label=avg_ge_target.Name)
pyp.plot(E_recoil_keV, L_factors_ar, label=avg_ar_target.Name)
pyp.plot(E_recoil_keV, L_factors_si, label=avg_si_target.Name)
pyp.plot(E_recoil_keV, L_factors_ne, label=avg_ne_target.Name)
pyp.plot(E_recoil_keV, L_factors_he, label=avg_he_target.Name)
pyp.plot(E_recoil_keV, L_factors_h, label=avg_h_target.Name)
pyp.plot(E_recoil_keV, L_factors_e, label=avg_e_target.Name)

pyp.legend(loc="lower right")
pyp.xlabel("Recoil energy [keV]")
pyp.ylabel("Lindhard factor")

## Create a plot of the angular
## dependence of recoil energy
## as a fraction of incoming
## particle energy (default 
## option is incoming neutron, 
## can change this with the 
## argument: 
## _incoming_particle_mass_GeV ).
## ==============================
angle_vals = n.linspace(start=-n.pi , stop=n.pi, num=1000)

frac_vals_xe = n.zeros(len(angle_vals))
frac_vals_ge = n.zeros(len(angle_vals))
frac_vals_ar = n.zeros(len(angle_vals))
frac_vals_si = n.zeros(len(angle_vals))
frac_vals_ne = n.zeros(len(angle_vals))
frac_vals_he = n.zeros(len(angle_vals))
frac_vals_h  = n.zeros(len(angle_vals))
frac_vals_e  = n.zeros(len(angle_vals))

for i in n.arange(len(angle_vals)):
	frac_vals_xe[i] = avg_xe_target.RecoilEnergyAngularDist(angle_vals[i])
	frac_vals_ge[i] = avg_ge_target.RecoilEnergyAngularDist(angle_vals[i])
	frac_vals_ar[i] = avg_ar_target.RecoilEnergyAngularDist(angle_vals[i])
	frac_vals_si[i] = avg_si_target.RecoilEnergyAngularDist(angle_vals[i])
	frac_vals_ne[i] = avg_ne_target.RecoilEnergyAngularDist(angle_vals[i])
	frac_vals_he[i] = avg_he_target.RecoilEnergyAngularDist(angle_vals[i])
	frac_vals_h[i]  = avg_h_target.RecoilEnergyAngularDist(angle_vals[i])
	frac_vals_e[i]  = avg_e_target.RecoilEnergyAngularDist(angle_vals[i])

pyp.figure()
pyp.plot(angle_vals, frac_vals_xe, label=avg_xe_target.Name)
pyp.plot(angle_vals, frac_vals_ge, label=avg_ge_target.Name)
pyp.plot(angle_vals, frac_vals_ar, label=avg_ar_target.Name)
pyp.plot(angle_vals, frac_vals_si, label=avg_si_target.Name)
pyp.plot(angle_vals, frac_vals_ne, label=avg_ne_target.Name)
pyp.plot(angle_vals, frac_vals_he, label=avg_he_target.Name)
pyp.plot(angle_vals, frac_vals_h, label=avg_h_target.Name)
pyp.plot(angle_vals, frac_vals_e, label=avg_e_target.Name)

pyp.legend(loc="lower right")
pyp.xlabel("Angle from Incident Particle [rad]")
pyp.ylabel("Recoil Energy Fraction")

## Create a plot of the maximum
## recoil energy of a nucleus
## for an incoming  particle with 
## a specific mass as a function
## of incomint particle energy. 
## ==============================
in_energy_vals  = n.linspace(start=0.0, stop=100.0, num=1000)

out_energy_vals_xe = n.zeros(len(in_energy_vals))
out_energy_vals_ge = n.zeros(len(in_energy_vals))
out_energy_vals_ar = n.zeros(len(in_energy_vals))
out_energy_vals_si = n.zeros(len(in_energy_vals))
out_energy_vals_ne = n.zeros(len(in_energy_vals))
out_energy_vals_he = n.zeros(len(in_energy_vals))
out_energy_vals_h  = n.zeros(len(in_energy_vals))
out_energy_vals_e  = n.zeros(len(in_energy_vals))

for i in n.arange(len(in_energy_vals)):
	out_energy_vals_xe[i] = avg_xe_target.RecoilEnergyMax_AnyParticle_keV(dmcon.m_neutron_GeV , in_energy_vals[i])
	out_energy_vals_ge[i] = avg_ge_target.RecoilEnergyMax_AnyParticle_keV(dmcon.m_neutron_GeV , in_energy_vals[i])
	out_energy_vals_ar[i] = avg_ar_target.RecoilEnergyMax_AnyParticle_keV(dmcon.m_neutron_GeV , in_energy_vals[i])
	out_energy_vals_si[i] = avg_si_target.RecoilEnergyMax_AnyParticle_keV(dmcon.m_neutron_GeV , in_energy_vals[i])
	out_energy_vals_ne[i] = avg_ne_target.RecoilEnergyMax_AnyParticle_keV(dmcon.m_neutron_GeV , in_energy_vals[i])
	out_energy_vals_he[i] = avg_he_target.RecoilEnergyMax_AnyParticle_keV(dmcon.m_neutron_GeV , in_energy_vals[i])
	out_energy_vals_h[i]  = avg_h_target.RecoilEnergyMax_AnyParticle_keV(dmcon.m_neutron_GeV , in_energy_vals[i])
	out_energy_vals_e[i]  = avg_e_target.RecoilEnergyMax_AnyParticle_keV(dmcon.m_neutron_GeV , in_energy_vals[i])

pyp.figure()
pyp.plot(in_energy_vals, out_energy_vals_xe, label=avg_xe_target.Name)
pyp.plot(in_energy_vals, out_energy_vals_ge, label=avg_ge_target.Name)
pyp.plot(in_energy_vals, out_energy_vals_ar, label=avg_ar_target.Name)
pyp.plot(in_energy_vals, out_energy_vals_si, label=avg_si_target.Name)
pyp.plot(in_energy_vals, out_energy_vals_ne, label=avg_ne_target.Name)
pyp.plot(in_energy_vals, out_energy_vals_he, label=avg_he_target.Name)
pyp.plot(in_energy_vals, out_energy_vals_h, label=avg_h_target.Name)
pyp.plot(in_energy_vals, out_energy_vals_e, label=avg_e_target.Name)

pyp.legend(loc='lower right')
pyp.xlabel("Incident Particle Energy [keV]")
pyp.ylabel("Maximum Recoil Energy [keV]")
pyp.title("Neutron, M="+str(dmcon.m_neutron_GeV)+" GeV")

## Create a plot of the maximum
## recoil energy of a nucleus
## for dark matter (uses the MW
## escape velocity for maximal
## incident energy).
## ==============================
DM_massGeV_vals = n.linspace(start=5.0, stop=1000.0, num=1000)

out_energy_vals_xe = n.zeros(len(DM_massGeV_vals))
out_energy_vals_ge = n.zeros(len(DM_massGeV_vals))
out_energy_vals_ar = n.zeros(len(DM_massGeV_vals))
out_energy_vals_si = n.zeros(len(DM_massGeV_vals))
out_energy_vals_ne = n.zeros(len(DM_massGeV_vals))
out_energy_vals_he = n.zeros(len(DM_massGeV_vals))
out_energy_vals_h  = n.zeros(len(DM_massGeV_vals))
out_energy_vals_e  = n.zeros(len(DM_massGeV_vals))

for i in n.arange(len(DM_massGeV_vals)):
	out_energy_vals_xe[i] = avg_xe_target.RecoilEnergyMax_DM_keV(DM_massGeV_vals[i])
	out_energy_vals_ge[i] = avg_ge_target.RecoilEnergyMax_DM_keV(DM_massGeV_vals[i])
	out_energy_vals_ar[i] = avg_ar_target.RecoilEnergyMax_DM_keV(DM_massGeV_vals[i])
	out_energy_vals_si[i] = avg_si_target.RecoilEnergyMax_DM_keV(DM_massGeV_vals[i])
	out_energy_vals_ne[i] = avg_ne_target.RecoilEnergyMax_DM_keV(DM_massGeV_vals[i])
	out_energy_vals_he[i] = avg_he_target.RecoilEnergyMax_DM_keV(DM_massGeV_vals[i])
	out_energy_vals_h[i]  = avg_h_target.RecoilEnergyMax_DM_keV(DM_massGeV_vals[i])
	out_energy_vals_e[i]  = avg_e_target.RecoilEnergyMax_DM_keV(DM_massGeV_vals[i])

pyp.figure()
pyp.plot(DM_massGeV_vals, out_energy_vals_xe, label=avg_xe_target.Name)
pyp.plot(DM_massGeV_vals, out_energy_vals_ge, label=avg_ge_target.Name)
pyp.plot(DM_massGeV_vals, out_energy_vals_ar, label=avg_ar_target.Name)
pyp.plot(DM_massGeV_vals, out_energy_vals_si, label=avg_si_target.Name)
pyp.plot(DM_massGeV_vals, out_energy_vals_ne, label=avg_ne_target.Name)
pyp.plot(DM_massGeV_vals, out_energy_vals_he, label=avg_he_target.Name)
pyp.plot(DM_massGeV_vals, out_energy_vals_h, label=avg_h_target.Name)
pyp.plot(DM_massGeV_vals, out_energy_vals_e, label=avg_e_target.Name)

pyp.legend(loc='lower right')
pyp.xlabel("Incident Dark Matter Mass [GeV"+r"/c$^2$"+"]")
pyp.ylabel("Maximum Recoil Energy [keV]")

## Create a plot of the minimum
## mass a detector of this target
## is sensitive to for a given
## energy threshold.
## ==============================
E_thresh_keV_vals = n.logspace(start=-2, stop=2, num=10000)

min_DM_mass_GeV_xe = n.zeros(len(E_thresh_keV_vals))
min_DM_mass_GeV_ge = n.zeros(len(E_thresh_keV_vals))
min_DM_mass_GeV_ar = n.zeros(len(E_thresh_keV_vals))
min_DM_mass_GeV_si = n.zeros(len(E_thresh_keV_vals))
min_DM_mass_GeV_ne = n.zeros(len(E_thresh_keV_vals))
min_DM_mass_GeV_he = n.zeros(len(E_thresh_keV_vals))
min_DM_mass_GeV_h  = n.zeros(len(E_thresh_keV_vals))
min_DM_mass_GeV_e  = n.zeros(len(E_thresh_keV_vals))

for i in n.arange(len(E_thresh_keV_vals)):
	min_DM_mass_GeV_xe[i] = avg_xe_target.MinimumDetectableMass_GeV(E_thresh_keV_vals[i])
	min_DM_mass_GeV_ge[i] = avg_ge_target.MinimumDetectableMass_GeV(E_thresh_keV_vals[i])
	min_DM_mass_GeV_ar[i] = avg_ar_target.MinimumDetectableMass_GeV(E_thresh_keV_vals[i])
	min_DM_mass_GeV_si[i] = avg_si_target.MinimumDetectableMass_GeV(E_thresh_keV_vals[i])
	min_DM_mass_GeV_ne[i] = avg_ne_target.MinimumDetectableMass_GeV(E_thresh_keV_vals[i])
	min_DM_mass_GeV_he[i] = avg_he_target.MinimumDetectableMass_GeV(E_thresh_keV_vals[i])
	min_DM_mass_GeV_h[i]  = avg_h_target.MinimumDetectableMass_GeV(E_thresh_keV_vals[i])
	min_DM_mass_GeV_e[i]  = avg_e_target.MinimumDetectableMass_GeV(E_thresh_keV_vals[i])

pyp.figure()
pyp.plot(E_thresh_keV_vals, min_DM_mass_GeV_xe, label=avg_xe_target.Name)
pyp.plot(E_thresh_keV_vals, min_DM_mass_GeV_ge, label=avg_ge_target.Name)
pyp.plot(E_thresh_keV_vals, min_DM_mass_GeV_ar, label=avg_ar_target.Name)
pyp.plot(E_thresh_keV_vals, min_DM_mass_GeV_si, label=avg_si_target.Name)
pyp.plot(E_thresh_keV_vals, min_DM_mass_GeV_ne, label=avg_ne_target.Name)
pyp.plot(E_thresh_keV_vals, min_DM_mass_GeV_he, label=avg_he_target.Name)
pyp.plot(E_thresh_keV_vals, min_DM_mass_GeV_h, label=avg_h_target.Name)
pyp.plot(E_thresh_keV_vals, min_DM_mass_GeV_e, label=avg_e_target.Name)

# pyp.xscale('log')
# pyp.yscale('log')

pyp.legend(loc='lower right')
pyp.xlabel("Detector threshold [keV]")
pyp.ylabel("Minimum Detectable Dark Matter Mass [GeV"+r"/c$^2$"+"]")

# ## ==============================
pyp.show()
