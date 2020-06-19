import numpy as n
import matplotlib.pyplot as pyp

## Import the modules
import DarkMatterUtilities.Constants as dmcon
import DarkMatterUtilities.DarkMatter as dmdm
import DarkMatterUtilities.Targets as dmt
import DarkMatterUtilities.FormFactors as dmff
import DarkMatterUtilities.Interactions as dmi

## Initialize a dark matter model with explicit Halo Model
WIMP_50GeV    = dmdm.DarkMatter(50.0 , 1.0e-45 , _HaloModelType=1)

## Initialize a simple xenon target
avg_xe_target = dmt.SimpleTarget( 131.293 , 54.0 , "Xe", 4.7808 )
# avg_xe_target = dmt.NaturalXenonTarget()

## Initialize a Helm form factor for the xenon target
ff_helm       = dmff.FormFactor(avg_xe_target, 1)
avg_xe_target.FormFactor = ff_helm

## Create an interaction class
interaction   = dmi.Interaction(avg_xe_target, WIMP_50GeV, "SI")

## Create a plot of the minimum
## velocity required to make a 
## recoil of energy E_r keV
## ==============================
E_range_keV   = n.linspace(start=0.0, stop=100.0, num=1000)

min_vel_1_ms  = n.zeros(len(E_range_keV))
min_vel_2_ms  = n.zeros(len(E_range_keV))

for i in n.arange(len(E_range_keV)):
	min_vel_1_ms[i] = interaction.MinimumVelocity_NoInelasticity_ms(E_range_keV[i])
	min_vel_2_ms[i] = interaction.MinimumVelocity_ms(E_range_keV[i])

pyp.figure()
pyp.plot(E_range_keV, min_vel_1_ms, label="No inelasticity")
pyp.plot(E_range_keV, min_vel_2_ms, label="Inelastic ("+r"$\delta=$"+str(WIMP_50GeV.delta_keV)+" keV")

pyp.legend(loc='lower right')
pyp.xlabel(r"Recoil Energy $E_r$ [keV]")
pyp.ylabel(r"Minimum Velocity to Create Recoil of $E_r$ [m/s]")

## Create a plot of the minimum
## velocity required to make a 
## recoil of energy E_r keV
## ==============================
E_range_keV   = n.linspace(start=0.0, stop=100.0, num=1000)

SI_xsec       = n.zeros(len(E_range_keV))

for i in n.arange(len(E_range_keV)):
	SI_xsec[i] = interaction.GetDifferentialCrossSection(E_range_keV[i])

pyp.figure()
pyp.plot(E_range_keV, SI_xsec, label="Working")

pyp.legend(loc='lower right')
pyp.xlabel("Recoil Energy [keV]")
pyp.ylabel(r"$d\sigma / dE_{r}$ [cm$^2/$keV]")


## Create a plot of the differential
## event rate per day per unit 
## mass (kg)
## ==============================
E_range_keV   = n.linspace(start=0.0, stop=1000.0, num=1000)

dRdE          = n.zeros(len(E_range_keV))

for i in n.arange(len(E_range_keV)):
	dRdE[i] = interaction.DifferentialRate(E_range_keV[i])

pyp.figure()
pyp.plot(E_range_keV, dRdE, label=avg_xe_target.Name)

pyp.xscale('log')
pyp.yscale('log')

pyp.ylim([1e-12 , 1e-3])

pyp.legend(loc='lower right')
pyp.xlabel("Recoil Energy [keV]")
pyp.ylabel(r"Differential Rate [keV$^{-1}\times$ kg$^{-1}\times$ day$^{-1}$]")

## ==============================
pyp.show()
