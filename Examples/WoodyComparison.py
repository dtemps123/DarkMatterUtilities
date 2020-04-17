import numpy as n
import matplotlib.pyplot as pyp
import math

from DMClasses      import *
from DMInteractions import *
from DMPlots        import *

pyp.rcParams['axes.grid'] = True
pyp.rcParams.update({'font.size': 16})
pyp.rc('text', usetex=True)
pyp.rc('font', family='serif')

DM_mass = 50.0		# GeV/c^2
DM_xsec = 1.0e-45	# cm^2

Xe_target_DT = Target( 131.293 , 54.0 , 7.0e3 , 1.0 , "Xe", 4.7808 )
Xe_target_DT.FF_type = 4	## Use my version of Helm Form Factor

Xe_target_DW = Target( 131.293 , 54.0 , 7.0e3 , 1.0 , "Xe", 4.7808 )
Xe_target_DW.FF_type = 5	## Use David's version of Helm Form Factor

DM_DT = DarkMatter(DM_mass, DM_xsec)
DM_DT.HaloModel = Halo(1)	## Use my version of SHM integral

DM_DW = DarkMatter(DM_mass, DM_xsec)
DM_DW.HaloModel = Halo(2)	## Use David's version of SHM integral

# get a recoil energy array [keV]
plot_Erange = n.logspace(start=-1.0, stop=3.0  , num=1000)
# plot_Erange = n.linspace(start= 0.1, stop=2.0e2, num=1000)
plot_qrn    = n.sqrt(2.0 * Xe_target_DT.NuclearMass_GeV * plot_Erange) * ( Xe_target_DT.FF_Rn / hbarc_MeV_fm) 
plot_vMins  = n.zeros(len(plot_Erange))

plot_ff_Helm   = n.zeros(len(plot_Erange))
plot_ff_HelmDW = n.zeros(len(plot_Erange))

plot_vI_SHM = n.zeros(len(plot_Erange))
plot_vI_DW  = n.zeros(len(plot_Erange))

plot_dru_DT = n.zeros(len(plot_Erange))
plot_dru_DW = n.zeros(len(plot_Erange))

for i in n.arange(len(plot_Erange)):
	plot_ff_Helm[i]   = Xe_target_DT.FormFactor(plot_Erange[i])
	plot_ff_HelmDW[i] = Xe_target_DW.FormFactor(plot_Erange[i])

	vMin_ms  = MinimumVelocity_ms(plot_Erange[i], Xe_target_DT, DM_DT)
	vMin_kms = vMin_ms / 1.0e3
	plot_vMins[i] = vMin_kms

	plot_vI_SHM[i] = DM_DT.HaloModel.GetHaloIntegral_ms(vMin_ms)
	plot_vI_DW[i]  = DM_DW.HaloModel.GetHaloIntegral_ms(vMin_ms)

	plot_dru_DT[i] = DifferentialRate(plot_Erange[i], Xe_target_DT, DM_DT) * 3600. * 24.
	plot_dru_DW[i] = DifferentialRate(plot_Erange[i], Xe_target_DW, DM_DW) * 3600. * 24.

## == Plot form factors as function of recoil energy
pyp.figure()
ax1 = pyp.gca()

# ax1.set_xscale('log')
ax1.set_yscale('log')

pyp.plot(plot_Erange , plot_ff_Helm  , color='m' , label="Helm Form Factor" )
pyp.plot(plot_Erange , plot_ff_HelmDW, color='g' , label="Woodward Form Factor" )

# pyp.xlim(LZ_WIMProi_keV)
pyp.xlabel("Recoil Energy [keV]")
pyp.ylabel("$|F(q)|^2$")
pyp.legend(loc='lower left')

## == Plot Velocity integral as function of vmin
pyp.figure()
ax1 = pyp.gca()

pyp.plot(plot_vMins , plot_vI_SHM , color='m' , label="SHM Integral" )
pyp.plot(plot_vMins , plot_vI_DW  , color='g' , label="DW Integral" )

pyp.xlabel("Minimum  velocity [km/s]")
pyp.ylabel("Velocity integral contribution")
pyp.legend(loc='lower left')

## == Plot Velocity integral as function of vmin
pyp.figure()
ax1 = pyp.gca()

ax1.set_yscale('log')

pyp.plot(plot_Erange , plot_dru_DT , color='m' , label="DT Calc" )
pyp.plot(plot_Erange , plot_dru_DW , color='g' , label="DW Calc" )

pyp.xlim([0.1 , 200])
pyp.xlabel("Recoil energy [keV]")
pyp.ylabel("Differential Rate [cts/day/kg/keV]")
pyp.legend(loc='lower left')

pyp.show()
