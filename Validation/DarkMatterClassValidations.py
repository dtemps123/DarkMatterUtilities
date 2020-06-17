import numpy as n
import matplotlib.pyplot as pyp

## Import the modules
import DarkMatterUtilities.Constants as dmcon
import DarkMatterUtilities.DarkMatter as dmdm

## Initialize a dark matter model with explicit Halo Model
WIMP_50GeV_Halo0 = dmdm.DarkMatter(50.0 , 1.0e-45 , _HaloModelType=0)
WIMP_50GeV_Halo1 = dmdm.DarkMatter(50.0 , 1.0e-45 , _HaloModelType=1)

## Create a plot of the velocity
## distribution
## ==============================
Vel_range_kms = n.linspace(start=0.0, stop=1000.0, num=1000)
Vel_range_ms  = Vel_range_kms * 1e3

PDF_val_0_ms  = n.zeros(len(Vel_range_ms))
PDF_val_1_ms  = n.zeros(len(Vel_range_ms))

for i in n.arange(len(Vel_range_ms)):
	PDF_val_0_ms[i] = WIMP_50GeV_Halo0.HaloModel.GetHaloPDF_ms(Vel_range_ms[i])
	PDF_val_1_ms[i] = WIMP_50GeV_Halo1.HaloModel.GetHaloPDF_ms(Vel_range_ms[i])

pyp.figure()
pyp.plot(Vel_range_kms, PDF_val_0_ms, label=WIMP_50GeV_Halo0.HaloModel.ModelName)
pyp.plot(Vel_range_kms, PDF_val_1_ms, label=WIMP_50GeV_Halo1.HaloModel.ModelName)

pyp.legend(loc='lower right')
pyp.xlabel("DM Velocity [km/s]")
pyp.ylabel("Halo Model PDF")
pyp.title(str(WIMP_50GeV_Halo0.Mass_GeV)+" GeV"+r"$/c^2$"+" WIMP")


## Create a plot of the velocity
## distribution integral
## ==============================
Vel_range_kms = n.linspace(start=0.0, stop=1000.0, num=1000)
Vel_range_ms  = Vel_range_kms * 1e3

CDF_val_0_ms  = n.zeros(len(Vel_range_ms))
CDF_val_1_ms  = n.zeros(len(Vel_range_ms))

for i in n.arange(len(Vel_range_ms)):
	CDF_val_0_ms[i] = WIMP_50GeV_Halo0.HaloModel.GetHaloIntegral_ms(Vel_range_ms[i])
	CDF_val_1_ms[i] = WIMP_50GeV_Halo1.HaloModel.GetHaloIntegral_ms(Vel_range_ms[i])

pyp.figure()
pyp.plot(Vel_range_kms, CDF_val_0_ms, label=WIMP_50GeV_Halo0.HaloModel.ModelName)
pyp.plot(Vel_range_kms, CDF_val_1_ms, label=WIMP_50GeV_Halo1.HaloModel.ModelName)

pyp.legend(loc='lower right')
pyp.xlabel("DM Velocity [km/s]")
pyp.ylabel("Halo Model CDF")
pyp.title(str(WIMP_50GeV_Halo0.Mass_GeV)+" GeV"+r"$/c^2$"+" WIMP")


## ==============================
pyp.show()
