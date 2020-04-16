import numpy as n
import matplotlib.pyplot as pyp

from DMClasses      import *
from DMInteractions import *
from DMPlots        import *

pyp.rcParams['axes.grid'] = True
pyp.rcParams.update({'font.size': 16})
pyp.rc('text', usetex=True)
pyp.rc('font', family='serif')

LZ_target = Target( 131.293 , 54.0 , 7.0e3 , 1.0 , "Xe", 4.7808 )
LZ_exposure_yr = 1000 / 365.25
LZ_WIMProi_keV = n.array([1.5 , 7.5])

# plot_mass_DM  = 50.0 # GeV/c^2
# plot_xsec_DM  = 1.0e-45 # cm^2
# plot_DM       = DarkMatter(plot_mass_DM, plot_xsec_DM)

# get a recoil energy array [keV]
plot_Erange = n.logspace(start=-1.0, stop=3.0, num=1000)
plot_qrn    = n.sqrt(2.0 * LZ_target.NuclearMass_GeV * plot_Erange) * ( LZ_target.FF_Rn / hbarc_MeV_fm) 

plot_ff_type0 = n.zeros(len(plot_Erange))
plot_ff_type1 = n.zeros(len(plot_Erange))
plot_ff_type2 = n.zeros(len(plot_Erange))
plot_ff_type3 = n.zeros(len(plot_Erange))
# plot_ff_type4 = n.zeros(len(plot_Erange))
plot_ff_Helm  = n.zeros(len(plot_Erange))
for i in n.arange(len(plot_Erange)):
	LZ_target.FF_type = 0
	plot_ff_type0[i] = LZ_target.FormFactor(plot_Erange[i])

	LZ_target.FF_type = 1
	plot_ff_type1[i] = LZ_target.FormFactor(plot_Erange[i])

	LZ_target.FF_type = 2
	plot_ff_type2[i] = LZ_target.FormFactor(plot_Erange[i])

	LZ_target.FF_type = 3
	plot_ff_type3[i] = LZ_target.FormFactor(plot_Erange[i])

	# LZ_target.FF_type = 4
	# plot_ff_type4[i] = LZ_target.FormFactor(plot_Erange[i])

	plot_ff_Helm[i]  = LZ_target.HelmFormFactor(plot_Erange[i])

## == Plot form factors as function of recoil energy
pyp.figure()
ax1 = pyp.gca()

ax1.set_xscale('log')
ax1.set_yscale('log')

pyp.plot(plot_Erange , plot_ff_type0 , color='r' , label="(0) Thin shell: exp[-$(q r_n)^{2/3} / 3$]" )
pyp.plot(plot_Erange , plot_ff_type1 , color='g' , label="(1) Thin shell: $[ \sin (q r_n) / (q r_n) ]^2$" )
pyp.plot(plot_Erange , plot_ff_type2 , color='b' , label="(2) Solid sphere: exp[-$(q r_n)^{2/3} / 5$]" )
pyp.plot(plot_Erange , plot_ff_type3 , color='c' , label="(3) Solid sphere: $\{ 3 [ \sin(q r_n) - q r_n \cos(q r_n)] / (q r_n)^3 \}^2$" )
# pyp.plot(plot_Erange , plot_ff_type4 , color='k' , label="(4) Two-parameter Fermi" )
pyp.plot(plot_Erange , plot_ff_Helm  , color='m' , label="Helm Form Factor" )

# pyp.xlim(LZ_WIMProi_keV)
pyp.xlabel("Recoil Energy [keV]")
pyp.ylabel("$|F(q)|^2$")
pyp.title("Helm Form Factor")
pyp.legend(loc='lower left')

## == Plot form factors as function of q*Rn to compare with Lewin&Smith
pyp.figure()
ax1 = pyp.gca()

ax1.set_yscale('log')

# pyp.plot(plot_qrn , plot_ff_type0 , 'r--' , label="(0) Thin shell: exp[-$(q r_n)^{2/3} / 3$]" )
pyp.plot(plot_qrn , plot_ff_type1 , 'r:'  , label="(1) Thin shell: $[ \sin (q r_n) / (q r_n) ]^2$" )
# pyp.plot(plot_qrn , plot_ff_type2 , 'b--' , label="(2) Solid sphere: exp[-$(q r_n)^{2/3} / 5$]" )
pyp.plot(plot_qrn , plot_ff_type3 , 'b:'  , label="(3) Solid sphere: $\{ 3 [ \sin(q r_n) - q r_n \cos(q r_n)] / (q r_n)^3 \}^2$" )
# pyp.plot(plot_qrn , plot_ff_type4 , 'g-'  , label="(4) Two-parameter Fermi" )
pyp.plot(plot_qrn , plot_ff_Helm  , color='m' , label="Helm Form Factor" )

pyp.xlim([0.0   , 10.0])
pyp.ylim([1.e-4 , 1.e0])

pyp.xlabel("$qr_n$")
pyp.ylabel("$|F(q)|^2$")
pyp.title("Form Factors")
pyp.legend(loc='lower left')

pyp.show()

