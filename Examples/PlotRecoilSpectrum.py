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
LZ_target.FF_type = 5

plot_mass_DM  = 10000.0 # GeV/c^2
plot_xsec_DM  = 1.0e-45 # cm^2
plot_DM       = DarkMatter(plot_mass_DM, plot_xsec_DM)
plot_DM.HaloModel.Model = 1

plot_Erange   = n.logspace(start=-1.0, stop=3.0, num=1000)
plot_diffrate_DT = n.zeros(len(plot_Erange))
plot_diffrate_DW = n.zeros(len(plot_Erange))
for i in n.arange(len(plot_Erange)):
	plot_diffrate_DT[i] = DifferentialRate(plot_Erange[i], LZ_target, plot_DM) * 3600. * 24.5
	plot_diffrate_DW[i] = DifferentialRate(plot_Erange[i], LZ_target, plot_DM) * 3600. * 24.5

## == Plot differential rate for specific mass
pyp.figure()
ax1 = pyp.gca()

ax1.set_xscale('log')
ax1.set_yscale('log')

pyp.plot(plot_Erange , plot_diffrate_DT , label="Temples")
pyp.plot(plot_Erange , plot_diffrate_DW , label="Woodward")

pyp.xlabel("Recoil Energy [keV]")
pyp.ylabel("Differential Rate [day$^{-1}$ kg$^{-1}$ keV$^{-1}$]")
pyp.title("Differential rate for M$_{\chi}=$"+str(plot_mass_DM)+" GeV/$c^2$, $\sigma_n=$"+str(plot_xsec_DM)+"cm$^2$")
pyp.legend(loc='lower left')

pyp.show()