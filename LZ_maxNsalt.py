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

# Mass points from X1T curve [GeV/c^2]
X1T_mass_vals = n.array([5.9748165212, 6.238511403, 6.3900049633, 6.6082962873, 6.7687695585, 6.9664896395, 7.3089407293, 7.594982817, 8.0838709799, 8.5630387015, 9.0706088686, 9.470934657, 9.9842932724, 10.6269810722, 11.0959963982, 11.9241374584, 12.6309352031, 13.6389274235, 14.4473685521, 15.4513124257, 16.4459116143, 17.1717418896, 17.9296062413, 18.9923753766, 20.3121505894, 21.6196413667, 23.1219846081, 24.6103440752, 26.4471179621, 28.4209780392, 30.5421556279, 32.9795250524, 35.4409262515, 37.9037069757, 40.7326206623, 43.7726681215, 47.2658781034, 50.3083771943, 54.0631037937, 57.2676693206, 61.8378269713, 65.8183207205, 69.3859096881, 74.5644730857, 79.364179569, 87.7787646427, 98.9670387567, 108.936010092, 121.647863831, 140.483662313, 155.378438954, 179.437036223, 206.228829275, 240.457760905, 275.037583294, 301.292893216, 342.9715435, 396.076815347, 450.867178626, 530.769712705, 607.098804041, 687.772066178, 782.913412248, 891.21591472, 985.707058888])
# Cross-section points from X1T curve [cm^2]
X1T_xsec_vals = n.array([2.522521417950E-44, 1.639377231470E-44, 1.272293053350E-44, 9.874052064050E-45, 7.663085474420E-45, 5.947191548850E-45, 4.117929423850E-45, 2.961826626240E-45, 1.924880042670E-45, 1.283089292410E-45, 9.112405290490E-46, 6.982914061860E-46, 5.351066731650E-46, 4.048922364610E-46, 3.264085061650E-46, 2.407971686880E-46, 1.965974958150E-46, 1.487568811140E-46, 1.214517200160E-46, 1.017043314990E-46, 8.516776085410E-47, 7.695529938910E-47, 6.953473996110E-47, 6.203839454820E-47, 5.535020912230E-47, 5.129697355250E-47, 4.754055201180E-47, 4.519036733520E-47, 4.295636490270E-47, 4.135364178480E-47, 4.135364178480E-47, 4.188112592130E-47, 4.241533835310E-47, 4.350429248700E-47, 4.462120398620E-47, 4.635056643380E-47, 4.876108879530E-47, 5.129697355250E-47, 5.396474034240E-47, 5.605622571490E-47, 5.972371151200E-47, 6.363114304030E-47, 6.609726451270E-47, 7.042168656810E-47, 7.408405929110E-47, 8.095746126720E-47, 8.959702262600E-47, 9.790969283150E-47, 1.083583505040E-46, 1.245698220460E-46, 1.361272133280E-46, 1.545221796950E-46, 1.776402126560E-46, 2.068218079720E-46, 2.347697847840E-46, 2.565513545050E-46, 2.912193200220E-46, 3.347885853050E-46, 3.848762398120E-46, 4.538169722100E-46, 5.151415798070E-46, 5.847530248900E-46, 6.637711136530E-46, 7.630777796250E-46, 8.338749304370E-46])
X1T_Npoints = len(X1T_mass_vals)

# Array to store results points in
LZ_Max_N_obs = n.zeros(X1T_Npoints)

# Calculate max # of salt LZ will observe
for i in n.arange(X1T_Npoints):
	this_DM = DarkMatter(X1T_mass_vals[i], X1T_xsec_vals[i])
	this_Nevts_per_yeaR = TruncatedIntegratedRate(LZ_WIMProi_keV[0], LZ_WIMProi_keV[1], LZ_target, this_DM)
	LZ_Max_N_obs[i] = this_Nevts_per_yeaR * LZ_exposure_yr

plot_mass_DM  = 50.0 # GeV/c^2
plot_xsec_DM  = 1.0e-45 # cm^2
plot_DM       = DarkMatter(plot_mass_DM, plot_xsec_DM)
plot_Erange   = n.logspace(start=-1.0, stop=3.0, num=1000)
plot_diffrate = n.zeros(len(plot_Erange))
for i in n.arange(len(plot_Erange)):
	plot_diffrate[i] = DifferentialRate(plot_Erange[i], LZ_target, plot_DM) * 3600. * 24.5

## == Plot differential rate for specific mass
pyp.figure()
ax1 = pyp.gca()

ax1.set_xscale('log')
ax1.set_yscale('log')

pyp.plot(plot_Erange , plot_diffrate ) # , 'ro')

# pyp.xlim(LZ_WIMProi_keV)
pyp.ylim([1e-12,1e-3])
pyp.xlabel("Recoil Energy [keV]")
pyp.ylabel("Differential Rate [day$^{-1}$ kg$^{-1}$ keV$^{-1}$]")
pyp.title("Differential rate for M$_{\chi}=$"+str(plot_mass_DM))

## == Plot X1T Result
pyp.figure()
ax1 = pyp.gca()

ax1.set_xscale('log')
ax1.set_yscale('log')

pyp.scatter(X1T_mass_vals , X1T_xsec_vals ) # , 'ro')

pyp.xlim([1.0,2.0e3])
pyp.xlabel("WIMP Mass [GeV/$c^2$]")
pyp.ylabel("X1T Limit Cross-section [cm$^2$]")

## == Plot LZ Max N DM events
pyp.figure()
ax1 = pyp.gca()

ax1.set_xscale('log')
ax1.set_yscale('log')

pyp.scatter(X1T_mass_vals , LZ_Max_N_obs ) # , 'ro')

pyp.xlim([1.0,2.0e3])
pyp.xlabel("WIMP Mass [GeV/$c^2$]")
pyp.ylabel("Max \# of DM events seen in LZ \n using XENON-1T limit")

pyp.show()

