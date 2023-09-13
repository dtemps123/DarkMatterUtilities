import numpy as n
import matplotlib.pyplot as pyp

## Import the modules
import DarkMatterUtilities.Constants   as dmcon
import DarkMatterUtilities.Targets     as dmt
import DarkMatterUtilities.FormFactors as dmff

## ==================================================== ##
## Define the plotting style options

## Turn the grid on for all plots
# pyp.rcParams['axes.grid'] = True

## Set the global font size for the plots
pyp.rcParams.update({'font.size': 18})

## Use the LaTeX engine to draw text
pyp.rc('text', usetex=True)

## Select the typeface
pyp.rc('font', family='serif')

## Get a list of the default colors in order
default_colors = pyp.rcParams['axes.prop_cycle'].by_key()['color']
## ==================================================== ##

## Initialize a simple xenon target
avg_xe_target = dmt.SimpleTarget( 131.293 , 54.0 , "Xe", 4.7808 )

## Initialize the form factors
# ff_unity   = dmff.FormFactor(avg_xe_target, 0)
ff_helm    = dmff.FormFactor(avg_xe_target, 1)
# ff_shell1  = dmff.FormFactor(avg_xe_target, 2)
ff_shell2  = dmff.FormFactor(avg_xe_target, 3)
# ff_sphere1 = dmff.FormFactor(avg_xe_target, 4)
ff_sphere2 = dmff.FormFactor(avg_xe_target, 5)
ff_tophat  = dmff.FormFactor(avg_xe_target, 6)
# ff_helmdw  = dmff.FormFactor(avg_xe_target, 7)
# ff_helm00  = dmff.FormFactor(avg_xe_target, 8)

## Create a plot of the form 
## factors against recoil energy
## ==============================
E_recoil_range = n.linspace(start=0.0, stop=700.0, num=100000)

# FF_vals_unity   = n.zeros(len(E_recoil_range))
FF_vals_helm    = n.zeros(len(E_recoil_range))
# FF_vals_shell1  = n.zeros(len(E_recoil_range))
FF_vals_shell2  = n.zeros(len(E_recoil_range))
# FF_vals_sphere1 = n.zeros(len(E_recoil_range))
FF_vals_sphere2 = n.zeros(len(E_recoil_range))
FF_vals_tophat  = n.zeros(len(E_recoil_range))
# FF_vals_helmdw  = n.zeros(len(E_recoil_range))
# FF_vals_helm00  = n.zeros(len(E_recoil_range))

## Do we want to renormalize the form factor such that at q=0, |F(q)|^2=1
norm = True

for i in n.arange(len(E_recoil_range)):
	# FF_vals_unity[i]   = ff_unity.EvaluateFormFactorSquared(E_recoil_range[i], _renormalize=norm)
	FF_vals_helm[i]    = ff_helm.EvaluateFormFactorSquared(E_recoil_range[i], _renormalize=norm)
	# FF_vals_shell1[i]  = ff_shell1.EvaluateFormFactorSquared(E_recoil_range[i], _renormalize=norm)
	FF_vals_shell2[i]  = ff_shell2.EvaluateFormFactorSquared(E_recoil_range[i], _renormalize=norm)
	# FF_vals_sphere1[i] = ff_sphere1.EvaluateFormFactorSquared(E_recoil_range[i], _renormalize=norm)
	FF_vals_sphere2[i] = ff_sphere2.EvaluateFormFactorSquared(E_recoil_range[i], _renormalize=norm)
	FF_vals_tophat[i]  = ff_tophat.EvaluateFormFactorSquared(E_recoil_range[i], _renormalize=norm)
	# FF_vals_helmdw[i]  = ff_helmdw.EvaluateFormFactorSquared(E_recoil_range[i], _renormalize=norm)
	# FF_vals_helm00[i]  = ff_helm00.EvaluateFormFactorSquared(E_recoil_range[i], _renormalize=norm)
	
# pyp.figure()
# pyp.plot(E_recoil_range, FF_vals_unity   , label=ff_unity.FF_Name)
# pyp.plot(E_recoil_range, FF_vals_helm    , label=ff_helm.FF_Name)
# pyp.plot(E_recoil_range, FF_vals_shell1  , label=ff_shell1.FF_Name)
# pyp.plot(E_recoil_range, FF_vals_shell2  , label=ff_shell2.FF_Name)
# pyp.plot(E_recoil_range, FF_vals_sphere1 , label=ff_sphere1.FF_Name)
# pyp.plot(E_recoil_range, FF_vals_sphere2 , label=ff_sphere2.FF_Name)
# pyp.plot(E_recoil_range, FF_vals_tophat  , label=ff_tophat.FF_Name)
# pyp.plot(E_recoil_range, FF_vals_helmdw  , label=ff_helmdw.FF_Name)
# pyp.plot(E_recoil_range, FF_vals_helm00  , label=ff_helm00.FF_Name)

# pyp.legend(loc='lower right')
# pyp.xlabel("Recoil energy [keV]")
# pyp.ylabel("Form Factor: "+r"$|F(q)|^2$")
# pyp.title("Form Factors")

## Create a plot of the form 
## factors against momentum transfer
## ==============================
qRn_range = n.zeros(len(E_recoil_range))
for i in n.arange(len(E_recoil_range)):
	qRn_range[i] = avg_xe_target.GetQRn(E_recoil_range[i])

pyp.figure()
# pyp.plot(qRn_range, FF_vals_unity   , label=ff_unity.FF_Name)
pyp.plot(qRn_range, FF_vals_helm    , label=ff_helm.FF_Name)
# pyp.plot(qRn_range, FF_vals_shell1  , label=ff_shell1.FF_Name)
pyp.plot(qRn_range, FF_vals_shell2  , label=ff_shell2.FF_Name)
# pyp.plot(qRn_range, FF_vals_sphere1 , label=ff_sphere1.FF_Name)
pyp.plot(qRn_range, FF_vals_sphere2 , label=ff_sphere2.FF_Name)
pyp.plot(qRn_range, FF_vals_tophat  , label=ff_tophat.FF_Name)
# pyp.plot(qRn_range, FF_vals_helmdw  , label=ff_helmdw.FF_Name)
# pyp.plot(qRn_range, FF_vals_helm00  , label=ff_helm00.FF_Name)

pyp.yscale('log')
pyp.ylim([1e-4, 4e0])

pyp.legend(loc='lower left')
pyp.xlabel(r"$q r_n$")
pyp.ylabel("Form Factor: "+r"$|F(q)|^2$")
pyp.title("Form Factors")

## ==============================
pyp.show()
