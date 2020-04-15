import numpy as n
import matplotlib.pyplot as pyp

from DMClasses      import *
from DMInteractions import *
from DMPlots        import *

xenon_100kg		= Target( 131.293 , 54.0 , 100.0 , 1.0 , "Xe", 4.7808 )
germanium_100kg	= Target(  72.630 , 32.0 , 100.0 , 1.0 , "Ge", 4.0576 )
argon_100kg		= Target(  39.948 , 18.0 , 100.0 , 1.0 , "Ar", 3.4274 )
silicon_100kg	= Target(  28.085 , 14.0 , 100.0 , 1.0 , "Si", 3.1224 )
neon_100kg		= Target(  20.180 , 10.0 , 100.0 , 1.0 , "Ne", 3.0055 )
helium_100kg	= Target(   4.0   ,  2.0 , 100.0 , 1.0 , "He", 1.6755 )
hydrogen_100kg	= Target(   1.0   ,  1.0 , 100.0 , 1.0 , "H" , 0.8783 )
electron_100kg  = Target( 0.00054858,1.0 , 100.0 , 1.0 , "e-", 2.81794092 )

xenon_10kg		= Target( 131.293 , 54.0 , 10.0 , 1.0 , "Xe", 4.7808 )
germanium_10kg	= Target(  72.630 , 32.0 , 10.0 , 1.0 , "Ge", 4.0576 )
argon_10kg		= Target(  39.948 , 18.0 , 10.0 , 1.0 , "Ar", 3.4274 )
silicon_10kg	= Target(  28.085 , 14.0 , 10.0 , 1.0 , "Si", 3.1224 )
neon_10kg		= Target(  20.180 , 10.0 , 10.0 , 1.0 , "Ne", 3.0055 )
helium_10kg		= Target(   4.0   ,  2.0 , 10.0 , 1.0 , "He", 1.6755 )
hydrogen_10kg	= Target(   1.0   ,  1.0 , 10.0 , 1.0 , "H" , 0.8783 )

WIMP_100MeV		= DarkMatter(   0.1 , 1e-42)

WIMP_1GeV 		= DarkMatter(   1.0 , 1e-45 )
WIMP_8GeV 		= DarkMatter(   8.0 , 1e-45 )
WIMP_10GeV 		= DarkMatter(  10.0 , 1e-45 )
WIMP_40GeV 		= DarkMatter(  40.0 , 1e-46 )
WIMP_80GeV 		= DarkMatter(  80.0 , 1e-46 )
WIMP_100GeV 	= DarkMatter( 100.0 , 1e-46 )

# Xenon recoils, compare WIMP masses
fig1 = WIMP_8GeV.Plot_VelocityDist()
Plot_Overlay_DM_Vmin( fig1, 3.0, xenon_100kg, WIMP_8GeV , 'r-' )
Plot_Overlay_DM_Vmin( fig1, 3.0, xenon_100kg, WIMP_80GeV, 'g-' )
fig1.gca()
pyp.title( "Xenon target, 3.0 keV threshold" )
pyp.legend(loc = 'upper right')
pyp.tight_layout()

# Fixed WIMP mass, compare targets
fig2 = WIMP_8GeV.Plot_VelocityDist()
Plot_Overlay_DM_Vmin( fig2, 3.0,    xenon_100kg, WIMP_8GeV, 'r-' )
Plot_Overlay_DM_Vmin( fig2, 3.0,    argon_100kg, WIMP_8GeV, 'g-' )
Plot_Overlay_DM_Vmin( fig2, 3.0,     neon_100kg, WIMP_8GeV, 'c-' )
Plot_Overlay_DM_Vmin( fig2, 3.0,   helium_100kg, WIMP_8GeV, 'k-' )
Plot_Overlay_DM_Vmin( fig2, 3.0, hydrogen_100kg, WIMP_8GeV, 'y-' )
fig2.gca()
pyp.title( "Various targets, 3.0 keV threshold" )
pyp.legend(loc = 'upper right')
pyp.tight_layout()

# Heavy WIMP, 100 kg detector rates
fig3, thresh_array = Plot_DM_Rate_axes(40.0)
Plot_Overlay_Detected_Rate( fig3, thresh_array,     xenon_100kg, WIMP_100GeV, 'm-' )
Plot_Overlay_Detected_Rate( fig3, thresh_array, germanium_100kg, WIMP_100GeV, 'r-' )
Plot_Overlay_Detected_Rate( fig3, thresh_array,     argon_100kg, WIMP_100GeV, 'y-' )
Plot_Overlay_Detected_Rate( fig3, thresh_array,   silicon_100kg, WIMP_100GeV, 'g-' )
Plot_Overlay_Detected_Rate( fig3, thresh_array,      neon_100kg, WIMP_100GeV, 'b-' )
# Plot_Overlay_Detected_Rate( fig3, thresh_array,    helium_100kg, WIMP_100GeV, 'k-' )
fig3.gca()
pyp.title( r"Rate (cts/100kg/yr) for $10^{-46}$ cm$^{2}$, 100 GeV WIMP" )
pyp.legend(loc = 'upper right')
pyp.tight_layout()

# Light WIMP, 10 kg detector rates
fig4, thresh_array = Plot_DM_Rate_axes(20.0)
Plot_Overlay_Detected_Rate( fig4, thresh_array,     xenon_10kg, WIMP_10GeV, 'm-' )
Plot_Overlay_Detected_Rate( fig4, thresh_array, germanium_10kg, WIMP_10GeV, 'r-' )
Plot_Overlay_Detected_Rate( fig4, thresh_array,     argon_10kg, WIMP_10GeV, 'y-' )
Plot_Overlay_Detected_Rate( fig4, thresh_array,   silicon_10kg, WIMP_10GeV, 'g-' )
Plot_Overlay_Detected_Rate( fig4, thresh_array,      neon_10kg, WIMP_10GeV, 'b-' )
# Plot_Overlay_Detected_Rate( fig4, thresh_array,    helium_10kg, WIMP_10GeV, 'k-' )
fig4.gca()
pyp.title( r"Rate (cts/10kg/yr) for $10^{-45}$ cm$^{2}$, 10 GeV WIMP" )
pyp.legend(loc = 'upper right')
pyp.tight_layout()

# Max recoil energy, various targets
fig5, dm_mass_array = Plot_MaxRecoilE_axes()
Plot_Overlay_TargetRecoils( fig5,  hydrogen_100kg , dm_mass_array, 'm-')
Plot_Overlay_TargetRecoils( fig5,    helium_100kg , dm_mass_array, 'b-')
Plot_Overlay_TargetRecoils( fig5, germanium_100kg , dm_mass_array, 'g-')
Plot_Overlay_TargetRecoils( fig5,     xenon_100kg , dm_mass_array, 'r-')
Plot_Overlay_TargetRecoils( fig5,  electron_100kg , dm_mass_array, 'c-')
fig5.gca()
pyp.title( "Maximum recoil energy, by target" )
pyp.legend(loc = 'upper left')
pyp.tight_layout()

# Minimum velocity for given DM, various targets
fig6 = Plot_MinVelocity_axes(3.0, WIMP_8GeV)
Plot_Overlay_MinVelocity(fig6, 3.0,    xenon_100kg, WIMP_8GeV, 'r.')
Plot_Overlay_MinVelocity(fig6, 3.0,    argon_100kg, WIMP_8GeV, 'b.')
Plot_Overlay_MinVelocity(fig6, 3.0,   helium_100kg, WIMP_8GeV, 'g.')
Plot_Overlay_MinVelocity(fig6, 3.0, hydrogen_100kg, WIMP_8GeV, 'c.')
fig6.gca()
pyp.legend(loc = 'lower right')
pyp.tight_layout()

# Minimum velocity for given DM, various targets
fig7 = Plot_MinVelocity_axes(0.1, WIMP_1GeV)
Plot_Overlay_MinVelocity(fig6, 0.1,    xenon_100kg, WIMP_1GeV, 'r.')
Plot_Overlay_MinVelocity(fig6, 0.1,    argon_100kg, WIMP_1GeV, 'b.')
Plot_Overlay_MinVelocity(fig6, 0.1,   helium_100kg, WIMP_1GeV, 'g.')
Plot_Overlay_MinVelocity(fig6, 0.1, hydrogen_100kg, WIMP_1GeV, 'c.')
fig7.gca()
pyp.legend(loc = 'lower right')
pyp.tight_layout()

# Minimum velocity for given DM, various targets
fig8 = Plot_MinVelocity_axes(3.0, WIMP_100GeV)
Plot_Overlay_MinVelocity(fig6, 3.0,    xenon_100kg, WIMP_100GeV, 'r.')
Plot_Overlay_MinVelocity(fig6, 3.0,    argon_100kg, WIMP_100GeV, 'b.')
Plot_Overlay_MinVelocity(fig6, 3.0,   helium_100kg, WIMP_100GeV, 'g.')
Plot_Overlay_MinVelocity(fig6, 3.0, hydrogen_100kg, WIMP_100GeV, 'c.')
fig8.gca()
pyp.legend(loc = 'upper right')
pyp.tight_layout()

# Lindhard factors
fig9, energies = Plot_LindhardFactor_axes()
Plot_Overlay_LindhardFactor(fig9, energies,    xenon_100kg, 'm-')
Plot_Overlay_LindhardFactor(fig9, energies,   helium_100kg, 'b-')
Plot_Overlay_LindhardFactor(fig9, energies, hydrogen_100kg, 'c-')
fig9.gca()
pyp.legend(loc = 'upper left')
pyp.tight_layout()

# Neutron scattering energy-angle dependence
fig10, theta = Plot_RecoilAngleDist_axes()
Plot_RecoilAngleDist(fig10,    xenon_100kg, theta, 'b-', 'b:')
Plot_RecoilAngleDist(fig10,   helium_100kg, theta, 'g-', 'g:')
Plot_RecoilAngleDist(fig10, hydrogen_100kg, theta, 'r-', 'r:')
fig10.gca()
pyp.legend(loc = 'upper right')
pyp.tight_layout()

pyp.show()