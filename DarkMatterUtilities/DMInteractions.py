from DMClasses import *
from scipy.integrate import quad

def DM_Nucleus_ReducedMass_kg(_target, _dm):
	# Returns the target nucelus -- dark matter system reduced mass in Kg
	_rm_GeV = _target.ReducedMass_Nucleus_GeV(_dm.Mass)
	return _rm_GeV * GeV_to_kg

def MinimumVelocity_ms(_Er_keV, _target, _dm):
	# Returns the minimum velocity [m/s] dark matter must have to create a target recoil of E_r [keV]
	_mu 	= DM_Nucleus_ReducedMass_kg(_target, _dm)
	_Er_J 	= _Er_keV * keV_to_J;
	_q    = n.sqrt(2.0 * _target.NuclearMass_kg * _Er_J)		# kg x m / s
	_frac = _q / (2.0*_mu)										# m / s 
	return _frac

def DifferentialRate(_Er_keV, _target, _dm):
	# Given a dark matter model and a detector model, find the differential rate as function of the recoil energy

	# get conversion factors
	_rho0_conv_fac    = GeV_to_J * 1.e6 / n.power(c_ms,2)
	_sigma_A_conv_fac = 1.0e-4
	_dm_mass_conv_fac = 1.e6/n.power(c_ms,2)
	_mu_A_conv_fac    = n.power(GeV_to_J,2)/n.power(n.power(c_ms,2),2)
	_total_conv_fac   = ( _rho0_conv_fac
						* _sigma_A_conv_fac
						/(_dm_mass_conv_fac*_mu_A_conv_fac) )

	# get DM number density
	_DM_num_dens_cm3  = _dm.Rho0 / _dm.Mass        # cm^-3
	
	# get reduced masses
	_mu_N_GeV = _target.ReducedMass_Nucleus_GeV(_dm.Mass)	# GeV
	_mu_n_GeV = _dm.Rmass_DM_proton_GeV							# GeV

	# get the total nuclear coupling
	_sigma_A_cm2   = _dm.Sigma * (_target.A**2) * n.power(_mu_N_GeV/_mu_n_GeV,2)

	# calculate coupling term
	_cpl_term  = 0.5 * _sigma_A_cm2 / n.power(_mu_N_GeV,2)	# m^2 / kg^2

	# get form factor
	_formfactor		= _target.FormFactor(_Er_keV)			# dimensionless

	# do velocity integral
	_vmin 			= MinimumVelocity_ms(_Er_keV, _target, _dm)								# (m/s)
	_vel_integral   = _dm.HaloModel.GetHaloIntegral_ms(_vmin)	# s / m

	# get product of terms in current units
	_unscaled_prod  = _DM_num_dens_cm3 * _cpl_term * _formfactor * _vel_integral
	# this has units of s / kg^2 / m^2

	_dru			= _total_conv_fac * _unscaled_prod
	return _dru

def IntegratedRate(_threshold_E_keV, _target, _dm):
	# Integrate the differential rate from threshold up to the maximum energy a DM particle can deposit
	_maxE 			= _target.RecoilEnergyMax_DM_keV(_dm.Mass)								# keV
	_rate			= quad( DifferentialRate, _threshold_E_keV, _maxE, args=(_target, _dm))	# Hz / kg / keV
	return _rate[0] * (365.25 * 24. * 3600.) * _target.TotalMass 							# Cts / total mass / year\\

def TruncatedIntegratedRate(_threshold_E_keV, _max_E_keV, _target, _dm):
	# Integrate the differential rate from threshold up to a specified maximum energy
	_rate			= quad( DifferentialRate, _threshold_E_keV, _max_E_keV, args=(_target, _dm))	# Hz / kg / keV
	return _rate[0] * (365.25 * 24. * 3600.) * _target.TotalMass 									# Cts / total mass / year
