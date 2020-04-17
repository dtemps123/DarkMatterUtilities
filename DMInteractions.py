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

def DifferentialRate_old(_Er_keV, _target, _dm):
	# Given a dark matter model and a detector model, find the differential rate as function of the recoil energy
	_DM_num_dens 	= _dm.Rho0 / _dm.Mass 													# cm^-3
	_coupling		= 0.5 * _dm.Sigma * (_target.A**2) / (_dm.Rmass_DM_proton**2) 			# cm^2  x  kg^-2
	_formfactor		= _target.FormFactor(_Er_keV)										# dimensionless
	_vmin 			= MinimumVelocity_ms(_Er_keV, _target, _dm)								# (m/s)
	_vel_integral   = _dm.HaloModel.GetHaloIntegral_ms(_vmin)												# m^-1  x  s
	_unitfactors	= 10. * c_ms**2 / kg_to_kev												# cm  x  m^-1  x  kg  x  keV^-1  x  m^2  x  s^-2
	_dru			= _DM_num_dens * _coupling * _formfactor * _vel_integral * _unitfactors	# Hz / kg / keV
	return _dru

def DifferentialRate(_Er_keV, _target, _dm):
	# Given a dark matter model and a detector model, find the differential rate as function of the recoil energy

	# get rho_0 in J/m^3
	_DM_rho_J_cm3 = (_dm.Rho0 * 1e6) * keV_to_J		# J / cm^3
	_DM_rho_J_m3  = _DM_rho_J_cm3 * 1.0e3			# J / m^3

	# get DM mass in kg
	_DM_mass_kg   = _dm.Mass * GeV_to_kg			# kg
	_DM_mass_J    = _DM_mass_kg * c_ms * c_ms		# J

	# get DM number density
	_DM_num_dens  =_DM_rho_J_m3 / _DM_mass_J        # m^-3
	
	# get reduced masses
	_mu_N_GeV = _target.ReducedMass_Nucleus_GeV(_dm.Mass)	# GeV
	_mu_n_GeV = _dm.Rmass_DM_proton							# GeV
	_mu_N_kg  = DM_Nucleus_ReducedMass_kg(_target, _dm)		# kg

	# get the total nuclear coupling
	_sigma_A_cm2   = _dm.Sigma * (_target.A**2) * n.power(_mu_N_GeV/_mu_n_GeV,2)
	_sigma_A_m2    = _sigma_A_cm2 / 100.					# m^2

	# calculate coupling term
	_cpl_term  = 0.5 * _sigma_A_m2 / n.power(_mu_N_kg,2)	# m^2 / kg^2

	# get form factor
	_formfactor		= _target.FormFactor(_Er_keV)			# dimensionless

	# do velocity integral
	_vmin 			= MinimumVelocity_ms(_Er_keV, _target, _dm)								# (m/s)
	_vel_integral   = _dm.HaloModel.GetHaloIntegral_ms(_vmin)	# s / m

	# get product of terms in current units
	_unscaled_prod  = _DM_num_dens * _cpl_term * _formfactor * _vel_integral
	# this has units of s / kg^2 / m^2

	_scale_facs     = c_ms * c_ms / kg_to_kev
	_dru			= _scale_facs * _unscaled_prod
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
