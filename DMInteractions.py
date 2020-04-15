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
	# _vmin2 	= 0.5 * _Er_J * _target.NuclearMass_kg  / _mu**2
	# return n.sqrt(_vmin2)

	# Alternate calculation
	_q    = n.sqrt(2.0 * _target.NuclearMass_kg * _Er_J)		# kg x m / s
	_frac = _q / (2.0*_mu)										# m / s 
	return _frac

def DifferentialRate(_Er_keV, _target, _dm):
	# Given a dark matter model and a detector model, find the differential rate as function of the recoil energy
	_DM_num_dens 	= _dm.Rho0 / _dm.Mass 													# cm^-3
	# _coupling		= _dm.Sigma * _target.A / (2.0 * _dm.Rmass_DM_proton**2) 				# cm^2  x  kg^-2
	_coupling		= 0.5 * _dm.Sigma * (_target.A**2) / (_dm.Rmass_DM_proton**2) 			# cm^2  x  kg^-2
	# _formfactor		= _target.FormFactor(_Er_keV)										# dimensionless
	_formfactor		= _target.HelmFormFactor_v2(_Er_keV)									# dimensionless
	_vmin 			= MinimumVelocity_ms(_Er_keV, _target, _dm)								# (m/s)
	if (_vmin > MW_esc_vel_ms):
		return 0
	_MBfactor		= quad(_dm.Velocity_Dist_ms, _vmin , MW_esc_vel_ms)						# m^-1  x  s
	_unitfactors	= 10. * c_ms**2 / kg_to_kev												# cm  x  m^-1  x  kg  x  keV^-1  x  m^2  x  s^-2
	_dru			= _DM_num_dens * _coupling * _formfactor * _MBfactor[0] * _unitfactors	# Hz / kg / keV
	return _dru

def IntegratedRate(_threshold_E_keV, _target, _dm):
	# Integrate the differential rate from threshold up to the maximum energy a DM particle can deposit
	_maxE 			= _target.RecoilEnergyMax_DM_keV(_dm.Mass)								# keV
	_rate			= quad( DifferentialRate, _threshold_E_keV, _maxE, args=(_target, _dm))	# Hz / kg / keV
	return _rate[0] * (365.25 * 24. * 3600.) * _target.TotalMass 							# Cts / total mass / year\\

def TruncatedIntegratedRate(_threshold_E_keV, _max_E_keV, _target, _dm):
	# Integrate the differential rate from threshold up to a specified maximum energy
	_rate			= quad( DifferentialRate, _threshold_E_keV, _max_E_keV, args=(_target, _dm))	# Hz / kg / keV
	return _rate[0] * (365.25 * 24. * 3600.) * _target.TotalMass 									# Cts / total mass / year\\