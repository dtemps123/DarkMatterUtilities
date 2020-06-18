import numpy
from scipy.integrate import quad
from scipy.special import erf, jv

from DarkMatterUtilities.Constants import *

class Interaction():
	InteractionType = "None"

	def __init__(self, _target, _dm, _int_type="SI"):
		self.Target = _target
		self.DarkMatter = _dm
		if ( (_int_type.upper()=="SI") or 
			 (_int_type.upper()=="SD") ):
			self.InteractionType = _int_type.upper()
		else:
			print("Invalid interaction type: "+_int_type.upper()+".")
			print("Please choose SI or SD.")

	def GetDifferentialCrossSection(self, _Er_keV):
		_x_sec_0 = 0.0
		if   (self.InteractionType == "SI"):
			x_sec_0 = self.dSigma_dER_SI(_Er_keV)
		elif (self.InteractionType == "SD"):
			x_sec_0 = self.dSigma_dEr_SD(_Er_keV)
		else:
			x_sec_0 = 0.0

		_q_MeV = numpy.sqrt( 2.0    * self.Target.NuclearMass_GeV * _Er_keV)		# MeV / c
		_frac  = numpy.power(_q_MeV / self.DarkMatter.q_ref_MeV ,2.0)
		_scale = numpy.power(_frac ,  self.DarkMatter.n_param)
		return x_sec_0 * _scale

	def DM_Nucleus_ReducedMass_GeV(self):
		# Returns the target nucelus -- dark matter system reduced mass in GeV/c^2
		_rm_GeV = self.Target.ReducedMass_Nucleus_GeV(self.DarkMatter.Mass_GeV)
		return _rm_GeV

	def DM_Nucleus_ReducedMass_kg(self):
		# Returns the target nucelus -- dark matter system reduced mass in Kg
		return self.DM_Nucleus_ReducedMass_GeV() * GeV_to_kg

	def MinimumVelocity_NoInelasticity_ms(self, _Er_keV):
		# Returns the minimum velocity [m/s] dark matter must have to create a target recoil of E_r [keV]
		# This is old, don't use it, MinimumVelocity_ms takes into account inelasticity
		_mu 	= self.DM_Nucleus_ReducedMass_kg()
		_Er_J 	= _Er_keV * keV_to_J;
		_q    = numpy.sqrt(2.0 * self.Target.NuclearMass_kg * _Er_J)		# kg x m / s
		_frac = _q / (2.0*_mu)												# m / s 
		return _frac

	def MinimumVelocity_ms(self, _Er_keV):
		_mu_GeV    = self.DM_Nucleus_ReducedMass_GeV()
		_M_N_GeV   = self.Target.NuclearMass_GeV
		_M_N_kg    = _M_N_GeV * GeV_to_kg
		_delta_keV = self.DarkMatter.delta_keV
		_fac1_keV  = _delta_keV + ( _M_N_GeV * _Er_keV / _mu_GeV )
		_fac1_J    = _fac1_keV * keV_to_J
		_Er_J      = _Er_keV * keV_to_J
		_fac2      = 1./numpy.sqrt(2.0 * _M_N_kg * _Er_J)
		return _fac1_J * _fac2

	def GetParameter_r(self):
		_num = 4.0 * self.DarkMatter.Mass_GeV * self.Target.NuclearMass_GeV
		_den = numpy.power(self.DarkMatter.Mass_GeV + self.Target.NuclearMass_GeV, 2.0)
		return _num / _den

	def SI_dSigma_dErHalo_cm2_keV(_Er_keV):
		_sigma_p_cm2 = self.DarkMatter.Sigma_cm2
		_mu_p_GeV    = self.DarkMatter.Rmass_DM_proton_GeV
		_mu_N_GeV    = self.DM_Nucleus_ReducedMass_GeV()
		_A           = self.Target.A

		_FFq2  = self.Target.FormFactor.EvaluateFormFactorSquared(_Er_keV, _renormalize=False)
		_r_val = GetParameter_r()

		_sigma_SI = ( 2.0 * _FFq2 * numpy.power(_A,2.) * 
			         _sigma_p_cm2 * numpy.power(_mu_N/_mu_p,2.) / 
			         (1e6 * self.DarkMatter.Mass_GeV * _r_val) )

		return _sigma_SI

	def SI_dSigma_dEr_cm2_keV(self, _Er_keV, _incoming_E_keV):
		_sigma_SI = self.SI_dSigma_dErHalo_cm2_keV(_Er_keV)
		_r_val    = self.GetParameter_r()
		
		if (_Er_keV > (_r_val * _incoming_E_keV)):
			return 0

		_prefac = 0.5 * (self.DarkMatter.Mass_GeV*1e6) / _incoming_E_keV
		return _prefac * _sigma_SI

	def dSigma_dER_SI(self, _Er_keV):
		_sigma_p_cm2 = self.DarkMatter.Sigma_cm2
		_mu_p_GeV    = self.DarkMatter.Rmass_DM_proton_GeV
		_mu_N_GeV    = self.DM_Nucleus_ReducedMass_GeV()
		_A           = self.Target.A

		_FF_norm = True
		_FFq2  = self.Target.FormFactor.EvaluateFormFactorSquared(_Er_keV, _renormalize=_FF_norm)

		# get the total nuclear coupling
		_sigma_A_cm2   = _sigma_p_cm2 * numpy.power(_A,2.) * numpy.power(_mu_N_GeV/_mu_p_GeV,2)

		# calculate coupling term
		_cpl_term  = 0.5 * _sigma_A_cm2 / numpy.power(_mu_N_GeV,2)

		# _E_max = self.Target.RecoilEnergyMax_DM_keV(self.DarkMatter.Mass_GeV)
		# if ( _Er_keV > _E_max ):
		# 	return 0.
		return _cpl_term * _FFq2

	def dSigma_dEr_SD(self, _Er_keV):
		_temp = 1
		return 0.0 

	def DifferentialRate(self, _Er_keV):
		# get conversion factors
		_rho0_conv_fac    = GeV_to_J * 1.e6 / numpy.power(c_ms,2)
		_sigma_A_conv_fac = 1.0e-4
		_dm_mass_conv_fac = 1.e6/numpy.power(c_ms,2)
		_mu_A_conv_fac    = numpy.power(GeV_to_J,2)/numpy.power(numpy.power(c_ms,2),2)
		_total_conv_fac   = ( _rho0_conv_fac
							* _sigma_A_conv_fac
							/(_dm_mass_conv_fac*_mu_A_conv_fac) )

		# get DM number density
		_DM_num_dens_cm3  = self.DarkMatter.Rho0 / self.DarkMatter.Mass_GeV     # cm^-3

		# do velocity integral
		_vmin 			= self.MinimumVelocity_ms(_Er_keV)						# (m/s)
		_vel_integral   = self.DarkMatter.HaloModel.GetHaloIntegral_ms(_vmin)	# s / m

		# get product of terms in current units
		_sigma = self.GetDifferentialCrossSection(_Er_keV)
		_unscaled_prod  = _DM_num_dens_cm3 * _sigma * _vel_integral

		_dru_keV_kg_sec	= _total_conv_fac * _unscaled_prod
		_dru_keV_kg_day = _dru_keV_kg_sec * (24.*3600.)
		return _dru_keV_kg_day		# units: keV^{-1} x kg^{-1} x day^{-1}

	def IntegratedRate(self, _threshold_E_keV):
		# Integrate the differential rate from threshold up to the maximum energy a DM particle can deposit
		_maxE_keV		= self.Target.RecoilEnergyMax_DM_keV(self.Dark.Mass_GeV)								# keV
		_rate			= quad( DifferentialRate, _threshold_E_keV, _maxE_keV)	# cts / day / kg
		return _rate[0]	# units: kg^{-1} x day^{-1}

	def TruncatedIntegratedRate(self, _threshold_E_keV, _max_E_keV):
		# Integrate the differential rate from threshold up to a specified maximum energy
		_rate			= quad( DifferentialRate, _threshold_E_keV, _max_E_keV)	# cts / day / kg
		return _rate[0] # units: kg^{-1} x day^{-1}
