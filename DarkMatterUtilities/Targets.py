import numpy

import DarkMatterUtilities.FormFactors
from DarkMatterUtilities.Constants import *

## This class defines a general simple target, which
## can be created by a user. The simple target assumes
## a single element, with it's average atomic mass
class SimpleTarget:
	N_isotopes      = 1
	Name 			= "Xenon"			# name of target
	A 				= 1.0				# amu "dimensionless"
	Z 				= 1.0				# amu "dimensionless"
	NuclearMass_GeV = A * amu_to_GeV	# nuclear mass in GeV
	NuclearMass_kg	= A * amu_to_kg		# nuclear mass in kg
	FF_Rn           = 1.0

	def __init__(self, _A, _Z, _Name, _FF_Rn):
		self.A = _A
		self.Z = _Z
		self.NuclearMass_GeV 	= _A * amu_to_GeV
		self.NuclearMass_kg		= _A * amu_to_kg
		self.Name 				= _Name
		self.FF_Rn				= _FF_Rn

	def ReducedMass_Nucleus_amu(self, _mass_amu):
		# Returns the reduced mass of the system consisting of target nucleus and a second mass specified in amu
		_numerator   = self.A * _mass_amu
		_denominator = self.A + _mass_amu
		return (_numerator / _denominator) 

	def ReducedMass_Nucleus_GeV(self, _mass_GeV):
		# Returns the reduced mass of the system consisting of target nucleus and a second mass specified in GeV
		_numerator   = self.NuclearMass_GeV * _mass_GeV
		_denominator = self.NuclearMass_GeV + _mass_GeV
		return (_numerator / _denominator) 

	def ReducedMass_Nucleus_kg(self, _mass_kg):
		# Returns the reduced mass of the system consisting of target nucleus and a second mass specified in kg
		_numerator   = self.NuclearMass_kg * _mass_kg
		_denominator = self.NuclearMass_kg + _mass_kg
		return (_numerator / _denominator) 

	def GetMomentumTransfer(self, _Er_keV):
		_q2 = 2.0 * self.NuclearMass_GeV * _Er_keV
		_q  = numpy.sqrt(_q2)
		return _q

	def GetQRn(self, _Er_keV):
		_q   = self.GetMomentumTransfer(_Er_keV)
		_qrn = _q * ( self.FF_Rn / hbarc_MeV_fm) 
		return _qrn

	def LindhardFactor(self, _Er_keV):
		# Determine the lindhard factor for a nuclear recoil of a specified recoil energy in keV
		_Z 	= self.Z
		_e 	= 11.5 * _Er_keV * numpy.power(_Z, -7./3.)
		_k 	= 0.133 * numpy.power(_Z, 2./3.) * numpy.power(self.A, -1./2.)
		_g 	= (3.0 * numpy.power(_e, 0.15)) + (0.7 * numpy.power(_e, 0.6)) + _e
		_LF	= (_k * _g) / ( 1. + (_k * _g))
		return _LF

	def RecoilEnergyAngularDist(self, _theta):
		# Given a fixed angle for the outgoing particle, what is the 
		# energy fraction deposited (this is recoil energy)
		_m1 	= m_neutron_GeV
		_m2 	= self.NuclearMass_GeV
		_Mfrac 	= _m1 * _m1 / numpy.power(_m1+_m2,2.)
		_term1 	= numpy.cos(_theta)
		_rtarg	= (_m2**2./_m1**2.) - numpy.power(numpy.sin(_theta),2.)
		_term2  = numpy.sqrt(_rtarg)
		_sqfac	= numpy.power(_term1 + _term2 , 2)
		return 1. - _Mfrac*_sqfac

	def RecoilEnergyMax_AnyParticle_keV(self, _incoming_mass_GeV, _incoming_E_keV):
		# The maximum energy recoil a particle of a given energy can produce.
		# Calculated from classical 2-body kinematics
		_m1 	= _incoming_mass_GeV
		_m2 	= self.NuclearMass_GeV
		_scale  = 4.0 * _m1 * _m2 / numpy.power(_m1+_m2,2)

		return _incoming_E_keV * _scale

	def RecoilEnergyMax_DM_keV(self, _dm_mass_GeV):
		# Calculates the maximum recoil energy a DM particle of specified mass moving at the Milky Way
		# escape velocity can produce in this detector
		_E_max_dm_keV	= 0.5 * (_dm_mass_GeV * GeV_to_kg) * (MW_esc_vel_ms)**2 * (1./keV_to_J)
		return self.RecoilEnergyMax_AnyParticle_keV(_dm_mass_GeV, _E_max_dm_keV)

	def MinimumDetectableMass_GeV(self, _threshold_keV):
		# For a given recoil energy threshold, what is the smallest mass particle that can produce a recoil 
		# of that energy.
		# Calculated from classical 2-body kinematics
		_threshold_GeV 	= _threshold_keV * 1e-6
		_M_GeV			= self.NuclearMass_GeV
		_beta 			= MW_esc_vel_ms / c_ms
		return numpy.sqrt(1./2. * _threshold_GeV * _M_GeV / _beta**2)

## This class defines a complex target consisting of 
## xenon isotopes in concentrations determined by 
## their natural abundance. Most results from methods
## are calculated as a weighted average over the 
## specified isotopic abundance.
class NaturalXenonTarget:
	## Long name of material
	Name 			  = "Xenon (natural abundance)"
	
	## Number of unique isotopes in the target
	N_isotopes        = 13

	## Array for atomic numbers of target
	Z_array           = 54.0 * numpy.ones(N_isotopes)
	Z_avg             = 0.0   ## calculated on initialization

	## Array for atomic numbers of target
	A_array           = numpy.arange(start=124., stop=137.)
	A_avg             = 0.0   ## calculated on initialization

	## Array for nuclear form factor radius parameter r_n
	FF_Rn_array       = 4.7808 * numpy.ones(N_isotopes)

	## Array for each isotope's abundance in target
	Abundance_array   = numpy.array([0.095 , 0.000 , 0.089 , 0.000 , 1.910 , 26.401 , 4.071 , 21.232 , 26.909 , 0.000 , 10.436 , 0.000 , 8.857 ]) / 100.

	## Create arrays for each isotope target and abundance fraction
	Isotope_list      = numpy.zeros(N_isotopes, dtype=object)
	
	def __init__(self):
		print ("Initializing target: "+str(self.Name))

		## Check all array sizes for consistency
		_all_good = True
		_all_good = _all_good and (self.N_isotopes==len(self.Z_array))
		_all_good = _all_good and (self.N_isotopes==len(self.A_array))
		_all_good = _all_good and (self.N_isotopes==len(self.FF_Rn_array))
		_all_good = _all_good and (self.N_isotopes==len(self.Abundance_array))
		_all_good = _all_good and (self.N_isotopes==len(self.Isotope_list))

		if not _all_good:
			print ("Error: One or more aray has length not equal to N_isotopes")
			return

		## Rescale abundances to sum to unity
		_abnd_sum = numpy.sum(self.Abundance_array)
		self.Abundance_array = self.Abundance_array / _abnd_sum

		## Calculate average A and Z
		self.Z_avg = self.WeightedAverage(self.Z_array)
		self.A_avg = self.WeightedAverage(self.A_array)

		## Populate the isotope list
		for i in numpy.arange(self.N_isotopes):
			_Z     = self.Z_array[i]
			_A     = self.A_array[i]
			_FF_Rn = self.FF_Rn_array[i]
			_Name  = "Xe"+str(int(_A))
			_abnd  = self.Abundance_array[i]

			print("Isotope: "+_Name,
				  "Z="+str(int(_Z)),
				  "A="+str(int(_A)),
				  str(100.*_abnd)+"%")

			self.Isotope_list[i] = SimpleTarget(_A, _Z, _Name, _FF_Rn)
			self.Isotope_list[i].FormFactor = DarkMatterUtilities.FormFactors.FormFactor(self.Isotope_list[i])

	def WeightedAverage(self, _values):
		_n_values  = len(_values)

		_weights   = self.Abundance_array
		_n_weights = len(_weights)

		_same_size    = (_n_values==_n_weights)
		_singular_val = (_n_values==1)

		if not (_same_size or _singular_val):
			print("Error: in WeightedAverage -- un-matched number of values and weights")
			return 0

		_wght_mean = numpy.average(_values, weights=_weights, axis=0)
		return _wght_mean

	def ReducedMass_Nucleus_amu(self, _mass_amu):
		_ans_array = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_ans_array[i] = self.Isotope_list[i].ReducedMass_Nucleus_amu(_mass_amu)
		return self.WeightedAverage(_ans_array)

	def ReducedMass_Nucleus_GeV(self, _mass_GeV):
		_ans_array = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_ans_array[i] = self.Isotope_list[i].ReducedMass_Nucleus_GeV(_mass_GeV)
		return self.WeightedAverage(_ans_array)

	def ReducedMass_Nucleus_kg(self, _mass_kg):
		_ans_array = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_ans_array[i] = self.Isotope_list[i].ReducedMass_Nucleus_kg(_mass_kg)
		return self.WeightedAverage(_ans_array)

	def GetMomentumTransfer(self, _Er_keV):
		_ans_array = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_ans_array[i] = self.Isotope_list[i].GetMomentumTransfer(_Er_keV)
		return self.WeightedAverage(_ans_array)

	def GetQRn(self, _Er_keV):
		_ans_array = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_ans_array[i] = self.Isotope_list[i].GetQRn(_Er_keV)
		return self.WeightedAverage(_ans_array)

	def LindhardFactor(self, _Er_keV):
		_ans_array = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_ans_array[i] = self.Isotope_list[i].LindhardFactor(_Er_keV)
		return self.WeightedAverage(_ans_array)

	def RecoilEnergyAngularDist(self, _theta):
		_ans_array = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_ans_array[i] = self.Isotope_list[i].RecoilEnergyAngularDist(_theta)
		return self.WeightedAverage(_ans_array)

	def RecoilEnergyMax_AnyParticle_keV(self, _incoming_mass_GeV, _incoming_E_keV):
		_ans_array = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_ans_array[i] = self.Isotope_list[i].RecoilEnergyMax_AnyParticle_keV(_incoming_mass_GeV, _incoming_E_keV)
		return self.WeightedAverage(_ans_array)

	def RecoilEnergyMax_DM_keV(self, _dm_mass_GeV):
		_ans_array = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_ans_array[i] = self.Isotope_list[i].RecoilEnergyMax_DM_keV(_dm_mass_GeV)
		return self.WeightedAverage(_ans_array)

	def MinimumDetectableMass_GeV(self, _threshold_keV):
		_ans_array = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_ans_array[i] = self.Isotope_list[i].MinimumDetectableMass_GeV(_threshold_keV)
		return self.WeightedAverage(_ans_array)

## This class defines a complex target consisting of 
## xenon isotopes in concentrations determined by 
## their natural abundance. Most results from methods
## are calculated as a weighted average over the 
## specified isotopic abundance.
class NaturalXenonTarget_2:
	## Long name of material
	Name 			= "Xenon (natural abundance)"
	
	## Number of unique isotopes in the target
	N_isotopes      = 13

	## Atomic number (# of protons) [amu]
	Z			= 54.0*numpy.ones(N_isotopes)

	## Atomic mass (# of protons + neutrons) [amu]
	# An array with an entry for each of 
	# the atomic masses that make up the target
	# Typically these are the isotopes present in 
	# natural abundance.
	A        = numpy.array([124.0 , 125.0 , 126.0 , 127.0 , 128.0 , 129.00 , 130.0 , 131.00 , 132.00 , 133.0 , 134.00 , 135.0 , 136.0 ])

	## Natural abundance fraction [fraction]
	# An array with an entry for each of 
	# the fractions (specified as a percent) 
	# present in target, correlated with Z
	# The entire array is divided by 100
	Abundace_array = numpy.array([0.095 , 0.000 , 0.089 , 0.000 , 1.910 , 26.401 , 4.071 , 21.232 , 26.909 , 0.000 , 10.436 , 0.000 , 8.857 ]) / 100.
	
	## Convert the atomic mass array into other useful units
	NuclearMass_GeV = A * amu_to_GeV	# nuclear mass in GeV
	NuclearMass_kg  = A * amu_to_kg		# nuclear mass in kg

	## Numbers for the form factor calculation
	FF_Rn		    	= 4.7808    				# nuclear form factor radius [fm]

	## Secondary numbers
	A_avg               = 0.0
	NuclearMass_GeV_avg = 0.0
	NuclearMass_kg_avg  = 0.0

	def __init__(self):
		print ("Initializing target: "+str(self.Name))

		self.A_avg               = self.WeightedAverage(self.A)
		self.NuclearMass_GeV_avg = self.WeightedAverage(self.NuclearMass_GeV)
		self.NuclearMass_kg_avg  = self.WeightedAverage(self.NuclearMass_kg)

		self.N_isotopes          = len(self.Abundace_array)

	def WeightedAverage(self, values):
		n_values  = len(values)

		weights   = self.Abundace_array
		n_weights = len(weights)

		same_size    = (n_values==n_weights)
		singular_val = (n_values==1)

		if not (same_size or singular_val):
			print("Error: in WeightedAverage -- un-matched number of values and weights")
			return 0

		wght_mean = numpy.average(values, weights=weights, axis=0)
		return wght_mean

	def ReducedMass_AvgNucleus_amu(self, _mass_amu):
		# Returns the reduced mass of the system consisting of target nucleus and a second mass specified in amu
		_numerator   = self.A_avg * _mass_amu
		_denominator = self.A_avg + _mass_amu
		return (_numerator / _denominator) 

	def ReducedMass_AvgNucleus_GeV(self, _mass_GeV):
		# Returns the reduced mass of the system consisting of target nucleus and a second mass specified in GeV
		_numerator   = self.NuclearMass_GeV_avg * _mass_GeV
		_denominator = self.NuclearMass_GeV_avg + _mass_GeV
		return (_numerator / _denominator) 

	def ReducedMass_AvgNucleus_kg(self, _mass_kg):
		# Returns the reduced mass of the system consisting of target nucleus and a second mass specified in kg
		_numerator   = self.NuclearMass_kg_avg * _mass_kg
		_denominator = self.NuclearMass_kg_avg + _mass_kg
		return (_numerator / _denominator) 

	def ReducedMass_Nucleus_amu(self, _mass_amu):
		# Returns the reduced mass of the system consisting of target nucleus and a second mass specified in amu
		_numerator   = self.A * _mass_amu
		_denominator = self.A + _mass_amu
		return (_numerator / _denominator) 

	def ReducedMass_Nucleus_GeV(self, _mass_GeV):
		# Returns the reduced mass of the system consisting of target nucleus and a second mass specified in GeV
		_numerator   = self.NuclearMass_GeV * _mass_GeV
		_denominator = self.NuclearMass_GeV + _mass_GeV
		return (_numerator / _denominator) 

	def ReducedMass_Nucleus_kg(self, _mass_kg):
		# Returns the reduced mass of the system consisting of target nucleus and a second mass specified in kg
		_numerator   = self.NuclearMass_kg * _mass_kg
		_denominator = self.NuclearMass_kg + _mass_kg
		return (_numerator / _denominator) 

	def GetMomentumTransfer(self, _Er_keV):
		_q_vals = numpy.zeros(shape=(self.N_isotopes,len(self.NuclearMass_GeV)))
		for i in numpy.arange(self.N_isotopes):
			_q2 = 2.0 * self.NuclearMass_GeV * _Er_keV
			_q_vals[i]  = numpy.sqrt(_q2)
		_q_wght_avg = self.WeightedAverage(_q_vals)
		return _q_wght_avg

	def GetQRn(self, _Er_keV):
		_q   = self.GetMomentumTransfer(_Er_keV)
		_qrn = _q * ( self.FF_Rn / hbarc_MeV_fm) 
		return _qrn

	def LindhardFactor(self, _Er_keV):
		# Determine the lindhard factor for a nuclear recoil of a specified recoil energy in keV
		_LF_vals = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_Z 	= self.Z[i]
			_e 	= 11.5 * _Er_keV * numpy.power(_Z, -7./3.)
			_k 	= 0.133 * numpy.power(_Z, 2./3.) * numpy.power(self.A[i], -1./2.)
			_g 	= (3.0 * numpy.power(_e, 0.15)) + (0.7 * numpy.power(_e, 0.6)) + _e
			_LF_vals[i]	= (_k * _g) / ( 1. + (_k * _g))
		_LF_wght_avg = self.WeightedAverage(_LF_vals)
		return _LF_wght_avg

	def RecoilEnergyAngularDist(self, _theta, _incoming_particle_mass_GeV=m_neutron_GeV):
		# Given a fixed angle for the outgoing particle, what is the 
		# energy fraction deposited (this is recoil energy)
		_frac_vals = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_m1 	= _incoming_particle_mass_GeV
			_m2 	= self.NuclearMass_GeV[i]
			_Mfrac 	= _m1 * _m1 / numpy.power(_m1+_m2,2.)
			_term1 	= numpy.cos(_theta)
			_rtarg	= (_m2**2./_m1**2.) - numpy.power(numpy.sin(_theta),2.)
			_term2  = numpy.sqrt(_rtarg)
			_sqfac	= numpy.power(_term1 + _term2 , 2)
			_frac_vals[i] = 1. - _Mfrac*_sqfac
		_frac_avg = self.WeightedAverage(_frac_vals)
		return _frac_avg

	def RecoilEnergyMax_AnyParticle_keV(self, _incoming_mass_GeV, _incoming_E_keV):
		# The maximum energy recoil a particle of a given energy can produce.
		# Calculated from classical 2-body kinematics
		_E_out_vals = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_m1 	= _incoming_mass_GeV
			_m2 	= self.NuclearMass_GeV[i]
			_scale  = 4.0 * _m1 * _m2 / numpy.power(_m1+_m2,2)
			_E_out_vals[i] = _incoming_E_keV * _scale
		_E_out_avg = self.WeightedAverage(_E_out_vals)
		return _E_out_avg

	def RecoilEnergyMax_DM_keV(self, _dm_mass_GeV):
		# Calculates the maximum recoil energy a DM particle of specified mass moving at the Milky Way
		# escape velocity can produce in this detector
		_E_max_dm_keV	= 0.5 * (_dm_mass_GeV * GeV_to_kg) * (MW_esc_vel_ms)**2 * (1./keV_to_J)
		return self.RecoilEnergyMax_AnyParticle_keV(_dm_mass_GeV, _E_max_dm_keV)

	def MinimumDetectableMass_GeV(self, _threshold_keV):
		# For a given recoil energy threshold, what is the smallest mass particle that can produce a recoil 
		# of that energy.
		# Calculated from classical 2-body kinematics
		_M_out_vals = numpy.zeros(self.N_isotopes)
		for i in numpy.arange(self.N_isotopes):
			_threshold_GeV 	= _threshold_keV * 1e-6
			_M_GeV			= self.NuclearMass_GeV[i]
			_beta 			= MW_esc_vel_ms / c_ms
			_M_out_vals[i]  = numpy.sqrt(1./2. * _threshold_GeV * _M_GeV / _beta**2)
		_M_out_avg = self.WeightedAverage(_M_out_vals)
		return _M_out_avg
