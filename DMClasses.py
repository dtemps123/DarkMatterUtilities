import numpy as n
import matplotlib.pyplot as pyp
from scipy.special import erf, jv

## -- Unit Conversions -- ##
keV_to_J		= 1.60218e-16
kg_to_kev		= 5.60958865e32
GeV_to_kg		= 1.78266191e-27
amu_to_GeV		= 0.93149431
amu_to_kg		= 1.6605402e-27

## -- Physical Constants -- ##
c_ms			= 3e8
hbarc_MeV_fm	= 197.323
m_proton_GeV	= 0.93827
m_proton_amu	= m_proton_GeV * (1./amu_to_GeV)	#1.007276
m_neutron_GeV	= 0.93957
m_neutron_amu	= m_neutron_GeV * (1./amu_to_GeV)	#1.008664
MW_esc_vel_ms	= 544e3
Solar_vel_ms	= 230e3
e_charge        = 1.602176634e-19

## -- Class Definitions -- ##

# Defines a detector made of a single material.
# Default is a 100 kg xenon target
class Target:
	Name 			= "Xenon"			# name of target
	A 				= 1.0				# amu "dimensionless"
	Z 				= 1.0				# amu "dimensionless"
	TotalMass 		= 100.0				# kg of detector
	ExposureTime 	= 1.0				# years of operation
	NuclearMass_GeV = A * amu_to_GeV	# nuclear mass in GeV
	NuclearMass_kg	= A * amu_to_kg		# nuclear mass in kg
	FF_type			= 1					# Which form factor to use
	FF_Rn			= 1.0				# nuclear form factor radius [fm]
	FF_alpha		= 1./3.				# nuclear form factor parameterization [dimensionless] only impacts FF type 2

	def __init__(self, _A, _Z, _TotalMass, _ExposureTime, _Name, _FF_Rn):
		self.A = _A
		self.Z = _Z
		self.TotalMass 			= _TotalMass
		self.ExposureTime 		= _ExposureTime
		self.NuclearMass_GeV 	= _A * amu_to_GeV
		self.NuclearMass_kg		= _A * amu_to_kg
		self.Name 				= _Name
		self.FF_Rn 				= _FF_Rn

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

	def FormFactor(self, _Er_keV):
		# Evaluates the dimensionless form factor for a given momentum transfer (recoil energy)
		# See Lewin and Smith (1996) for description
		# Helm form factor:
		_alpha 	= self.FF_alpha 															# dimensionless
		_rn 	= self.FF_Rn 																# fm
		_s 		= 1.0																		# fm
		_q 		= n.sqrt(2.0 * self.NuclearMass_GeV * _Er_keV)								# MeV / c
		_q_fm   = n.sqrt(2.0 * self.NuclearMass_GeV * _Er_keV * 1e-6) / 0.197				# fm^-1 
		_qrn 	= _q * ( _rn / hbarc_MeV_fm) 											 	# dimensionless
		_qs		= _q * ( _s  / hbarc_MeV_fm)												# dimensionless

		if   ( self.FF_type == 0 ):
			## Lewin & Smith -- thin shell: exp[-(q r_n)^(2/3) / 3]
			return n.exp(-1.0 * n.power(_qrn,2./3.) / 3.0 )
		
		elif ( self.FF_type == 1 ):
			## Lewin & Smith -- thin shell: [ sin(q r_n) / (q r_n) ]^2
			## Confirmed to match spectrum from L&S
			return n.power( n.sin(_qrn) / _qrn , 2.)
		
		elif ( self.FF_type == 2 ):
			## Lewin & Smith -- solid sphere: exp[-(q r_n)^(2/3) / 5]
			return n.exp(-1.0 * n.power(_qrn,2./3.) / 5.0 )

		elif ( self.FF_type == 3 ):
			## Lewin & Smith -- solid sphere: { 3 [ sin(q r_n) - q r_n cos(q r_n)] / (q r_n)^3 }^2
			## Confirmed to match spectrum from L&S
			_arg1 = n.sin(_qrn) - (_qrn * n.cos(_qrn)) 
			_arg2 = _arg1 / n.power(_qrn,3.)
			return n.power(3.0 * _arg2,2.)
			# return 3.0*_arg2 * n.exp( - (_qs**2)/2.)

		elif ( self.FF_type == 4):
			## Duda & Kemper -- Fermi Two parameter distribution (values in appendix using Xe131)
			_a_val = 0.523	# fm
			_c_val = 5.6384	# fm
			_rho_c = self.Z * e_charge / ( _c_val * n.log(n.exp(_a_val/_c_val)+1) )
			_r_val = 1.0 / _q_fm
			return _rho_c / ( n.exp((_r_val-_c_val)/_a_val) + 1 )


		elif ( self.FF_type == 5 ):
			return n.exp( -_alpha * (_qrn**2) )
		elif ( self.FF_type == 6 ):
			return 3.0*( n.sin(_qrn) - _qrn * n.cos(_qrn)) / (_qrn**3) * n.exp( - (_qs**2)/2.)
		else:
			if _Er_keV <= 1e5:
				return 1.0
			else:
				return 0

	def HelmFormFactor(self, _Er_keV):
		# [arXiv:1412.6091] Vietze et all 2014
		_J		= 3./2.
		_c0		= 3.0e-2
		_c 		= 1.23 * (self.A**(1./3.)) - 0.6						# fm
		_a 		= 0.52													# fm
		_s 		= 1.0													# fm
		_rn		= n.sqrt(_c**2 + (7./3.)*n.pi**2*_a**2 - 5.0*_s**2)		# fm

		_Ss0	= (self.A**2) * (_c0**2) * (2.0*_J + 1.0) / (4.0*n.pi)	# dimensionless

		_q 		= n.sqrt(2.0 * _Er_keV * self.NuclearMass_GeV)			# MeV / c
		_qs		= (_s  / hbarc_MeV_fm) * _q								# dimensionless
		_qrn 	= (_rn / hbarc_MeV_fm) * _q								# dimensionless

		_expfac	= n.exp(-1.0 * _qs)										# dimensionless
		_jterm	= ( 3.0 * jv(1, _qrn) / _qrn )**2						# dimensionless

		# _expfac	= n.exp(-1.0 * _qs**2 / 2.0)						# dimensionless
		# _jterm	= ( 3.0 * jv(1, _qrn) / _qrn )						# dimensionless

		return _Ss0 * _jterm * _expfac

	def HelmFormFactor_v2(self, _Er_keV):
		# [arXiv:0608035] Duda et al 2007
		_a 		= 0.52													# fm
		_s 		= 0.9													# fm
		_c      = (1.23 * n.power(self.A,1./3.)) - 0.60					# fm
		_R1     = n.sqrt(  (_c**2) 
			             + ( (7./3.)*(n.pi**2)*(_a**2) ) 
			             - ( 5.0*(_s**2) )   )							# fm

		_q 		= n.sqrt(2.0 * _Er_keV * self.NuclearMass_GeV)			# MeV / c
		_qs		= (_s  / hbarc_MeV_fm) * _q								# dimensionless
		_qR1 	= (_R1 / hbarc_MeV_fm) * _q								# dimensionless

		_expfac	= n.exp(-1.0 * _qs**2)									# dimensionless
		_jterm	= ( 3.0 * jv(1, _qR1) / _qR1 )**2						# dimensionless

		return _jterm * _expfac

	def LindhardFactor(self, _Er_keV):
		# Determine the lindhard factor for a nuclear recoil of a specified recoil energy in keV
		_Z 	= self.Z
		_e 	= 11.5 * _Er_keV * n.power(_Z, -7./3.)
		_k 	= 0.133 * n.power(_Z, 2./3.) * n.power(self.A, -1./2.)
		_g 	= (3.0 * n.power(_e, 0.15)) + (0.7 * n.power(_e, 0.6)) + _e
		_LF	= (_k * _g) / ( 1. + (_k * _g))
		return _LF

	def RecoilEnergyAngularDist(self, _theta):
		# Given a fixed angle for the outgoing particle, what is the 
		# energy fraction deposited (this is recoil energy)
		_m1 	= m_neutron_GeV
		_m2 	= self.NuclearMass_GeV
		_Mfrac 	= _m1 * _m1 / n.power(_m1+_m2,2.)
		_term1 	= n.cos(_theta)
		_rtarg	= (_m2**2./_m1**2.) - n.power(n.sin(_theta),2.)
		_term2  = n.sqrt(_rtarg)
		_sqfac	= n.power(_term1 + _term2 , 2)
		return 1. - _Mfrac*_sqfac

	def RecoilEnergyMax_AnyParticle_keV(self, _incoming_mass_GeV, _incoming_E_keV):
		# The maximum energy recoil a particle of a given energy can produce.
		# Calculated from classical 2-body kinematics
		_m1 	= _incoming_mass_GeV
		_m2 	= self.NuclearMass_GeV
		_scale  = 4.0 * _m1 * _m2 / n.power(_m1+_m2,2)

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
		_M_GeV			= _target.NuclearMass_GeV
		_beta 			= MW_esc_vel_ms / c_ms
		return n.sqrt(1./2. * _threshold_GeV * _M_GeV / _beta**2)

# Defines a WIMP-like dark matter particle
# Default is a 100 GeV/c^2 WIMP with local density of 0.3 GeV / cm^3
class DarkMatter:
	Rho0 	= 0.3				# GeV / cm^3
	Mass 	= 100.0				# GeV
	Sigma 	= 1e-46				# cm^2
	Rmass_DM_proton = 1.0		# kg

	def __init__(self, _Mass, _Sigma):
		self.Mass 				= _Mass 		# Mass of DM particle in GeV/c^2
		self.Sigma 				= _Sigma 		# Spin-independent DM-proton elastic scattering cross-section in cm^2
		self.Rmass_DM_proton	= (( _Mass * m_proton_GeV ) / ( _Mass + m_proton_GeV)) * GeV_to_kg

	def MaxwellBoltzmann_PDF_ms(self, _v):
		# Returns the value of the PDF (units of (m/s)^-1) given a velocity in m/s
		# c.f. Anton N Baushev 2012 J. Phys.: Conf. Ser. 375 012048 [eq 1]
		# Note this is independent of all dark matter particle properties
		_V 			= MW_esc_vel_ms
		_v0			= Solar_vel_ms
		_norm 		= 1./ ( erf(_V/_v0) - (2./n.sqrt(n.pi))*(_V/_v0)* n.exp( - _V**2 / _v0**2 ) )
		_exp_fac	= n.exp( - _v**2 / _v0**2 )
		_fac1		= 4. / n.sqrt(n.pi)
		_fac2		= _v**2 / _v0**3
		return _norm * _fac1 * _fac2 * _exp_fac

	def Velocity_Dist_ms(self, _v):
		# Turns the MB PDF into a distribution over which to integrate
		dist = self.MaxwellBoltzmann_PDF_ms(_v) / _v
		return dist

	def Plot_VelocityDist(self):
		_fig = pyp.figure()
		_vel_array = n.linspace( 1e-3	, 600	, num=600)
		_MB_dist = self.MaxwellBoltzmann_PDF_ms(_vel_array*1e3)
		_DM_vel_dist = _MB_dist / _vel_array

		_vel_array_esc = n.append( _vel_array[_vel_array <= MW_esc_vel_ms/1e3] , (MW_esc_vel_ms/1e3)+1. )
		_DM_vel_dist_esc = n.append( _DM_vel_dist[_vel_array <= MW_esc_vel_ms/1e3] , 1e-10 )

		pyp.semilogy(_vel_array, _DM_vel_dist, 'b:', label=r'MB $v_\mathrm{esc}=600$ km/s')
		pyp.semilogy(_vel_array_esc, _DM_vel_dist_esc, 'b-', label=r'MB $v_\mathrm{esc}=544$ km/s')

		pyp.xlim([0.,600.])
		pyp.ylim([1e-10,1e-7])

		pyp.xlabel("DM velocity [km/s]")
		pyp.ylabel("Probability Density [a.u.]")

		return _fig

