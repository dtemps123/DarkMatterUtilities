import numpy
from scipy.integrate import quad
from scipy.special import erf, jv

from DarkMatterUtilities.Constants import *

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
	FF_type			= 4					# Which form factor to use
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
		_q 		= numpy.sqrt(2.0 * self.NuclearMass_GeV * _Er_keV)								# MeV / c
		_q_fm   = numpy.sqrt(2.0 * self.NuclearMass_GeV * _Er_keV * 1e-6) / 0.197				# fm^-1 
		_qrn 	= _q * ( _rn / hbarc_MeV_fm) 											 	# dimensionless
		_qs		= _q * ( _s  / hbarc_MeV_fm)												# dimensionless

		if   ( self.FF_type == 0 ):
			## Lewin & Smith -- thin shell: exp[-(q r_n)^(2/3) / 3]
			return numpy.exp(-1.0 * numpy.power(_qrn,2./3.) / 3.0 )
		
		elif ( self.FF_type == 1 ):
			## Lewin & Smith -- thin shell: [ sin(q r_n) / (q r_n) ]^2
			## Confirmed to match spectrum from L&S
			return numpy.power( numpy.sin(_qrn) / _qrn , 2.)
		
		elif ( self.FF_type == 2 ):
			## Lewin & Smith -- solid sphere: exp[-(q r_n)^(2/3) / 5]
			return numpy.exp(-1.0 * numpy.power(_qrn,2./3.) / 5.0 )

		elif ( self.FF_type == 3 ):
			## Lewin & Smith -- solid sphere: { 3 [ sin(q r_n) - q r_n cos(q r_n)] / (q r_n)^3 }^2
			## Confirmed to match spectrum from L&S
			_arg1 = numpy.sin(_qrn) - (_qrn * numpy.cos(_qrn)) 
			_arg2 = _arg1 / numpy.power(_qrn,3.)
			return numpy.power(3.0 * _arg2,2.)
			# return 3.0*_arg2 * numpy.exp( - (_qs**2)/2.)
		
		elif (self.FF_type == 4 ):
			return self.HelmFormFactor(_Er_keV)
		
		elif (self.FF_type == 5 ):
			return self.HelmFormFactor_DW(_Er_keV)
		
		else:
			if _Er_keV <= 1e5:
				return 1.0
			else:
				return 0

	def HelmFormFactor(self, _Er_keV):
		# [arXiv:0608035] Duda et al 2007 (consistent with DMCalc implementation)
		_a 	 = 0.52														# fm
		_s 	 = 0.9														# fm
		_c	 = (1.23 * numpy.power(self.A,1./3.)) - 0.60					# fm
		_R1	 = numpy.sqrt(  (_c**2) 
						 + ( (7./3.)*(numpy.pi**2)*(_a**2) ) 
						 - ( 5.0*(_s**2) )   )							# fm

		_q 		= numpy.sqrt(2.0 * _Er_keV * self.NuclearMass_GeV)			# MeV / c
		_qs		= (_s  / hbarc_MeV_fm) * _q								# dimensionless
		_qR1 	= (_R1 / hbarc_MeV_fm) * _q								# dimensionless

		_expfac	= numpy.exp(-1.0 * _qs**2)									# dimensionless
		_jterm	= ( 3.0 * jv(1, _qR1) / _qR1 )**2						# dimensionless

		return _jterm * _expfac

	def HelmFormFactor_DW(self, _Er_keV):
		_s 	 = 1.0														# fm
		_R	 = (1.20 * numpy.power(self.A,1./3.))					# fm
		_r   = numpy.power((numpy.power(_R,2)-5*numpy.power(_s,2)),0.5) # m
		
		_q 		= numpy.sqrt(2.0 * _Er_keV * self.NuclearMass_GeV)			# MeV / c
		_qs		= (_s / hbarc_MeV_fm) * _q								# dimensionless
		_qr  	= (_r / hbarc_MeV_fm) * _q								# dimensionless

		_expfac	= numpy.exp(-numpy.power(_qs,2))							# dimensionless
		_jterm	= (numpy.sin(_qr)/(numpy.power(_qr,2)))-(numpy.cos(_qr)/_qr)						# dimensionless

		return pow(((3*_jterm)/(_qr)),2)*_expfac

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
		_M_GeV			= _target.NuclearMass_GeV
		_beta 			= MW_esc_vel_ms / c_ms
		return numpy.sqrt(1./2. * _threshold_GeV * _M_GeV / _beta**2)

class Halo:
	Model = 0
	fVc_km_s = 220

	def __init__(self, _HaloModel):
		self.Model = _HaloModel

	def GetHaloPDF_ms(self, _v_ms):
		if (self.Model == 0):
			return self.MaxwellBoltzmann_PDF_ms(_v_ms)
		elif (self.Model == 1):
			return self.StandardHaloModel_PDF_ms(_v_ms)
		elif (self.Model == 2):
			return 1.0
		else:
			return self.StandardHaloModel_PDF_ms(_v_ms)

	def GetHaloIntegral_ms(self, _vmin_ms):
		if (self.Model == 0):
			return self.MaxwellBoltzmann_Integral_ms(_vmin_ms)
		elif (self.Model == 1):
			return self.StandardHaloModel_Integral_ms(_vmin_ms)
		elif (self.Model == 2):
			return self.WoodwardHalo_Integral_ms(_vmin_ms)
		else:
			return self.StandardHaloModel_Integral_ms(_vmin_ms)

	def MaxwellBoltzmann_PDF_ms(self, _v_ms):
		# Returns the value of the PDF (units of (m/s)^-1) given a velocity in m/s
		# c.f. Anton N Baushev 2012 J. Phys.: Conf. Ser. 375 012048 [eq 1]
		# Note this is independent of all dark matter particle properties
		_V 			= MW_esc_vel_ms
		_v0			= Solar_vel_ms
		_norm 		= 1./ ( erf(_V/_v0) - (2./numpy.sqrt(numpy.pi))*(_V/_v0)* numpy.exp( - _V**2 / _v0**2 ) )
		_exp_fac	= numpy.exp( - _v_ms**2 / _v0**2 )
		_fac1		= 4. / numpy.sqrt(numpy.pi)
		_fac2		= _v_ms**2 / _v0**3
		MB_pdf      = _norm * _fac1 * _fac2 * _exp_fac
		return MB_pdf / _v_ms

	def MaxwellBoltzmann_Integral_ms(self, _vmin_ms):
		if (_vmin_ms > MW_esc_vel_ms):
			return 0
		
		_MBfactor   = quad(self.MaxwellBoltzmann_PDF_ms, _vmin_ms , MW_esc_vel_ms)
		return _MBfactor[0]	# m^-1  x  s 

	def StandardHaloModel_PDF_ms(self, _v_ms):
		velocity  = _v_ms / 1e3  ## get this in km/s
		v_esc_kms = MW_esc_vel_ms / 1e3
		v_sol_kms = Solar_vel_ms  / 1e3

		_k0 = numpy.power(numpy.pi * self.fVc_km_s*self.fVc_km_s, 1.5);
		_k1 = ( _k0 * erf(v_esc_kms/self.fVc_km_s)  
		       - (2./numpy.sqrt(numpy.pi)) * numpy.sqrt(v_esc_kms/self.fVc_km_s) * numpy.exp(-numpy.power(v_esc_kms/self.fVc_km_s,2)) )
		flab = velocity * (numpy.exp(-numpy.power((velocity - v_sol_kms)/self.fVc_km_s,2)) -
		                   numpy.exp(-numpy.power((velocity + v_sol_kms)/self.fVc_km_s,2)))
		return flab / (_k1) * (numpy.pi *self.fVc_km_s*self.fVc_km_s/v_sol_kms)

	def StandardHaloModel_Integral_ms(self, _vmin_ms):
		## McCabe [arXiv:1005.0579]
		x_esc	= (MW_esc_vel_ms*1e-3) / self.fVc_km_s
		x_min	= (_vmin_ms*1e-3)	   / self.fVc_km_s
		x_earth = (Solar_vel_ms*1e-3)  / self.fVc_km_s

		zeta = 0 ## this will be in (km / s)^-1
		beta = 1 ## this controls the escape velocity cut: beta=0 (hard cut), beta=1 (soft cut)	

		norm = (  numpy.power(numpy.pi,1.5)
				* numpy.power(self.fVc_km_s,3.) 
				* (  erf(x_esc) - 4./numpy.sqrt(numpy.pi)*numpy.exp(-x_esc*x_esc) 
				   * (0.5*x_esc + beta*numpy.power(x_esc,3.)/3.)
				  )
			   )

		if  ( x_earth+x_min < x_esc ):
			zeta = ( numpy.power(numpy.pi,1.5)*numpy.power(self.fVc_km_s,2.)/(2.*norm*x_earth)
					 * ( erf(x_min+x_earth) - erf(x_min-x_earth)
					 - 4.*x_earth/numpy.sqrt(numpy.pi)*numpy.exp(-x_esc*x_esc)
					 * (1 + beta*(x_esc*x_esc - x_earth*x_earth/3. - x_min*x_min)) )
				   )

		elif( x_min >= numpy.abs(x_esc-x_earth) and x_min < (x_earth+x_esc) ):
			zeta = ( numpy.power(numpy.pi,1.5)*numpy.power(self.fVc_km_s,2.)/(2.*norm*x_earth)
					 * ( erf(x_esc) + erf(x_earth-x_min)
					 - 2./numpy.sqrt(numpy.pi)*numpy.exp(-x_esc*x_esc)
					 *(x_esc + x_earth - x_min - beta/3.*(x_earth-2.*x_esc-x_min)
					 *numpy.power(x_esc+x_earth-x_min,2.)) )
				   )

		elif( x_earth>x_min+x_esc ):
			zeta = 1./(self.fVc_km_s*x_earth)

		elif( (x_earth + x_esc) < x_min ):
			zeta = 0  
		 	
		zeta_ms = zeta / 1e3
		return zeta_ms

	def WoodwardHalo_Integral_ms(self, _vmin_ms):
		# astrophysical constants
		_vmin_kms = _vmin_ms / 1.e3
		v_esc   =MW_esc_vel_ms/1.e3 # km /s
		v_0     =self.fVc_km_s # km /s
		v_earth =230. # km /s
		y       =v_earth/v_0
		z       =v_esc/v_0
		contrib =0
		x       = _vmin_kms/v_0   
		Nesc_z = (erf(z)-((2*z*numpy.exp(-numpy.power(z,2)))/numpy.power(numpy.pi,0.5)))
		if x<abs(y-z):
			contrib=(1./(2*Nesc_z*v_0*y))*((erf(x+y))-(erf(x-y))-((4./(numpy.power(numpy.pi,0.5)))*y*numpy.exp(-numpy.power(z,2))))
		if x>abs(y-z) and x<(y+z):
			contrib=(1./(2*Nesc_z*v_0*y))*((erf(z))-(erf(x-y))-((2./(numpy.power(numpy.pi,0.5)))*(y+z-x)*numpy.exp(-numpy.power(z,2))))    
		if x>(y+z):
			contrib=0

		return contrib / 1.e3 ## to get m/s

# Defines a WIMP-like dark matter particle
# Default is a 100 GeV/c^2 WIMP with local density of 0.3 GeV / cm^3
class DarkMatter:
	Rho0 	= 0.3				# GeV / cm^3
	Mass 	= 100.0				# GeV
	Sigma 	= 1e-46				# cm^2
	Rmass_DM_proton_kg  = 1.0		# kg
	Rmass_DM_proton_GeV = 1.0		# kg
	HaloModel = Halo(1)

	def __init__(self, _Mass, _Sigma):
		self.Mass 				 = _Mass 		# Mass of DM particle in GeV/c^2
		self.Sigma 				 = _Sigma 		# Spin-independent DM-proton elastic scattering cross-section in cm^2
		self.Rmass_DM_proton_GeV = (( _Mass * m_proton_GeV ) / ( _Mass + m_proton_GeV))
		self.Rmass_DM_proton_kg	 = self.Rmass_DM_proton_GeV * GeV_to_kg 