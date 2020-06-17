import numpy
from scipy.integrate import quad
from scipy.special import erf, jv

from DarkMatterUtilities.Constants import *

## This class defines a the veloicity distribution
## of the DM and its integral. There are 2 models implemented
## (0) a Maxwell-Boltzmann model, and
## (1) the standard halo model
class Halo:
	Model     = 0
	ModelName = "ModelName"
	fVc_km_s  = 220

	def __init__(self, _HaloModel):
		self.Model = _HaloModel

	def GetHaloPDF_ms(self, _v_ms):
		if (self.Model == 0):
			self.ModelName = "Maxwell-Boltzmann"
			return self.MaxwellBoltzmann_PDF_ms(_v_ms)
		# elif (self.Model == 1):
		# 	return self.StandardHaloModel_PDF_ms(_v_ms)
		# elif (self.Model == 2):
		# 	return 1.0
		else:
			self.ModelName = "Standard Halo"
			return self.StandardHaloModel_PDF_ms(_v_ms)

	def GetHaloIntegral_ms(self, _vmin_ms):
		if (self.Model == 0):
			self.ModelName = "Maxwell-Boltzmann"
			return self.MaxwellBoltzmann_Integral_ms(_vmin_ms)
		# elif (self.Model == 1):
		# 	return self.StandardHaloModel_Integral_ms(_vmin_ms)
		# elif (self.Model == 2):
		# 	return self.WoodwardHalo_Integral_ms(_vmin_ms)
		else:
			self.ModelName = "Standard Halo"
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

	def StandardHaloModel_Integral_DW_ms(self, _vmin_ms):
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

## This class defines a general dark matter
## particle with a local energy density specifed
## in (GeV/c^2/cm^3), a mass (specified in 
## GeV/c^2) and a per-nucleon cross section
## (specified in cm^2)
class DarkMatter:
	Rho0 	    = 0.3				# GeV / cm^3
	Mass_GeV 	= 100.0				# GeV
	Sigma_cm2 	= 1e-46				# cm^2
	Rmass_DM_proton_kg  = 1.0		# kg
	Rmass_DM_proton_GeV = 1.0		# GeV
	HaloModel   = Halo(1)

	def __init__(self, _Mass, _Sigma, _HaloModelType=1):
		self.Mass_GeV 	    	 = _Mass 		# Mass of DM particle in GeV/c^2
		self.Sigma_cm2 		     = _Sigma 		# Spin-independent DM-proton elastic scattering cross-section in cm^2
		self.Rmass_DM_proton_GeV = (( _Mass * m_proton_GeV ) / ( _Mass + m_proton_GeV))
		self.Rmass_DM_proton_kg	 = self.Rmass_DM_proton_GeV * GeV_to_kg 
		self.HaloModel           = Halo(_HaloModelType)