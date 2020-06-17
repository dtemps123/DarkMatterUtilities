import numpy
from scipy.special import erf, jv

from DarkMatterUtilities.Constants import *

class FormFactor:
	FF_type  = 4
	FF_Rn    = 1.0

	## Get the nuclear mass of the target
	## This will be either a single number
	## for a simple target or an array for
	## a complex target
	NuclearMass_GeV = 1.
	A               = 1.

	## Name the form factor
	FF_Name  = "formfactor"

	def __init__(self, _target, _FF_type=1):
		self.FF_type         = _FF_type
		self.FF_Rn           = _target.FF_Rn
		self.NuclearMass_GeV = _target.NuclearMass_GeV
		self.A               = _target.A

		if   ( self.FF_type == 0 ):
			self.FF_Name  = "Unity"
		
		elif ( self.FF_type == 1 ):
			self.FF_Name  = "Helm"
		
		elif ( self.FF_type == 2 ):
			self.FF_Name  = "Shell1"

		elif ( self.FF_type == 3 ):
			self.FF_Name  = "Shell2"
		
		elif (self.FF_type == 4 ):
			self.FF_Name  = "Sphere1"
		
		elif (self.FF_type == 5 ):
			self.FF_Name  = "Sphere2"	

		elif (self.FF_type == 6 ):
			self.FF_Name  = "TopHat"	

		elif (self.FF_type == 7 ):
			self.FF_Name  = "Helm(DW)"	

		elif (self.FF_type == 8 ):
			self.FF_Name  = "Helm(Original)"	

		else:
			self.FF_Name  = "Unity"

	def EvaluateFormFactorSquared(self, _Er_keV, _renormalize=False):
		_scale = 1
		if   ( self.FF_type == 0 ):
			return self.UnityFormFactor(_Er_keV)
		
		elif ( self.FF_type == 1 ):
			if (_renormalize):
				_scale = self.HelmFormFactor(1e-10)
			return self.HelmFormFactor(_Er_keV) /_scale
		
		elif ( self.FF_type == 2 ):
			if (_renormalize):
				_scale = self.LS_ThinShell_1(1e-10)
			return self.LS_ThinShell_1(_Er_keV) / _scale

		elif ( self.FF_type == 3 ):
			if (_renormalize):
				_scale = self.LS_ThinShell_2(1e-10)
			return self.LS_ThinShell_2(_Er_keV) / _scale
		
		elif (self.FF_type == 4 ):
			if (_renormalize):
				_scale = self.LS_SolidSphere_1(1e-10)
			return  self.LS_SolidSphere_1(_Er_keV) / _scale
		
		elif (self.FF_type == 5 ):
			if (_renormalize):
				_scale = self.LS_SolidSphere_2(1e-10)
			return self.LS_SolidSphere_2(_Er_keV) / _scale

		elif (self.FF_type == 6 ):
			if (_renormalize):
				_scale = self.TopHat(1e-10)
			return self.TopHat(_Er_keV) / _scale

		elif ( self.FF_type == 7 ):
			if (_renormalize):
				_scale = self.HelmFormFactor_DW(1e-10)
			return self.HelmFormFactor_DW(_Er_keV) / _scale

		elif ( self.FF_type == 8 ):
			if (_renormalize):
				_scale = self.HelmFormFactor_Original(1e-10)
			return self.HelmFormFactor_Original(_Er_keV) / _scale

		else:
			return 	self.UnityFormFactor(_Er_keV)
		
	def UnityFormFactor(self, _Er_keV):
		if _Er_keV <= 1e5:
			return 1.0
		else:
			return 0

	def HelmFormFactor(self, _Er_keV):
		# [arXiv:0608035] Duda et al 2007 (consistent with DMCalc implementation)
		_a 	 = 0.52														# fm
		_s 	 = 0.9														# fm
		_c	 = (1.23 * numpy.power(self.A,1./3.)) - 0.60				# fm
		
		_R1	 = numpy.sqrt(  (_c**2) 
						 + ( (7./3.)*(numpy.pi**2)*(_a**2) ) 
						 - ( 5.0*(_s**2) )   )							# fm

		_q   = self.RecoilE_to_Qifm(_Er_keV)
		_qR  = _q * _R1
		_qs  = _q * _s

		# if (_qR < 1e-12):
		# 	return 1.0

		_expfac	= numpy.exp(-1.0 * _qs**2)									# dimensionless
		_jterm	= ( 3.0 * jv(1, _qR) / _qR )						# dimensionless
		_ff_val = _jterm * _expfac / 2.0

		return numpy.power(_ff_val,2.)

	def LS_ThinShell_1(self, _Er_keV):
		## Lewin & Smith -- thin shell: exp[-(q r_n)^(2/3) / 3]
		_rn 	= self.FF_Rn 																# fm
		_q 		= numpy.sqrt(2.0 * self.NuclearMass_GeV * _Er_keV)							# MeV / c
		_qrn 	= _q * ( _rn / hbarc_MeV_fm) 											 	# dimensionless
		
		_pow    = 2./3.
		return numpy.exp(-1.0 * numpy.power(_qrn,_pow) / 3.0 )

	def LS_ThinShell_2(self, _Er_keV):
		## Lewin & Smith -- thin shell: [ sin(q r_n) / (q r_n) ]^2
		## Confirmed to match spectrum from L&S
		_rn 	= self.FF_Rn 																# fm
		_q 		= numpy.sqrt(2.0 * self.NuclearMass_GeV * _Er_keV)							# MeV / c
		_qrn 	= _q * ( _rn / hbarc_MeV_fm) 											 	# dimensionless
		
		return numpy.power( numpy.sin(_qrn) / _qrn , 2.)

	def LS_SolidSphere_1(self, _Er_keV):
		## Lewin & Smith -- solid sphere: exp[-(q r_n)^(2/3) / 5]
		_rn 	= self.FF_Rn 																# fm
		_q 		= numpy.sqrt(2.0 * self.NuclearMass_GeV * _Er_keV)							# MeV / c
		_qrn 	= _q * ( _rn / hbarc_MeV_fm) 											 	# dimensionless
		
		_pow    = 2./3.
		return numpy.exp(-1.0 * numpy.power(_qrn,_pow) / 5.0 )

	def LS_SolidSphere_2(self, _Er_keV):
		## Lewin & Smith -- solid sphere: { 3 [ sin(q r_n) - q r_n cos(q r_n)] / (q r_n)^3 }^2
		## Confirmed to match spectrum from L&S
		_rn 	= self.FF_Rn 																# fm
		_q 		= numpy.sqrt(2.0 * self.NuclearMass_GeV * _Er_keV)							# MeV / c
		_qrn 	= _q * ( _rn / hbarc_MeV_fm) 											 	# dimensionless
		
		_arg1 = numpy.sin(_qrn) - (_qrn * numpy.cos(_qrn)) 
		_arg2 = _arg1 / numpy.power(_qrn,3.)
		return numpy.power(3.0 * _arg2,2.)

	def TopHat(self, _Er_keV):
		_a 	 = 0.52														# fm
		_s 	 = 0.9														# fm
		_c	 = (1.23 * numpy.power(self.A,1./3.)) - 0.60				# fm
		
		_R1	 = numpy.sqrt(  (_c**2) 
						 + ( (7./3.)*(numpy.pi**2)*(_a**2) ) 
						 - ( 5.0*(_s**2) )   )							# fm

		_q   = self.RecoilE_to_Qifm(_Er_keV)
		_qR  = _q * _R1
		_jterm	= ( 3.0 * jv(1, _qR) / _qR )
		return numpy.power(_jterm,2.)

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

	def HelmFormFactor_Original(self, _Er_keV):
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

		return _jterm * _expfac / 2

	def RecoilE_to_Qifm(self, _Er_keV):
		_term1 = 2. * self.NuclearMass_GeV * _Er_keV 
		_term2 = _Er_keV * _Er_keV / 1e6
		_arg   = (_term1 + _term2) / 1e6
		_sqrt  = numpy.sqrt(_arg) / invGeV_to_fm
		return _sqrt	## [fm]^{-1}
		