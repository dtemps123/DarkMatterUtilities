from DMClasses      import *
from DMInteractions import *
import matplotlib.pyplot as pyp

def Plot_DM_Rate_axes(_max_thresh):
	_fig = pyp.figure()
	_thresh_array = n.linspace( 1e-3	, _max_thresh	, num=100)

	pyp.semilogy(_thresh_array, 10.*n.ones(len(_thresh_array)), 'k--')

	pyp.xlim([0.0 , _max_thresh])
	# pyp.ylim([1e-2, 4.0])

	pyp.xlabel("Energy threshold [keV]")
	pyp.ylabel("Rate")

	return _fig, _thresh_array

def Plot_Overlay_Detected_Rate(_fig, _thresh_array, _target, _dm, _c_opt):
	_rate_pts = n.zeros(len(_thresh_array))
	for i in n.arange(len(_rate_pts)):
		_rate_pts[i] = IntegratedRate( _thresh_array[i] , _target, _dm)

	_fig.gca()
	_label = _target.Name
	pyp.semilogy(_thresh_array, _rate_pts, _c_opt, label=_label)

def Plot_DMVelocityDist(_dm):
	_fig = pyp.figure()
	_vel_array_kms = n.linspace( 1e-3	, 600	, num=600)
	_HaloModel_dist = _dm.HaloModel.GetHaloPDF_ms(_vel_array_kms*1e3)

	model_label = ""
	if (_dm.HaloModel.Model == 0):
		model_label = "MB "
		plot_y_lims = n.array([1e-13,1e-10])
	if (_dm.HaloModel.Model == 1):
		model_label = "SHM "
		plot_y_lims = n.array([1e-6,1e-2])

	_vel_array_esc = n.append( _vel_array_kms[_vel_array_kms <= MW_esc_vel_ms/1e3] , (MW_esc_vel_ms/1e3)+1. )
	_DM_vel_dist_esc = n.append( _HaloModel_dist[_vel_array_kms <= MW_esc_vel_ms/1e3] , plot_y_lims[0] )

	pyp.semilogy(_vel_array_kms, _HaloModel_dist , 'b:', label=model_label+r'$v_\mathrm{esc}=600$ km/s')
	pyp.semilogy(_vel_array_esc, _DM_vel_dist_esc, 'b-', label=model_label+r'$v_\mathrm{esc}=544$ km/s')

	pyp.xlim([0.,600.])
	pyp.ylim(plot_y_lims)

	pyp.xlabel("DM velocity [km/s]")
	pyp.ylabel("Probability Density [a.u.]")

	return _fig  

def Plot_Overlay_DM_Vmin(_fig, _threshold_keV, _target, _dm, _c_opt):
	_ax1 = _fig.gca()
	_ylims = _ax1.get_ylim()
	_vmin = MinimumVelocity_ms(_threshold_keV, _target, _dm) / 1e3
	_vmin_Distval = _dm.HaloModel.GetHaloPDF_ms(_vmin*1e3)
	_label = r'$v_\mathrm{min}$ for $M_\chi=$' + str(_dm.Mass) + ' GeV, ' + _target.Name
	pyp.semilogy( [_vmin, _vmin] , [_ylims[0] , _vmin_Distval] , _c_opt, label=_label )

def Plot_MaxRecoilE_axes():
	_fig = pyp.figure()
	_dm_mass_array = n.linspace( 1e-2	, 1e2	, num=1e4)

	pyp.loglog(_dm_mass_array, 1e3*n.ones(len(_dm_mass_array)), 'k--')

	pyp.xlim([1e-2 , 1e2])
	pyp.ylim([1e-4 , 50.])

	pyp.xlabel(r"WIMP mass [GeV/$c^2$]")
	pyp.ylabel("Maximum recoil energy [keV]")

	return _fig, _dm_mass_array

def Plot_Overlay_TargetRecoils( _fig, _target, _dm_mass_array, _c_opt):
	_recoil_pts = n.zeros(len(_dm_mass_array))
	for i in n.arange(len(_recoil_pts)):
		_recoil_pts[i] = _target.RecoilEnergyMax_DM_keV(_dm_mass_array[i] )

	_fig.gca()
	_label = _target.Name
	pyp.loglog(_dm_mass_array , _recoil_pts, _c_opt, label=_label)

def Plot_MinVelocity_axes(_threshold_keV, _dm):
	_fig = pyp.figure()
	_target_mass_array = n.linspace( 1.0 , 150. , num=150)
	_vmin_array = n.zeros(len(_target_mass_array))

	for i in n.arange(len(_target_mass_array)):
		_trgt = Target( _target_mass_array[i] , 1.0 , 1.0 , 1.0 , "tmp" , 1.5 )
		_vmin = MinimumVelocity_ms( _threshold_keV , _trgt , _dm )
		_vmin_array[i] = _vmin / 1e3

	_label = r"$v_\mathrm{min}$"
	# _title = "Min. DM velocity to produce a recoil of " + str(_threshold_keV) + r" keV, $m_\chi=$" + str(_dm.Mass) + r" GeV/$c^2$"
	_title = r"$E_\mathrm{R}=$" + str(_threshold_keV) + r" keV, $m_\chi=$" + str(_dm.Mass) + r" GeV/$c^2$"

	pyp.plot(_target_mass_array, _vmin_array, 'k-', label=_label)
	pyp.plot(_target_mass_array, MW_esc_vel_ms*n.ones(len(_target_mass_array))/1e3, 'k--', label="MW $v_\mathrm{esc}$")

	pyp.xlabel(r"Target nuclear mass ($A$) [amu]")
	pyp.ylabel("Minimum velocity [km/s]")

	pyp.title(_title)

	return _fig

def Plot_Overlay_MinVelocity(_fig, _threshold_keV, _target, _dm, _c_opt):
	_x_pts = n.array([ _target.A ])
	_y_pts = n.array([ MinimumVelocity_ms(_threshold_keV, _target, _dm)/1e3])

	_fig.gca()
	_label = _target.Name
	pyp.plot(_x_pts, _y_pts, _c_opt, label=_label)

def Plot_LindhardFactor_axes():
	_fig = pyp.figure()

	_Er_array = n.linspace(5e-1, 20., num=100)
	_Lfactors = -1.0 * n.ones(len(_Er_array))

	pyp.semilogx(_Er_array,_Lfactors,'k--')

	pyp.xlim([5e-1 , 20.])
	pyp.ylim([0.0  , 1.0])

	pyp.xlabel("Recoil energy [keV]")
	pyp.ylabel("Lindhard factor")
	pyp.title("Fraction of energy that produces signal")

	return _fig, _Er_array

def Plot_Overlay_LindhardFactor(_fig, _Er_array, _target, _c_opt):
	_Lfactors = n.zeros(len(_Er_array))

	for i in n.arange(len(_Er_array)):
		_Lfactors[i] = _target.LindhardFactor(_Er_array[i])

	_label = _target.Name
	_fig.gca()
	pyp.semilogx(_Er_array,_Lfactors,_c_opt,label=_label)

def Plot_RecoilAngleDist_axes():
	_fig = pyp.figure()

	_theta_arr = n.linspace( -n.pi , n.pi , num=10000)
	_dist = -1.0 * n.ones(len(_theta_arr))

	pyp.plot(_theta_arr,_dist,'k--')

	pyp.xlim([-n.pi , n.pi])
	pyp.ylim([-0.1  , 1.1])

	pyp.xlabel("Angle of outgoing neutron")
	pyp.ylabel(r"$E_{R} / E_0$")
	pyp.title("Recoil energy fraction angular distribution")

	return _fig, _theta_arr

def Plot_RecoilAngleDist(_fig, _target, _theta_arr, _c_opt1, _c_opt2):
	_incoming_E_keV = 1.0 # Check in units of incoming energy
	_Efrac		= _target.RecoilEnergyAngularDist(_theta_arr)
	_Emax	    = _target.RecoilEnergyMax_AnyParticle_keV(m_neutron_GeV, _incoming_E_keV) * n.ones(len(_Efrac))
	
	if (_target.Name == "H"):
		_Efrac[ _theta_arr >= n.pi/2. ] = n.nan
		_Efrac[ _theta_arr <= -n.pi/2. ] = n.nan

	_label = _target.Name
	_fig.gca()
	pyp.plot(_theta_arr, _incoming_E_keV*_Efrac ,_c_opt1, label=_label)
	pyp.plot(_theta_arr, _incoming_E_keV*_Emax  ,_c_opt2)#,label=_label+" maximum recoil")

