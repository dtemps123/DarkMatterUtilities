import numpy as n
import matplotlib.pyplot as pyp

## Import the modules
import DarkMatterUtilities.Constants as dmcon
import DarkMatterUtilities.DarkMatter as dmdm
import DarkMatterUtilities.Targets as dmt
import DarkMatterUtilities.FormFactors as dmff
import DarkMatterUtilities.Interactions as dmi

def ComparisonPlot_SI_elastic():
	## Create elastic dark matter models with various masses
	SI_sigma_per_nucleon_cm2 = 1.0e-45
	DM_mass_array_GeV        = n.array([5.0 , 10.0 , 50.0 , 1.0e2 , 1.0e3 , 1.0e4])
	N_DM_models              = len(DM_mass_array_GeV)
	DM_model_array           = n.zeros(N_DM_models, dtype=object)

	for i in n.arange(N_DM_models):
		this_model = dmdm.DarkMatter(DM_mass_array_GeV[i], SI_sigma_per_nucleon_cm2)
		DM_model_array[i] = this_model

	## Create a simple xenon target with a Helm form factor
	avg_xe_target = dmt.SimpleTarget( 131.293 , 54.0 , "Xe", 4.7808 )
	ff_helm       = dmff.FormFactor(avg_xe_target, 1)
	avg_xe_target.FormFactor = ff_helm

	## Create interactions for each of the DM models
	interaction_array        = n.zeros(N_DM_models, dtype=object)

	for i in n.arange(N_DM_models):
		this_int = dmi.Interaction(avg_xe_target, DM_model_array[i], _int_type="SI")
		interaction_array[i] = this_int

	## Calculate the differential scattering rate for each interaction
	E_recoil_range_keV = n.logspace(start=-1.0, stop=3.0, num=1000)
	output_data_array  = n.zeros(N_DM_models, dtype=n.ndarray)

	for i in n.arange(N_DM_models):
		this_int = interaction_array[i]
		output_data_array[i] = n.zeros(len(E_recoil_range_keV))
		for j in n.arange(len(E_recoil_range_keV)):
			this_Er_keV = E_recoil_range_keV[j]
			output_data_array[i][j] = this_int.DifferentialRate(this_Er_keV)

	## Create a plot of the models
	pyp.figure()

	for i in n.arange(N_DM_models):
		pyp.plot(E_recoil_range_keV, output_data_array[i], label=str(DM_model_array[i].Mass_GeV)+r" GeV$/c^2$")

	pyp.legend(loc='upper right')

	pyp.xscale('log')
	pyp.yscale('log')

	pyp.ylim([1.e-12 , 1.e-3])

	pyp.xlabel(r"Dark Matter Mass [GeV$/c^2$]")
	pyp.ylabel(r"Differential Rate [keV$^{-1}\times$ kg$^{-1}\times$ day$^{-1}$]")
	pyp.title(r"SI Elastic WIMP: $\sigma=$"+str(SI_sigma_per_nucleon_cm2)+r" cm$^2$")

def ComparisonPlot_SI_inelastic():
	## Create elastic dark matter models with various masses
	SI_sigma_per_nucleon_cm2 = 1.0e-45
	DM_mass_array_GeV        = n.array([50.0 , 50.0 ,  50.0 , 1.0e4 , 1.0e4 , 1.0e4])
	DM_delta_array_keV       = n.array([ 0.0 , 10.0 , 100.0 , 0.0   , 10.0  , 100.0])
	N_DM_models              = len(DM_mass_array_GeV)
	DM_model_array           = n.zeros(N_DM_models, dtype=object)

	for i in n.arange(N_DM_models):
		this_model = dmdm.DarkMatter(DM_mass_array_GeV[i], SI_sigma_per_nucleon_cm2, _delta_keV=DM_delta_array_keV[i])
		DM_model_array[i] = this_model

	## Create a simple xenon target with a Helm form factor
	avg_xe_target = dmt.SimpleTarget( 131.293 , 54.0 , "Xe", 4.7808 )
	ff_helm       = dmff.FormFactor(avg_xe_target, 1)
	avg_xe_target.FormFactor = ff_helm

	## Create interactions for each of the DM models
	interaction_array        = n.zeros(N_DM_models, dtype=object)

	for i in n.arange(N_DM_models):
		this_int = dmi.Interaction(avg_xe_target, DM_model_array[i], _int_type="SI")
		interaction_array[i] = this_int

	## Calculate the differential scattering rate for each interaction
	E_recoil_range_keV = n.logspace(start=-1.0, stop=3.0, num=1000)
	output_data_array  = n.zeros(N_DM_models, dtype=n.ndarray)

	for i in n.arange(N_DM_models):
		this_int = interaction_array[i]
		output_data_array[i] = n.zeros(len(E_recoil_range_keV))
		for j in n.arange(len(E_recoil_range_keV)):
			this_Er_keV = E_recoil_range_keV[j]
			output_data_array[i][j] = this_int.DifferentialRate(this_Er_keV)

	## Create a plot of the models
	pyp.figure()

	for i in n.arange(N_DM_models):
		pyp.plot(E_recoil_range_keV, output_data_array[i], 
			label=str(DM_model_array[i].Mass_GeV)+r" GeV$/c^2$"+"\t"+r"$\delta=$"+str(DM_model_array[i].delta_keV)+r" keV")

	pyp.legend(loc='lower left')

	pyp.xscale('log')
	pyp.yscale('log')

	pyp.ylim([1.e-12 , 1.e-3])

	pyp.xlabel(r"Dark Matter Mass [GeV$/c^2$]")
	pyp.ylabel(r"Differential Rate [keV$^{-1}\times$ kg$^{-1}\times$ day$^{-1}$]")
	pyp.title(r"SI Inelastic WIMP: $\sigma=$"+str(SI_sigma_per_nucleon_cm2)+r" cm$^2$")

ComparisonPlot_SI_elastic()
# ComparisonPlot_SI_inelastic()
pyp.show()