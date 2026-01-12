"""
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
"""

#import Release.openCOSMORS # when running or debugging from within visual studio please use "import openCOSMORS" instead
import openCOSMORS
import numpy as np
import copy as cp
import os

options = {

    'sw_skip_COSMOSPACE_errors' : 0,                                # whether to skip COSMOSPACE errors in the case a parameter makes the equations unsolvable
                                                                    # this is handy when running optimizations of the parameters of COSMO-RS
                                                                    # you still need to catch the error and handle it approprietlywhether to skip COSMOSPACE errors in the case a parameter makes the equations unsolvable
    # input switches
    'sw_SR_COSMOfiles_type': 'ORCA_COSMO_TZVPD',                # ['Turbomole_COSMO_TZVP', 'Turbomole_COSMO_TZVPD_FINE', 'ORCA_COSMO_TZVPD']
    'sw_SR_combTerm': 1, # 0 No combinatorial term
                         # 1 to use the combinatorial term by Staverman-Guggenheim
                         # 2 to use the combinatorial term by Klamt (2003)
                         # 3 to use modified Staverman-Guggenheim combinatorial term with exponential scaling */

    # optional calculation switches
    'sw_SR_alwaysReloadSigmaProfiles': 0,
    'sw_SR_alwaysCalculateSizeRelatedParameters': 1,
    'sw_SR_useSegmentReferenceStateForInteractionMatrix': 0,        # [0, 1] 
                                                                    #       0 : conductor
                                                                    #       1 : pure segment
    'sw_SR_calculateContactStatisticsAndAdditionalProperties': 2,   # [0, 1, 2]
                                                                    #       0 : do not calculate anything additionaly
                                                                    #       1 : calculate contact statistics and average surface energies
                                                                    #       2 : calculate contact statistics, average surface energies and partial molar properties
    'sw_SR_partialInteractionMatrices' : [],                        # partial interaction matrices to be calculated as partial mlar properties
                                                                    # examples are ['E_mf', 'G_hb'], these must however also be added on the C++ side
    
                                                                  
    # segment descriptor switches
    'sw_SR_atomicNumber': 1,                                        # [0, 1] : differentiate between atomic numbers
    'sw_SR_misfit': 2,                                              # [0, 1, 2]
                                                                    #       0: do not use misfit correlation
                                                                    #       1 : use misfit correlation on all molecules
                                                                    #       2 : use misfit correlation only on neutral molecules
    'sw_SR_differentiateHydrogens' : 0,                             # [0, 1] : differentiate between hydrogen atoms depending on the heteroatom they are bound to
    'sw_SR_differentiateMoleculeGroups' : 0,  
    'sw_SR_polarizabilities' : 7                      # [0,1,3,5,6,7,8] : polarizability projections modes. Parametrizations and notations are from https://doi.org/10.1016/j.ces.2025.122170. 1 - Eq.(16)+𝑤1(Eq.(9)), 6 - Eq.(16)+𝑤5(Eq.(11)), 7 - Eq.(24)+𝑤5(Eq.(11)), 8 - Eq.(25)+𝑤5(Eq.(11)).  

}

# parameters (parameters below are the latest parametrization with polarizability projections: Eq.(24)+𝑤5(Eq.(11))) from https://doi.org/10.1016/j.ces.2025.122170. 
parameters = {  
    'Aeff': 4.90825,
    'alpha': 7876000.0,
    'CHB': 49318000.0,
    'CHBT': 1.5,
    'SigmaHB': 0.009953,
    'Rav': 0.5,
    'RavCorr': 1,
    'fCorr': 2.4,
    'comb_SG_z_coord': 0.0,
    'comb_SG_A_std': 1.0,
    'comb_modSG_exp': 2.0/3.0,
    'comb_lambda0': 0.463,
    'comb_lambda1': 0.42,
    'comb_lambda2': 0.065,
    'comb_SGG_lambda': 0.773,
    'comb_SGG_beta': 0.778,
    'm_vdW':29.567 ,
    'E_F_corr':346.82,

    'radii': {},

    'exp': {}
}

#  COSMOfiles
components = ['1112tetrachloroethane.orcacosmo',
                  'methanol.orcacosmo',
                  'water.orcacosmo']

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# structure for calculations
calculations = []

# add calculations
# calculation 1
calculation = {
    'concentrations' : np.array([[0.2, 0.5, 0.3]]),
    'temperatures' : np.array([298.15]),
    'components' : components,
    'reference_state_types' : np.array([3]),
    'reference_state_concentrations' : np.array([[]]),
    'component_indices' : [0, 1, 2]
}
calculations.append(calculation)

# calculation 2
calculation2 = cp.deepcopy(calculation)
calculation2['concentrations'] = np.array([[0.0, 1.0, 0.0]])
calculations.append(calculation2)

# calculation 3
calculation3 = cp.deepcopy(calculation)
calculation3['concentrations'] = np.array([[0.2, 0.5, 0.3]])
calculation3['reference_state_types'] = np.array([0])

calculations.append(calculation3)

# calculation 4
calculation4 = cp.deepcopy(calculation)
calculation4['concentrations'] = np.array([[0.2, 0.5, 0.3]])
calculation4['reference_state_concentrations'] = np.array([[0, 1.0, 0.]])
calculation4['reference_state_types'] = np.array([2])
calculations.append(calculation4)

# fill missing fields for each calculation
def fill_missing_calculation_structures(calculations, options):

    for i, calculation in enumerate(calculations):

        number_of_components = calculation['concentrations'].shape[1]
        number_of_concentrations = calculation['concentrations'].shape[0]

        calculation['ln_gamma_x_SR_combinatorial_calc'] = np.zeros((number_of_concentrations, number_of_components))
        calculation['ln_gamma_x_SR_residual_calc'] = np.zeros((number_of_concentrations, number_of_components))
        calculation['ln_gamma_x_SR_calc'] = np.zeros((number_of_concentrations, number_of_components))
        calculation['index'] = i

        if options['sw_SR_calculateContactStatisticsAndAdditionalProperties'] > 0:
            # one for Aij + number of  partial interaction matrices
            number_of_interaction_matrices = 1 + len(options['sw_SR_partialInteractionMatrices'])
            
            calculation['contact_statistics'] = np.zeros((number_of_concentrations, number_of_components, number_of_components))
            calculation['average_surface_energies'] = np.zeros((number_of_concentrations, number_of_interaction_matrices, number_of_components, number_of_components))

            if options['sw_SR_calculateContactStatisticsAndAdditionalProperties'] == 2:
                calculation['partial_molar_energies'] = np.zeros((number_of_concentrations, number_of_interaction_matrices, number_of_components))
    return calculations

calculations = fill_missing_calculation_structures(calculations, options)

# load the needed molecules, load the calculations and then do the calculation
try:
    openCOSMORS.loadMolecules(options, parameters, components)
    openCOSMORS.loadCalculations(calculations)
    calculations = openCOSMORS.calculate(parameters, calculations, False)

    for i in range(len(calculations)):
        print(i + 1)
        
        print('ln_gamma                      ', calculations[i]['ln_gamma_x_SR_residual_calc'])
        print('ln_gamma                      ', calculations[i]['ln_gamma_x_SR_combinatorial_calc'])
        if options['sw_SR_calculateContactStatisticsAndAdditionalProperties'] > 0:
            print('average_surface_energies      ', calculations[i]['average_surface_energies'][0,0,:,:].sum(1))
            
        if options['sw_SR_calculateContactStatisticsAndAdditionalProperties'] > 1:
            print('partial_molar_energies        ', calculations[i]['partial_molar_energies'][0,0,:])
except Exception as err:
    print('ERROR:')
    print(type(err))
    print(err.args)
