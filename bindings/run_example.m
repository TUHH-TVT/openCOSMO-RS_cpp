%{
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @authors: Simon Mueller & Andrés González de Castilla, 2022
%}

function run_example(mh)
    % general notes:
    % NOTE: MATLAB C++ Data API will not read 'text' as type MATLABString, but it will read "text" as
    % MATLABString. The mex code expects MATLABstrings, so please use the correct quotation marks below.
    
    if ~exist('mh','var')
        mh = mexhost;
    end
    
    % this is to ensure no allocated memory is still going to influence the execution.
    clear mex;

    %% options

    options.sw_skip_COSMOSPACE_errors = 0;                                  % whether to skip COSMOSPACE errors in the case a parameter makes the equations unsolvable
                                                                            % this is handy when running optimizations of the parameters of COSMO-RS
                                                                            % you still need to catch the error and handle it approprietly

    % input switches
    options.sw_SR_COSMOfiles_type = "ORCA_COSMO_TZVPD";                 % ['Turbomole_COSMO_TZVP', 'Turbomole_COSMO_TZVPD_FINE', 'ORCA_COSMO_TZVPD']
    options.sw_SR_combTerm = 2;

    % optional calculation switches
    options.sw_SR_alwaysReloadSigmaProfiles = 1;                            % Whether to reload the sigma profiles on every calculation
    options.sw_SR_alwaysCalculateSizeRelatedParameters = 1;                 % Always calculate the combinatorial term
    options.sw_SR_useSegmentReferenceStateForInteractionMatrix = 0;         % 0  =  conductor | 1  =  pure segment
    options.sw_SR_calculateContactStatisticsAndAdditionalProperties = 2;    % 0  =  do not calculate
                                                                            % 1  =  calculate contact statistics and average surface energies
                                                                            % 2  =  calculate contact statistics; average surface energies and partial molar properties
    options.sw_SR_partialInteractionMatrices = {};                          % e.g. {'E_mf', 'G_hb'}, these need to be specified in the C++ function "calculateInteractionMatrix"

    % segment descriptor switches
    options.sw_SR_atomicNumber = 0;                                         % [0, 1] : differentiate between atomic numbers
    options.sw_SR_misfit = 2;                                               % 0  =  do not use misfit correlation
                                                                            % 1  =  use misfit correlation on all molecules (ERROR)
                                                                            % 2  =  use misfit correlation only on neutral molecules
    options.sw_SR_differentiateHydrogens = 0;                               % [0, 1] : differentiate between hydrogen atoms depending on the heteroatom they are bound to
    options.sw_SR_differentiateMoleculeGroups = 0;                          % [0, 1] : differentiate between molecule groups

    %% parameters
    parameters.Aeff =  6.25;
    parameters.ln_alpha =  0;
    parameters.ln_CHB =  0;
    parameters.CHBT =  1.5;
    parameters.SigmaHB =  0.0085;
    parameters.Rav =  0.5;
    parameters.RavCorr =  1;
    parameters.fCorr =  2.4;
    parameters.comb_SG_z_coord =  10;
    parameters.comb_SG_A_std =  79.53;
    parameters.comb_modSG_exp =  2.0/3.0;
    parameters.comb_lambda0 =  0.463;
    parameters.comb_lambda1 =  0.42;
    parameters.comb_lambda2 =  0.065;
    parameters.comb_SGG_lambda =  0.773;
    parameters.comb_SGG_beta =  0.778;

    parameters.radii = ones(1,118)*NaN;

    parameters.exp = struct([]);

    %% COSMOfiles
    components = {"1112tetrachloroethane.orcacosmo", ...
                    "methanol.orcacosmo", ...
                    "water.orcacosmo"};
    cd(fileparts(which(mfilename)));
    
    %% tructure for calculations

    calculations = struct("concentrations", [],...
        "temperatures", [],...
        "components", [],...
        "reference_state_types", [],...
        "reference_state_concentrations", struct([]),...
        "ln_gamma_x_SR_combinatorial_calc",[],...
        "ln_gamma_x_SR_residual_calc", [],...
        "ln_gamma_x_SR_calc", [],...
        "component_indices", []);

    %% add calculations

    % calculation 1
    calculations(1).components = [1, 2, 3];
    calculations(1).temperatures = [298.15];
    calculations(1).concentrations = [0.2 0.5 0.3];
    calculations(1).reference_state_types = [3];
    calculations(1).reference_state_concentrations{end + 1} = [];

    % calculation 2
    calculations(2).components = [1, 2, 3];
    calculations(2).temperatures = [298.15];
    calculations(2).concentrations = [0 1 0];
    calculations(2).reference_state_types = [3];
    calculations(2).reference_state_concentrations{end + 1} = [];

    % calculation 3
    calculations(3).components = [1, 2, 3];
    calculations(3).temperatures = [298.15];
    calculations(3).concentrations = [0.2 0.5 0.3];
    calculations(3).reference_state_types = [0];
    calculations(3).reference_state_concentrations{end + 1} = [];

    % calculation 4
    calculations(4).components = [1, 2, 3];
    calculations(4).temperatures = [298.15];
    calculations(4).concentrations = [0.2 0.5 0.3];
    calculations(4).reference_state_types = [2];
    calculations(4).reference_state_concentrations{end + 1} = [0 1 0];

    %% fill missing fields for each calculation  
    function [calculations] = fill_missing_calculation_structures(calculations, options)
        for i = 1:length(calculations)
            x = size(calculations(i).concentrations);
            number_of_concentrations = x(1);
            number_of_components = x(2);
            calculations(i).ln_gamma_x_SR_combinatorial_calc = zeros(number_of_concentrations, number_of_components, 'single');
            calculations(i).ln_gamma_x_SR_residual_calc = zeros(number_of_concentrations, number_of_components, 'single');
            calculations(i).ln_gamma_x_SR_calc = zeros(number_of_concentrations, number_of_components, 'single');
            calculations(i).component_indices = calculations(i).components - 1;
            calculations(i).index = i - 1;

            if options.sw_SR_calculateContactStatisticsAndAdditionalProperties > 0
                % one for Aij + number of partial interaction matrices
                % defined above
                number_of_interaction_matrices = 1 + length(options.sw_SR_partialInteractionMatrices);
                calculations(i).contact_statistics = zeros(number_of_concentrations, number_of_components, number_of_components, 'single');
                calculations(i).average_surface_energies = zeros(number_of_concentrations, number_of_interaction_matrices, number_of_components, number_of_components, 'single');

                if options.sw_SR_calculateContactStatisticsAndAdditionalProperties == 2
                    calculations(i).partial_molar_energies = zeros(number_of_concentrations, number_of_interaction_matrices, number_of_components, 'single');
                end
            end
        end
    end

    calculations = fill_missing_calculation_structures(calculations, options);
    
    % load the needed molecules, load the calculations and then do the calculation
    feval(mh,"openCOSMORS", options, parameters, components);
    feval(mh,"openCOSMORS", options, parameters, calculations);
    calculations = feval(mh,"openCOSMORS", options, parameters, calculations);

    for i = 1:length(calculations)
        disp(i);
        disp('ln_gamma');
        disp(calculations(i).ln_gamma_x_SR_residual_calc);
        
        if options.sw_SR_calculateContactStatisticsAndAdditionalProperties > 0
            disp('average_surface_energies');
            disp(reshape(sum(calculations(i).average_surface_energies, 4), [1,3]));
            if options.sw_SR_calculateContactStatisticsAndAdditionalProperties > 1
                disp('partial_molar_energies');
                disp(reshape(calculations(i).partial_molar_energies(1,1,:), [1, 3]));
            end
        end
    end

end