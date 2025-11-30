
#pragma once
#include <fstream>
#include <iostream>
#include "nlohmann/json.hpp"
#include "general.hpp"
#include "core_functions.hpp"

using json = nlohmann::json;
using namespace std;


void displayOnCLI(std::string message) {
    cout << message;
}

void displayTimeOnCLI(std::string message, unsigned long durationInMicroseconds) {
    cout << message << ": " << std::to_string(durationInMicroseconds) << " microseconds\n";
}

void initializeOnCLI() {

    display = displayOnCLI;
    displayTime = displayTimeOnCLI;

    initialize(param, false);
    param.sw_dGsolv_calculation_strict = 0;
    warnings = std::vector<std::string>();
    n_ex = 3;
}


void loadParametersOnCLI(const json& parameters) {

    param.Aeff = parameters["Aeff"].template get<double>();
    param.alpha = exp(parameters["ln_alpha"].template get<double>());

    param.CHB = exp(parameters["ln_CHB"].template get<double>());
    param.CHBT = parameters["CHBT"].template get<double>();
    param.SigmaHB = parameters["SigmaHB"].template get<double>();

    param.Rav = parameters["Rav"].template get<double>();

    if (param.sw_misfit > 0) {
        param.fCorr = parameters["fCorr"].template get<double>();
        param.RavCorr = parameters["RavCorr"].template get<double>();
    }

    if (param.sw_combTerm == 1 || param.sw_combTerm == 3) {
        param.comb_SG_A_std = parameters["comb_SG_A_std"].template get<double>();
        param.comb_SG_z_coord = parameters["comb_SG_z_coord"].template get<double>();

        if (param.sw_combTerm == 3)
            param.comb_modSG_exp = parameters["comb_modSG_exp"].template get<double>();
    }

    if (param.sw_combTerm == 2 || param.sw_combTerm == 5) {
        param.comb_lambda0 = parameters["comb_lambda0"].template get<double>();
        param.comb_lambda1 = parameters["comb_lambda1"].template get<double>();
        param.comb_lambda2 = parameters["comb_lambda2"].template get<double>();
    }

    if (param.sw_combTerm == 4) {
        param.comb_SGG_lambda = parameters["comb_SGG_lambda"].template get<double>();
        param.comb_SGG_beta = parameters["comb_SGG_beta"].template get<double>();
    }

    if (parameters.contains("dGsolv_eta")) {
        param.dGsolv_eta = parameters["dGsolv_eta"].template get<double>();
        param.dGsolv_omega_ring = parameters["dGsolv_omega_ring"].template get<double>();

        json radii = parameters["dGsolv_tau"];

        for (json::iterator it = radii.begin(); it != radii.end(); ++it)
            param.dGsolv_tau[std::stoi(it.key())] = it.value().template get<double>();

        json dGsolv_numberOfAtomsInRing = parameters["dGsolv_numberOfAtomsInRing"];
        for (int it : dGsolv_numberOfAtomsInRing)
            param.dGsolv_numberOfAtomsInRing.push_back(it);

        json dGsolv_E_gas = parameters["dGsolv_E_gas"];
        for (double it : dGsolv_E_gas)
            param.dGsolv_E_gas.push_back(it);
    }

    if (parameters.contains("radii")) {
        json radii = parameters["radii"];

        for (json::iterator it = radii.begin(); it != radii.end(); ++it)
            param.R_i[std::stoi(it.key())] = it.value().template get<double>();
    }
    // experimental parameters for prototyping
    if (parameters.contains("exp")) {
        json exp = parameters["exp"];
        for (json::iterator it = exp.begin(); it != exp.end(); ++it)
            param.exp_param[it.key()] = it.value().template get<double>();
    }

}

void loadMoleculesOnCLI(const json& options, const json& parameters, const json& componentPaths) {

    // options

    if (param.sw_calculateContactStatisticsAndAdditionalProperties != 0) {
        const json& partialInteractionMatrices = options["sw_SR_partialInteractionMatrices"];
        param.numberOfPartialInteractionMatrices = int(partialInteractionMatrices.size());
    }
    else {
        param.numberOfPartialInteractionMatrices = 0;
    }
    if (param.sw_misfit < 0 || param.sw_misfit > 2) {
        throw std::runtime_error("sw_SR_misfit should have one of the following values: [0, 1, 2].");
    }

    // parameters
    loadParametersOnCLI(parameters);

    for (const auto& componentPath : componentPaths) {
        molecule newMolecule = loadNewMolecule(param, componentPath.template get<std::string>());
        molecules.push_back(std::make_shared<molecule>(newMolecule));
    }

    if (molecules.size() == 0) {
        throw std::runtime_error("Please load at least one molecule.");
    }
}

void loadCalculationsOnCLI(const json& calculationsOnCLI) {

    const size_t numCalcs = calculationsOnCLI.size();

    if (numCalcs == 0) {
        throw std::runtime_error("Please specify at least one calculation.");
    }

    for (int i = 0; i < numCalcs; i++) {

        const json& calculationDict = calculationsOnCLI[i];

        // array of component indices
        const json& componentList = calculationDict["component_indices"];
        int numberOfComponents = int(componentList.size());

        calculation newCalculation(numberOfComponents);

        for (int j = 0; j < numberOfComponents; j++) {
            std::shared_ptr<molecule> thisMolecule = molecules[componentList[j].get<int>()];

            for (int k = 0; k < thisMolecule->segments.size(); k++) {
                newCalculation.segments.add((unsigned short)j, thisMolecule->segments.SegmentTypeGroup[k],
                    thisMolecule->segments.SegmentTypeSigma[k],
                    thisMolecule->segments.SegmentTypeSigmaCorr[k],
                    thisMolecule->segments.SegmentTypeHBtype[k],
                    thisMolecule->segments.SegmentTypeAtomicNumber[k],
                    thisMolecule->segments.SegmentTypeAtomicPolariz[k],
                    thisMolecule->segments.SegmentTypeAreas[k][0]);
            }

            newCalculation.components.push_back(thisMolecule);
        }
        newCalculation.segments.sort();
        newCalculation.segments.shrink_to_fit();

        // concentrations and temperatures
        auto temperatures = calculationDict["temperatures"].template get<std::vector<double>>();
        auto concentrations = calculationDict["concentrations"].template get<std::vector<std::vector<double>>>();

        for (int j = 0; j < concentrations.size(); j++) {
            std::vector<float> rowConcentration;

            float tempSumOfConcentrations = 0;
            for (int k = 0; k < numberOfComponents; k++) {
                float val = static_cast<float>(concentrations[j][k]);
                tempSumOfConcentrations += val;
                rowConcentration.push_back(val);
            }

            if (abs(1.0f - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
                throw std::runtime_error("For calculation number " + std::to_string(i) + ", the concentrations do not add up to unity. residual concentration: " + std::to_string(abs(1.0f - tempSumOfConcentrations)));
            }

            float temperature = static_cast<float>(temperatures[j]);

            newCalculation.temperatures.push_back(temperature);
            newCalculation.concentrations.push_back(rowConcentration);
        }

        std::vector<std::vector<double>> referenceStateConcentrations;
        if (calculationDict.contains("reference_state_concentrations")) {
            referenceStateConcentrations = calculationDict["reference_state_concentrations"].template get<std::vector<std::vector<double>>>();

            if (referenceStateConcentrations.size() != newCalculation.concentrations.size()) {
                throw std::runtime_error("concentrations and referenceStateConcentrations of calculation number " + std::to_string(i) + " have different sizes.");
            }
        }

        newCalculation.originalNumberOfCalculations = static_cast<unsigned short>(newCalculation.concentrations.size());

        // reference states
        auto referenceStateTypes = calculationDict["reference_state_types"].get<std::vector<int>>();
        for (int j = 0; j < referenceStateTypes.size(); j++) {
            int referenceStateType = referenceStateTypes[j];

            newCalculation.referenceStateType.push_back(static_cast<unsigned short>(referenceStateType));

            float tempSumOfConcentrations = 0;
            if (calculationDict.contains("reference_state_concentrations")) {
                for (int k = 0; k < referenceStateConcentrations[j].size(); k++) {
                    tempSumOfConcentrations += static_cast<float>(referenceStateConcentrations[j][k]);
                }
            }

            if (referenceStateType == 0) { // Pure component
                if (tempSumOfConcentrations != 0) {
                    throw std::runtime_error("A reference state concentration was specified for a calculation with reference state PureComponents, this does not make sense.");
                }

                std::vector<int> thisReferenceStateCalculationIndices;
                for (int k = 0; k < numberOfComponents; k++) {
                    std::vector<float> referenceStateConcentration;
                    for (int m = 0; m < numberOfComponents; m++) {
                        referenceStateConcentration.push_back(k == m ? 1.0f : 0.0f);
                    }

                    float temperature = static_cast<float>(temperatures[j]);
                    int referenceStateCalculationIndex = static_cast<int>(newCalculation.addOrFindArrayIndexForConcentration(referenceStateConcentration, temperature));
                    thisReferenceStateCalculationIndices.push_back(referenceStateCalculationIndex);
                }
                newCalculation.referenceStateCalculationIndices.push_back(thisReferenceStateCalculationIndices);
            }
            else if (referenceStateType == 1) { // untested: Pure component only neutral
                if (tempSumOfConcentrations != 0) {
                    throw std::runtime_error("A reference state concentration was specified for a calculation with reference state PureComponentsOnlyNeutral, this does not make sense.");
                }

                std::vector<int> thisReferenceStateCalculationIndices;
                for (int k = 0; k < numberOfComponents; k++) {
                    if (newCalculation.components[k]->moleculeCharge == 0) {
                        std::vector<float> referenceStateConcentration;
                        for (int m = 0; m < numberOfComponents; m++) {
                            referenceStateConcentration.push_back(k == m ? 1.0f : 0.0f);
                        }

                        float temperature = static_cast<float>(temperatures[j]);
                        int referenceStateCalculationIndex = static_cast<int>(newCalculation.addOrFindArrayIndexForConcentration(referenceStateConcentration, temperature));
                        thisReferenceStateCalculationIndices.push_back(referenceStateCalculationIndex);
                    }
                    else {
                        thisReferenceStateCalculationIndices.push_back(-1);
                    }
                }
                newCalculation.referenceStateCalculationIndices.push_back(thisReferenceStateCalculationIndices);
            }
            else if (referenceStateType == 2) { // Reference mixture
                std::vector<float> referenceStateConcentration;
                for (int m = 0; m < numberOfComponents; m++) {
                    referenceStateConcentration.push_back(static_cast<float>(referenceStateConcentrations[j][m]));
                }

                if (referenceStateConcentrations[j].size() != newCalculation.components.size()) {
                    throw std::runtime_error("A reference state concentration was specified with the wrong amount of concentrations.");
                }

                if (abs(1.0f - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
                    throw std::runtime_error("For calculation number " + std::to_string(i) + ", the reference concentrations do not add up to unity. residual concentration: " + std::to_string(abs(1.0f - tempSumOfConcentrations)));
                }

                float temperature = static_cast<float>(temperatures[j]);
                int referenceStateCalculationIndex = static_cast<int>(newCalculation.addOrFindArrayIndexForConcentration(referenceStateConcentration, temperature));

                std::vector<int> thisReferenceStateCalculationIndices;
                for (int m = 0; m < numberOfComponents; m++) {
                    thisReferenceStateCalculationIndices.push_back(referenceStateCalculationIndex);
                }
                newCalculation.referenceStateCalculationIndices.push_back(thisReferenceStateCalculationIndices);
            }
            else if (referenceStateType == 3 || referenceStateType == 4) { // COSMO or COSMO for solvation energy calculation
                if (tempSumOfConcentrations != 0) {
                    throw std::runtime_error("A reference state concentration was specified for a calculation with reference state COSMO, this does not make sense.");
                }

                std::vector<int> thisReferenceStateCalculationIndices(numberOfComponents, -1);
                newCalculation.referenceStateCalculationIndices.push_back(thisReferenceStateCalculationIndices);
            }
            else {
                throw std::runtime_error("An unknown reference state type was given.");
            }
        }

        // bind to matrices for it to work with the rest o the code

        newCalculation.lnGammaCombinatorial_data = Eigen::MatrixXf(
            int(newCalculation.originalNumberOfCalculations),
            int(newCalculation.components.size()));

        newCalculation.lnGammaCombinatorial_data.setZero();

        new (&newCalculation.lnGammaCombinatorial) Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            newCalculation.lnGammaCombinatorial_data.data(),
            int(newCalculation.originalNumberOfCalculations),
            int(newCalculation.components.size()));

        newCalculation.lnGammaResidual_data = Eigen::MatrixXf(
            int(newCalculation.originalNumberOfCalculations),
            int(newCalculation.components.size()));

        newCalculation.lnGammaResidual_data.setZero();

        new (&newCalculation.lnGammaResidual) Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            newCalculation.lnGammaResidual_data.data(),
            int(newCalculation.originalNumberOfCalculations),
            int(newCalculation.components.size()));

        newCalculation.lnGammaTotal_data = Eigen::MatrixXf(
            int(newCalculation.originalNumberOfCalculations),
            int(newCalculation.components.size()));

        newCalculation.lnGammaTotal_data.setZero();

        new (&newCalculation.lnGammaTotal) Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            newCalculation.lnGammaTotal_data.data(),
            int(newCalculation.originalNumberOfCalculations),
            int(newCalculation.components.size()));

        newCalculation.dGsolv_data = Eigen::MatrixXf(
            int(newCalculation.originalNumberOfCalculations),
            int(newCalculation.components.size()));

        newCalculation.dGsolv_data.setZero();

        new (&newCalculation.dGsolv) Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            newCalculation.dGsolv_data.data(),
            int(newCalculation.originalNumberOfCalculations),
            int(newCalculation.components.size()));

        if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {

            newCalculation.contactStatistics_data = Eigen::Tensor<float, 3, Eigen::RowMajor>(
                int(newCalculation.originalNumberOfCalculations),
                int(newCalculation.components.size()),
                int(newCalculation.components.size()));

            newCalculation.contactStatistics_data.setZero();

            new (&newCalculation.contactStatistics) Eigen::TensorMap<Eigen::Tensor<float, 3, Eigen::RowMajor>>(newCalculation.contactStatistics_data.data(),
                int(newCalculation.originalNumberOfCalculations),
                int(newCalculation.components.size()),
                int(newCalculation.components.size()));


            newCalculation.averageSurfaceEnergies_data = Eigen::Tensor<float, 4, Eigen::RowMajor>(
                int(newCalculation.originalNumberOfCalculations),
                int(param.numberOfPartialInteractionMatrices) + 1, // +1 because A_int is the first one
                int(newCalculation.components.size()),
                int(newCalculation.components.size()));

            newCalculation.averageSurfaceEnergies_data.setZero();

            new (&newCalculation.averageSurfaceEnergies) Eigen::TensorMap<Eigen::Tensor<float, 4, Eigen::RowMajor>>(newCalculation.averageSurfaceEnergies_data.data(),
                int(newCalculation.originalNumberOfCalculations),
                int(param.numberOfPartialInteractionMatrices) + 1, // +1 because A_int is the first one
                int(newCalculation.components.size()),
                int(newCalculation.components.size()));

            if (param.sw_calculateContactStatisticsAndAdditionalProperties == 2) {

                newCalculation.partialMolarEnergies_data = Eigen::Tensor<float, 3, Eigen::RowMajor>(
                    int(newCalculation.originalNumberOfCalculations),
                    int(param.numberOfPartialInteractionMatrices) + 1, // +1 because A_int is the first one
                    int(newCalculation.components.size()));

                newCalculation.partialMolarEnergies_data.setZero();

                new (&newCalculation.partialMolarEnergies) Eigen::TensorMap<Eigen::Tensor<float, 3, Eigen::RowMajor>>(newCalculation.partialMolarEnergies_data.data(),
                    int(newCalculation.originalNumberOfCalculations),
                    int(param.numberOfPartialInteractionMatrices) + 1, // +1 because A_int is the first one
                    int(newCalculation.components.size()));
            }

        }

        newCalculation.number = (int)i;
        finishCalculationInitiation(newCalculation);
        // the following std::move statement is very important as otherwise the Eigen::Map
        // are pointing to matrices deleted after exiting this function
        calculations.push_back(std::move(newCalculation));
    }
}


int main(int argc, char** argv)
{
    try {

        initializeOnCLI();
        std::string inputFilePath;
        std::string outputFilePath;
        if (argc < 2) {
            throw std::runtime_error("The required input json file path was not given.");
        }
        else if (argc == 2 || argc == 3) {
            inputFilePath = argv[1];
            if (!endsWith(inputFilePath, ".json")) {
                throw std::runtime_error("The required input json file path has to end in '.json'");
            }
            if (argc == 3) {
                outputFilePath = argv[2];
                if (!endsWith(inputFilePath, ".json")) {
                    throw std::runtime_error("The output json file path has to end in '.json'");
                }
            }
            else {
                outputFilePath = replace(inputFilePath, ".json", "_out.json");
            }
        }
        else {
            throw std::runtime_error("The openCOSMO-RS binary accepts two or three input arguments. More were given.");
        }

        std::ifstream f(inputFilePath);
        if (f.fail())
            throw std::runtime_error("The required input json file path was not found. Does it exists? Is the path correct?");
        json inputFileData = json::parse(f);

        loadMoleculesOnCLI(inputFileData, inputFileData, inputFileData["componentPaths"]);

        loadCalculationsOnCLI(inputFileData["calculations"]);

        std::vector<int> calculationIndices = {};

        for (int i = 0; i < inputFileData["calculations"].size(); i++)
            calculationIndices.push_back(i);

        calculate(calculationIndices);

        json outputJson = json::object();
        outputJson["dGsolv"] = json::array();

        if (warnings.size() > 0) {
            display("\nWARNINGS: \n");
            warnings.insert(warnings.begin(), "Some issues may lead to the calculated solvation energies having larger deviations than originally reported:");

            for (int i_warning = 0; i_warning < warnings.size(); i_warning++) {
                display(warnings[i_warning] + "\n");
            }
            display("\n");
        }
        outputJson["warnings"] = warnings;

        for (int calculationIndex = 0; calculationIndex < inputFileData["calculations"].size(); calculationIndex++) {
            json dGsolv_thisCalculation;
            for (int i = 0; i < calculations[calculationIndex].lnGammaTotal.rows(); i++) {
                std::vector<float> dGsolv;
                for (int j = 0; j < calculations[calculationIndex].lnGammaTotal.cols(); j++) {
                    dGsolv.push_back(calculations[calculationIndex].dGsolv(i, j));
                }
                json dGsolv_vec(dGsolv);
                dGsolv_thisCalculation.push_back(dGsolv_vec);
                
            }
            outputJson["dGsolv"].push_back(dGsolv_thisCalculation);
        }


        std::ofstream o(outputFilePath);
        o << std::setw(4) << outputJson << std::endl;

    }
    catch (const std::exception& e) {
        display("An error ocurred executing openCOSMO-RS:\n\n");
        display(e.what());
        display("\n");
        return 1;
    }

    return 0;
}