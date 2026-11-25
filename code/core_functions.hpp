/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/


#pragma once
#include "interaction_matrix.hpp"
#include "contact_statistics.hpp"
#include "COSMOfile_functions.hpp"
#include "simd.hpp"
#include <stdexcept>
#include <functional>
#include <iostream>


#if defined(MEASURE_TIME) 
#include <atomic>
#include <chrono>
std::atomic<unsigned long> oneIteration_total_ms = 0;
std::chrono::high_resolution_clock::time_point oneIteration_last;
std::atomic<unsigned long> calculateTau_total_ms = 0;
std::atomic<unsigned long> calculateCOSMOSPACE_total_ms = 0;
std::atomic<unsigned long> calculateGammasForMolecules_total_ms = 0;
std::atomic<unsigned long> calculateContactStatistics_total_ms = 0;
std::atomic<unsigned long> rescaleSegments_total_ms = 0;
std::atomic<unsigned long> calculateCombinatorial_total_ms = 0;
std::atomic<unsigned long> calculateResidual_total_ms = 0;
std::atomic<unsigned long> addContributions_total_ms = 0;

void startCalculationMeasurement() {

    oneIteration_last = std::chrono::high_resolution_clock::now();

    // by commenting the following lines you achieve a cummulative sum of the time over more than one iteration.
    rescaleSegments_total_ms = 0;
    calculateCombinatorial_total_ms = 0;
    calculateResidual_total_ms = 0;
    calculateTau_total_ms = 0;
    calculateCOSMOSPACE_total_ms = 0;
    calculateContactStatistics_total_ms = 0;
    calculateGammasForMolecules_total_ms = 0;
    addContributions_total_ms = 0;
    oneIteration_total_ms = 0;
}

void stopCalculationMeasurement() {

    oneIteration_total_ms += (const unsigned long)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - oneIteration_last).count());

    // print accumulated times
    displayTime("rescaleSegments_total_ms               ", rescaleSegments_total_ms);
    displayTime("calculateCombinatorial_total_ms        ", calculateCombinatorial_total_ms);
    displayTime("calculateResidual_total_ms             ", calculateResidual_total_ms);
    displayTime("   calculateTau_total_ms               ", calculateTau_total_ms);
    displayTime("   calculateCOSMOSPACE_total_ms        ", calculateCOSMOSPACE_total_ms);
    if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {
        displayTime("   calculateContactStatistics_total_ms ", calculateContactStatistics_total_ms);
    }
    displayTime("   calculateGammasForMolecules_total_ms", calculateGammasForMolecules_total_ms);
    displayTime("addContributions_total_ms              ", addContributions_total_ms);
    displayTime("oneIteration_total_ms                  ", oneIteration_total_ms);
    display("\n");
}
#endif

void initialize(parameters& param, bool initializeParameters = true, bool initializeMolecules = true, bool initializeCalculations = true, bool showBinarySpecs = true) {


    if (initializeParameters) {
        n_ex = 0;
        param.ChargeRaster.clear();
        param.exp_param.clear();
        param.R_i = std::vector<double>(NUMBER_OF_ELEMENT_SLOTS, 0.0);
        param.R_i_COSMO = std::vector<double>(NUMBER_OF_ELEMENT_SLOTS, 0.0);
        param.HBClassElmnt = std::vector<int>(NUMBER_OF_ELEMENT_SLOTS, 0);

        // only set if polarizabilities are used, reset to not keep values from a previous run
        param.E_F_corr = 0.0;
        param.m_vdW = 0.0;

        // only set if solvation energies are calculated, reset to not keep values from a previous run
        param.dGsolv_tau = std::vector<double>(NUMBER_OF_ELEMENT_SLOTS, 0.0);
        param.dGsolv_E_gas.clear();
        param.dGsolv_numberOfAtomsInRing.clear();

        // Initialize hydrogen bond classes of the elements HBClassElmnt
        // 0 : only non HB  | 1 : potential donor  | 2 : potential acceptor | 3 : potential donor or acceptor
        // classify all hydrogens and some metals as potential donors and all others as potential acceptors.

        for (int atomic_number = 0; atomic_number < param.HBClassElmnt.size(); atomic_number++) {
            if (atomic_number <= 100) param.HBClassElmnt[atomic_number] = 2;
            else if (atomic_number > 100) param.HBClassElmnt[atomic_number] = 1;   // all hydrogens
        }

        // set some values manually
        param.HBClassElmnt[1] = 1;   // hydrogen
        param.HBClassElmnt[3] = 1;   // li 
        param.HBClassElmnt[4] = 1;   // be 
        param.HBClassElmnt[11] = 1;  // na 
        param.HBClassElmnt[12] = 1;  // mg 
        param.HBClassElmnt[13] = 1;  // al 
        param.HBClassElmnt[19] = 1;  // k* 
        param.HBClassElmnt[20] = 1;  // ca 
        param.HBClassElmnt[24] = 1;  // cr 
        param.HBClassElmnt[26] = 1;  // fe 
        param.HBClassElmnt[27] = 1;  // co 
        param.HBClassElmnt[29] = 1;  // cu 
        param.HBClassElmnt[30] = 1;  // zn 
        param.HBClassElmnt[37] = 1;  // rb 
        param.HBClassElmnt[38] = 1;  // sr 
        param.HBClassElmnt[48] = 1;  // cd 
        param.HBClassElmnt[55] = 1;  // cs 
        param.HBClassElmnt[56] = 1;  // ba 


        // initialize charge raster
        int steps = (int)((param.sigmaMax - param.sigmaMin) / param.sigmaStep + 1 + 0.00001);
        for (int i = 0; i < steps; i++) {
            param.ChargeRaster.push_back(param.sigmaMin + param.sigmaStep * i);
        }
    }

    if (initializeMolecules) {
        molecules.clear();
        monoatomicCationsResized = false;
        sigmaProfilesUsedInCalculation = false;
    }

    if (initializeCalculations)
        calculations.clear();

    if (showBinarySpecs)
        display("\nBINARY SPECS\n-------------------------\n" + compilation_mode + "\n" + OPENMP_parallelization + "\n" + vectorization_level + "\n" + precision + "\n-------------------------\n\n");
}

// dGsolv_E_gas and dGsolv_numberOfAtomsInRing are given per loaded molecule
void validateSolvationEnergyParameters(const parameters& param, int numberOfMolecules) {
    if (param.dGsolv_E_gas.size() == 0 && param.dGsolv_numberOfAtomsInRing.size() == 0) {
        return;
    }
    if (int(param.dGsolv_E_gas.size()) != numberOfMolecules) {
        throw std::runtime_error("dGsolv_E_gas has " + std::to_string(param.dGsolv_E_gas.size()) + " entries, but " + std::to_string(numberOfMolecules) +
            " molecules are loaded. One gas phase energy per loaded molecule (in the order of the component paths) is needed.");
    }
    if (int(param.dGsolv_numberOfAtomsInRing.size()) != numberOfMolecules) {
        throw std::runtime_error("dGsolv_numberOfAtomsInRing has " + std::to_string(param.dGsolv_numberOfAtomsInRing.size()) + " entries, but " + std::to_string(numberOfMolecules) +
            " molecules are loaded. One entry per loaded molecule (in the order of the component paths) is needed.");
    }
}

// checks the sizes of the input of one calculation before it is read
void validateCalculationSizes(int calculationNumber, int numberOfComponents, const std::vector<int>& componentIndices, int numberOfLoadedMolecules,
    int numberOfRows, int numberOfColumns, int numberOfTemperatures, int numberOfReferenceStateTypes, int numberOfReferenceStateRows) {

    const std::string prefix = "Calculation number " + std::to_string(calculationNumber) + ": ";

    if (numberOfComponents == 0) {
        throw std::runtime_error(prefix + "'component_indices' is empty.");
    }
    for (int index : componentIndices) {
        if (index < 0 || index >= numberOfLoadedMolecules) {
            throw std::runtime_error(prefix + "the component index " + std::to_string(index) + " is invalid, " + std::to_string(numberOfLoadedMolecules) + " molecules were loaded.");
        }
    }
    if (numberOfRows == 0) {
        throw std::runtime_error(prefix + "'concentrations' does not contain any row.");
    }
    if (numberOfRows > 65535) {
        throw std::runtime_error(prefix + "'concentrations' has " + std::to_string(numberOfRows) + " rows, at most 65535 rows (including the reference states) are possible per calculation. Please split the calculation.");
    }
    if (numberOfColumns != numberOfComponents) {
        throw std::runtime_error(prefix + "'concentrations' has " + std::to_string(numberOfColumns) + " columns, but " + std::to_string(numberOfComponents) + " components were given in 'component_indices'.");
    }
    if (numberOfTemperatures != numberOfRows) {
        throw std::runtime_error(prefix + "'temperatures' has " + std::to_string(numberOfTemperatures) + " entries, but 'concentrations' has " + std::to_string(numberOfRows) + " rows.");
    }
    if (numberOfReferenceStateTypes != numberOfRows) {
        throw std::runtime_error(prefix + "'reference_state_types' has " + std::to_string(numberOfReferenceStateTypes) + " entries, but 'concentrations' has " + std::to_string(numberOfRows) + " rows.");
    }
    if (numberOfReferenceStateRows != numberOfRows) {
        throw std::runtime_error(prefix + "'reference_state_concentrations' has " + std::to_string(numberOfReferenceStateRows) + " rows, but 'concentrations' has " + std::to_string(numberOfRows) + " rows.");
    }
}

// returns the index of the charge raster point left of sigma, throws if sigma is outside of the charge raster
int getLeftChargeRasterIndex(const parameters& param, double sigma, const std::string& description, const std::string& moleculeName) {

    // range check in floating point before converting to int (converting NaN or huge values to int is undefined)
    const double position = std::floor((sigma - param.sigmaMin) / param.sigmaStep);
    const int numberOfRasterPoints = int(param.ChargeRaster.size());

    int index = -1;
    if (std::isfinite(position) && position >= 0 && position <= double(numberOfRasterPoints - 2)) {
        index = int(position);
    }

    if (index < 0) {
        throw std::runtime_error("The " + description + " of a segment of molecule '" + moleculeName + "' is " + std::to_string(sigma) +
            " e/A^2, which is outside of the charge raster [" + std::to_string(param.sigmaMin) + ", " + std::to_string(param.sigmaMax) + "] e/A^2.");
    }

    return index;
}

void averageAndClusterSegments(parameters& param, molecule& _molecule, int approximateNumberOfSegmentTypes = 0) {

    // save reallocation time by specifying the approximate segment type number
    // this is possible whenever the sigma profile is reloaded and the previous number
    // is known
    if (approximateNumberOfSegmentTypes != 0) {
        // as default a vector resizes to the double of the current size when needed
        // here we are taking 110% of the approximateSegmentTypeNumber hopefully  
        // hindering a reallocation
        approximateNumberOfSegmentTypes = int(1.1 * approximateNumberOfSegmentTypes);
        _molecule.segments.reserve(approximateNumberOfSegmentTypes);
    }

    // average the segments
    int numberOfSegments = int(_molecule.segmentAreas.size());
    int numberOfAtoms = int(_molecule.atomAtomicNumbers.size());

    Eigen::VectorXd segmentRadiiSquared = _molecule.segmentAreas / PI;
    double RavSquared = param.Rav * param.Rav;
    double RavCorrSquared = param.RavCorr * param.RavCorr;


    Eigen::VectorXd averagedSigmas = Eigen::VectorXd::Zero(numberOfSegments);
    Eigen::VectorXd averagedSigmaCorrs;

    bool calculateMisfitCorrelation = false;
    if (param.sw_misfit == 1) {
        calculateMisfitCorrelation = true;
        averagedSigmaCorrs = Eigen::VectorXd::Zero(numberOfSegments);
    }
    else if (param.sw_misfit == 2 && (_molecule.moleculeCharge == 0 || numberOfAtoms > 1)) {
        calculateMisfitCorrelation = true;
        averagedSigmaCorrs = Eigen::VectorXd::Zero(numberOfSegments);
    }

    double distanceSegmentSegmentSquared;
    double temporaryValue = 0, multiplyingFactor = 0;

    for (int segmentIndexI = 0; segmentIndexI < numberOfSegments; segmentIndexI++) {

        double runningTotalSigmas = 0, runningTotalSigmaCorrs = 0;

        for (int segmentIndexJ = 0; segmentIndexJ < numberOfSegments; segmentIndexJ++) {
            distanceSegmentSegmentSquared = (_molecule.segmentPositions(segmentIndexI, Eigen::indexing::all) - _molecule.segmentPositions(segmentIndexJ, Eigen::indexing::all)).array().square().matrix().sum();

            temporaryValue = segmentRadiiSquared(segmentIndexJ) + RavSquared;
            multiplyingFactor = ((segmentRadiiSquared(segmentIndexJ) * RavSquared) / temporaryValue) * exp(-distanceSegmentSegmentSquared / temporaryValue);

            runningTotalSigmas += multiplyingFactor;
            averagedSigmas(segmentIndexI) = averagedSigmas(segmentIndexI) + _molecule.segmentSigmas(segmentIndexJ) * multiplyingFactor;

            if (calculateMisfitCorrelation == true) {
                temporaryValue = segmentRadiiSquared(segmentIndexJ) + RavCorrSquared;
                multiplyingFactor = ((segmentRadiiSquared(segmentIndexJ) * RavCorrSquared) / temporaryValue) * exp(-distanceSegmentSegmentSquared / temporaryValue);

                runningTotalSigmaCorrs += multiplyingFactor;
                averagedSigmaCorrs(segmentIndexI) = averagedSigmaCorrs(segmentIndexI) + _molecule.segmentSigmas(segmentIndexJ) * multiplyingFactor;
            }
        }
        averagedSigmas(segmentIndexI) = averagedSigmas(segmentIndexI) / runningTotalSigmas;

        if (calculateMisfitCorrelation == true)
            averagedSigmaCorrs(segmentIndexI) = averagedSigmaCorrs(segmentIndexI) / runningTotalSigmaCorrs;
    }

    bool calculateSolvationEnergies = param.dGsolv_E_gas.size() > 0;
    if (calculateSolvationEnergies){
        if (_molecule.qmMethod != "DFT_CPCM_BP86_def2-TZVP+def2-TZVPD_SP" && _molecule.qmMethod != "DFT_BP86_def2-TZVPD_SP"){
            if (param.sw_dGsolv_calculation_strict == 1) {
                throw std::runtime_error("The QSPR model for the molar volume only works for the quantum chemistry method DFT_BP86_def2-TZVPD_SP");
            }
            else {
                addWarning(" - The QSPR model for the molar volume was parametrized using a different quantum chemistry method than the one you are using. Recommended method: DFT_BP86_def2-TZVPD_SP");
            }
        }

        int numberOfSiAtoms = 0;
        int numberOfHAtoms = 0;
        int numberOfOAtoms = 0;
        for (int i = 0; i < numberOfAtoms; i++){
            if (_molecule.atomAtomicNumbers[i] == 14)
                numberOfSiAtoms += 1;
            else if (_molecule.atomAtomicNumbers[i] == 1 || _molecule.atomAtomicNumbers[i] > 100)
                numberOfHAtoms += 1;
            else if (_molecule.atomAtomicNumbers[i] == 8)
                numberOfOAtoms += 1;
        }

        if (numberOfAtoms == 3 && numberOfHAtoms == 2 && numberOfOAtoms == 1){
            _molecule.molarVolumeAt25C = 18.06863632;
        } else {
            Eigen::ArrayXd averagedSigmasSquared = averagedSigmas.array() * averagedSigmas.array();
            double secondSigmaMoment = (averagedSigmasSquared * _molecule.segmentAreas.array()).sum() * 10000;
            double fourthSigmaMoment = (averagedSigmasSquared * averagedSigmasSquared * _molecule.segmentAreas.array()).sum() * 100000000;

            _molecule.molarVolumeAt25C = 0.9430785419976806 * numberOfAtoms \
                + 0.6977322963011842 * _molecule.Area \
                - 0.3161763939689293 * secondSigmaMoment \
                + 0.032441059832647084 * fourthSigmaMoment \
                + 8.113026329415828 * numberOfSiAtoms \
                - 0.07066832029215675;
        }
    }

    // determine the hydrogen bonding type based on the atomic number

    for (int i = 0; i < numberOfSegments; i++) {

        int atomicNumber = _molecule.segmentAtomicNumber[i];
        // 0 : non HB  | 1 : potential donor  | 2 : potential acceptor | 3 : potential donor or acceptor
        switch (param.HBClassElmnt[atomicNumber]) {
        case 0:
            _molecule.segmentHydrogenBondingType(i) = 0;
            break;
        case 1:
            _molecule.segmentHydrogenBondingType(i) = averagedSigmas(i) < 0 ? 1 : 0;
            break;
        case 2:
            _molecule.segmentHydrogenBondingType(i) = averagedSigmas(i) < 0 ? 0 : 2;
            break;
        case 3:
            _molecule.segmentHydrogenBondingType(i) = averagedSigmas(i) < 0 ? 1 : 2;
            break;
        default:
            throw std::runtime_error("Unknown HB classification value used.");
        }
    }

    // cluster segments into segment types
    double sigmaLeft = -1;
    double sigmaRight = -1;
    double AsigmaLeft = -1;
    double AsigmaRight = -1;

    double sigmaCorrLeft = -1;
    double sigmaCorrRight = -1;

    double AsigmaLeftSigmaCorrLeft = -1;
    double AsigmaLeftSigmaCorrRight = -1;
    double AsigmaRightSigmaCorrLeft = -1;
    double AsigmaRightSigmaCorrRight = -1;

    int ind_SigmaCorr_left = 0;

    for (int j = 0; j < numberOfSegments; j++) {

        int ind_Sigma_left = getLeftChargeRasterIndex(param, averagedSigmas(j), "averaged screening charge density (sigma)", _molecule.name);

        sigmaLeft = param.ChargeRaster[ind_Sigma_left];
        sigmaRight = param.ChargeRaster[ind_Sigma_left + 1];

        AsigmaRight = _molecule.segmentAreas(j) * (averagedSigmas(j) - sigmaLeft) / param.sigmaStep;
        AsigmaLeft = _molecule.segmentAreas(j) * (sigmaRight - averagedSigmas(j)) / param.sigmaStep;


        if (calculateMisfitCorrelation == true) {
            ind_SigmaCorr_left = getLeftChargeRasterIndex(param, averagedSigmaCorrs(j), "correlation screening charge density (sigmaCorr)", _molecule.name);

            sigmaCorrLeft = param.ChargeRaster[ind_SigmaCorr_left];
            sigmaCorrRight = param.ChargeRaster[ind_SigmaCorr_left + 1];

            AsigmaLeftSigmaCorrRight = AsigmaLeft * (averagedSigmaCorrs(j) - sigmaCorrLeft) / param.sigmaStep;
            AsigmaLeftSigmaCorrLeft = AsigmaLeft * (sigmaCorrRight - averagedSigmaCorrs(j)) / param.sigmaStep;
            AsigmaRightSigmaCorrRight = AsigmaRight * (averagedSigmaCorrs(j) - sigmaCorrLeft) / param.sigmaStep;
            AsigmaRightSigmaCorrLeft = AsigmaRight * (sigmaCorrRight - averagedSigmaCorrs(j)) / param.sigmaStep;
        }


        unsigned short atomicNumber = _molecule.segmentAtomicNumber(j);
        if (param.sw_atomicNumber == 0) {
            atomicNumber = 0;
        }

        // for monoatomic ions or if misfit correlation is deactivated
        if (_molecule.moleculeGroup == 3 || _molecule.moleculeGroup == 5 || calculateMisfitCorrelation == false) {

            _molecule.segments.add(0, _molecule.moleculeGroup, sigmaLeft, 0.0f, _molecule.segmentHydrogenBondingType(j), atomicNumber, _molecule.segmentAtomicPolariz(j), AsigmaLeft);
            _molecule.segments.add(0, _molecule.moleculeGroup, sigmaRight, 0.0f, _molecule.segmentHydrogenBondingType(j), atomicNumber, _molecule.segmentAtomicPolariz(j), AsigmaRight);
        }
        else {
            _molecule.segments.add(0, _molecule.moleculeGroup, sigmaLeft, sigmaCorrLeft, _molecule.segmentHydrogenBondingType(j), atomicNumber, _molecule.segmentAtomicPolariz(j), AsigmaLeftSigmaCorrLeft);
            _molecule.segments.add(0, _molecule.moleculeGroup, sigmaLeft, sigmaCorrRight, _molecule.segmentHydrogenBondingType(j), atomicNumber, _molecule.segmentAtomicPolariz(j), AsigmaLeftSigmaCorrRight);
            _molecule.segments.add(0, _molecule.moleculeGroup, sigmaRight, sigmaCorrLeft, _molecule.segmentHydrogenBondingType(j), atomicNumber, _molecule.segmentAtomicPolariz(j), AsigmaRightSigmaCorrLeft);
            _molecule.segments.add(0, _molecule.moleculeGroup, sigmaRight, sigmaCorrRight, _molecule.segmentHydrogenBondingType(j), atomicNumber, _molecule.segmentAtomicPolariz(j), AsigmaRightSigmaCorrRight);
        }
    }
}


inline Eigen::MatrixXd extractAtomicPolarizabilityTensor(molecule& _molecule, int atomIndex) {

    Eigen::MatrixXd polarizabilityTensors = _molecule.atomPolarizabilityTensors;
    Eigen::MatrixXd polarizabilityTensor(3,3);
    Eigen::VectorXd polarizabilityTensorLine = polarizabilityTensors.row(atomIndex);

    polarizabilityTensor(0,0) = polarizabilityTensorLine(0);
    polarizabilityTensor(1,1) = polarizabilityTensorLine(1);
    polarizabilityTensor(2,2) = polarizabilityTensorLine(2);
    polarizabilityTensor(0,1) = polarizabilityTensorLine(3);
    polarizabilityTensor(1,0) = polarizabilityTensorLine(3);
    polarizabilityTensor(0,2) = polarizabilityTensorLine(4);
    polarizabilityTensor(2,0) = polarizabilityTensorLine(4);
    polarizabilityTensor(1,2) = polarizabilityTensorLine(5);
    polarizabilityTensor(2,1) = polarizabilityTensorLine(5);

    return polarizabilityTensor;
}

void projectSegments(parameters& param, molecule& _molecule) {
    // project polarizability tensors onto segments and scale them

    int numberOfSegments = int(_molecule.segmentAreas.size());
    int numberOfAtoms = int(_molecule.atomAtomicNumbers.size());

    double polarizability_projection;

    Eigen::MatrixXd surfaceCoordinatesInAu = _molecule.segmentPositions / 0.529177249;
    Eigen::VectorXd surfaceAreasInAu = _molecule.segmentAreas / pow(0.529177249, 2);
    Eigen::VectorXd surfaceRadiiSquared = _molecule.segmentAreas / PI;
    Eigen::MatrixXd atomCoordinatesInAu = _molecule.atomPositions / 0.529177249;
    Eigen::VectorXi segmentAtomIndex = _molecule.segmentAtomIndices;
    Eigen::MatrixXd polarizabilityTensor(3, 3);
    Eigen::VectorXd areaOfAtom = Eigen::VectorXd::Zero(numberOfAtoms);
    Eigen::VectorXd segmentsperAtom = Eigen::VectorXd::Zero(numberOfAtoms);
    Eigen::VectorXd distanceSegmentAtom(3);

    if (param.sw_SR_polarizabilities == 0) {
        for (int segmentIndexI = 0; segmentIndexI < numberOfSegments; segmentIndexI++) {
            _molecule.segmentAtomicPolariz(segmentIndexI) = 0;
        }
    }
    else {
        const int mode = param.sw_SR_polarizabilities;
        if (mode != 1 && mode != 3 && mode != 5 && mode != 6 && mode != 7 && mode != 8) {
            throw std::runtime_error("sw_SR_polarizabilities = " + std::to_string(mode) + " is not supported. Possible values are: 0, 1, 3, 5, 6, 7, 8.");
        }

        if (_molecule.atomPolarizabilityTensors.rows() != numberOfAtoms || _molecule.atomPolarizabilityTensors.cols() != 6) {
            throw std::runtime_error("The COSMO file of molecule '" + _molecule.name + "' does not contain atomic polarizabilities for all " + std::to_string(numberOfAtoms) +
                " atoms, which are needed for sw_SR_polarizabilities > 0. Please use ORCA COSMO files containing the section '#ATOMIC POLARIZABILITIES'.");
        }

        for (int segmentIndexI = 0; segmentIndexI < numberOfSegments; segmentIndexI++) {
            int atomIndex = segmentAtomIndex(segmentIndexI);
            areaOfAtom(atomIndex) = areaOfAtom(atomIndex) + surfaceAreasInAu(segmentIndexI);
            segmentsperAtom(atomIndex) = segmentsperAtom(atomIndex) + 1;
        }


        for (int segmentIndexI = 0; segmentIndexI < numberOfSegments; segmentIndexI++) {
            int atomIndex = segmentAtomIndex(segmentIndexI);
            for (int k = 0; k < 3; k++) {
                distanceSegmentAtom(k) = std::abs(surfaceCoordinatesInAu(segmentIndexI, k) - atomCoordinatesInAu(atomIndex, k));
            }
            double distanceSegmentAtom_norm = distanceSegmentAtom(0) + distanceSegmentAtom(1) + distanceSegmentAtom(2);
            distanceSegmentAtom = distanceSegmentAtom / distanceSegmentAtom_norm;
            polarizabilityTensor = extractAtomicPolarizabilityTensor(_molecule, atomIndex);

            Eigen::RowVector3d distanceSegmentAtom_row = distanceSegmentAtom.transpose();
            Eigen::RowVector3d tmp_vec = distanceSegmentAtom_row * polarizabilityTensor;

            polarizability_projection = tmp_vec.norm();

            // The projections are clustered by rounding on purpose (Grigorash et al., Chem. Eng. Sci. 319 (2026) 122170,
            // Sect. 3.1.3 and Table 3): grid step 10 for unscaled projections and the scaling w5 (modes 5 to 8),
            // grid step 0.1 for the scalings w1 (mode 1) and w3 (mode 3). The parameters were fitted with these grids.
            // round_and_truncate(x, n) rounds to multiples of 10^n.
            if (param.sw_SR_polarizabilities == 5) {
                polarizability_projection = round_and_truncate(polarizability_projection, 1);
            }
            else if (param.sw_SR_polarizabilities == 1) {
                polarizability_projection = round_and_truncate(polarizability_projection / segmentsperAtom(atomIndex), -1);
            }
            else if (param.sw_SR_polarizabilities == 3) {
                polarizability_projection = round_and_truncate(polarizability_projection * surfaceAreasInAu(segmentIndexI) / areaOfAtom(atomIndex), -1);
            }
            _molecule.segmentAtomicPolariz(segmentIndexI) = polarizability_projection;
        }

        if (param.sw_SR_polarizabilities == 6 || param.sw_SR_polarizabilities == 7 || param.sw_SR_polarizabilities == 8) {
            Eigen::VectorXd segmentAtomicPolarizCopy = _molecule.segmentAtomicPolariz;
            for (int segmentIndexI = 0; segmentIndexI < numberOfSegments; segmentIndexI++) {
                Eigen::VectorXd distanceSegmentAllSegments = Eigen::VectorXd::Zero(numberOfSegments);
                Eigen::VectorXd averagingFactors = Eigen::VectorXd::Zero(numberOfSegments);
                double RavSquared = param.Rav * param.Rav;
                distanceSegmentAllSegments(Eigen::indexing::all) = 0.52917721067121 * (surfaceCoordinatesInAu(segmentIndexI, Eigen::indexing::all).replicate(numberOfSegments, 1) - surfaceCoordinatesInAu).array().square().rowwise().sum().sqrt();
                averagingFactors = ((surfaceRadiiSquared * RavSquared).array() / (surfaceRadiiSquared.array() + RavSquared)) * exp(-1 * distanceSegmentAllSegments.array().square() / (surfaceRadiiSquared.array() + RavSquared));
                averagingFactors = averagingFactors / averagingFactors.sum();
                polarizability_projection = (segmentAtomicPolarizCopy.array() * averagingFactors.array()).sum();
                _molecule.segmentAtomicPolariz(segmentIndexI) = round_and_truncate(polarizability_projection, 1);
            }
        }
    }
}

molecule loadNewMolecule(parameters& param, std::string componentPath) {

    molecule newMolecule;
    if (param.sw_COSMOfiles_type == "Turbomole_COSMO_TZVP" || param.sw_COSMOfiles_type == "Turbomole_COSMO_TZVPD_FINE") {
        newMolecule = getMoleculeFromTurbomoleCOSMOfile(componentPath);
    }
    else if (param.sw_COSMOfiles_type == "ORCA_COSMO_TZVPD") {
        newMolecule = getMoleculeFromORCACOSMOfile(componentPath);
    }
    else {
        throw std::runtime_error("No method for reading COSMOfiles has been implemented for the following type: " + param.sw_COSMOfiles_type);
    }

    std::string componentName = componentPath;
    std::vector<std::string> parts = split(componentName, '\\');
    parts = split(parts[parts.size() - 1], '/');
    parts = split(parts[parts.size() - 1], '.');
    newMolecule.name = trim(parts[0]);

    // if sigma profile not yet loaded
    // check the data read from the COSMO file
    {
        const int numberOfAtoms = int(newMolecule.atomAtomicNumbers.size());
        const int numberOfSegments = int(newMolecule.segmentAreas.size());

        if (numberOfAtoms == 0) {
            throw std::runtime_error("No atoms were found in the COSMO file: " + componentPath);
        }
        if (numberOfSegments == 0) {
            throw std::runtime_error("No surface segments were found in the COSMO file: " + componentPath);
        }
        if (!std::isfinite(newMolecule.Area) || newMolecule.Area <= 0 || !std::isfinite(newMolecule.Volume) || newMolecule.Volume <= 0) {
            throw std::runtime_error("The cavity area and volume could not be read from the COSMO file: " + componentPath);
        }
        for (int i = 0; i < numberOfSegments; i++) {
            const int atomIndex = newMolecule.segmentAtomIndices(i);
            if (atomIndex < 0 || atomIndex >= numberOfAtoms) {
                throw std::runtime_error("Surface segment " + std::to_string(i) + " refers to the atom index " + std::to_string(atomIndex) +
                    ", but the molecule has " + std::to_string(numberOfAtoms) + " atoms. COSMO file: " + componentPath);
            }
            if (!(newMolecule.segmentAreas(i) > 0) || !std::isfinite(newMolecule.segmentSigmas(i))) {
                throw std::runtime_error("Surface segment " + std::to_string(i) + " has an invalid area or charge. COSMO file: " + componentPath);
            }
        }
    }

    if (newMolecule.segments.size() == 0) {
        double sumOfScreeningCharge = (newMolecule.segmentAreas.array() * newMolecule.segmentSigmas.array()).matrix().sum();
        newMolecule.moleculeCharge = (signed char)(std::round(-1.0f * sumOfScreeningCharge));

        int numberOfAtoms = int(newMolecule.atomAtomicNumbers.size());

        // Store atomic radii and check for consistency
        for (int atomIndex = 0; atomIndex < numberOfAtoms; atomIndex++) {
            int AtomicNumber = newMolecule.atomAtomicNumbers(atomIndex);
            if (param.R_i_COSMO[AtomicNumber] != 0 && newMolecule.atomRadii(atomIndex) != param.R_i_COSMO[AtomicNumber]) {
                throw std::runtime_error("Inconsistent radii set for atomic number " + std::to_string(AtomicNumber) + " was found.");
            }
            else if (param.R_i_COSMO[AtomicNumber] == 0 && newMolecule.atomRadii(atomIndex) != 0) {
                param.R_i_COSMO[AtomicNumber] = newMolecule.atomRadii(atomIndex);
            }
        }

        // Calculate distance between atoms
        Eigen::MatrixXd distanceAtomAtomSquared(numberOfAtoms, numberOfAtoms);

        for (int atomIndex = 0; atomIndex < numberOfAtoms; atomIndex++) {
            Eigen::MatrixXd thisAtomPosition = newMolecule.atomPositions(atomIndex, Eigen::indexing::all);
            distanceAtomAtomSquared(Eigen::indexing::all, atomIndex) = (newMolecule.atomPositions - thisAtomPosition.replicate(numberOfAtoms, 1)).array().square().rowwise().sum();
        }

        // set moleculeGroup
        int moleculeGroup = 0;

        if (param.sw_differentiateMoleculeGroups == 1) {

            moleculeGroup = -1;

            if (numberOfAtoms == 3) {
                int numberOfFoundOxygens = 0;
                int numberOfFoundHydrogens = 0;
                for (int atomIndex = 0; atomIndex < numberOfAtoms; atomIndex++) {
                    if (newMolecule.atomAtomicNumbers(atomIndex) == 1) {
                        numberOfFoundHydrogens += 1;
                    }
                    else if (newMolecule.atomAtomicNumbers(atomIndex) == 8) {
                        numberOfFoundOxygens += 1;
                    }
                }

                if (numberOfFoundOxygens == 1 && numberOfFoundHydrogens == 2) {
                    moleculeGroup = 2;
                }
            }

            if (moleculeGroup == -1) {
                if (newMolecule.moleculeCharge == 0) {
                    moleculeGroup = numberOfAtoms == 1 ? 0 : 1;
                }
                else if (newMolecule.moleculeCharge > 0) {
                    moleculeGroup = numberOfAtoms == 1 ? 3 : 4;
                }
                else if (newMolecule.moleculeCharge < 0) {
                    moleculeGroup = numberOfAtoms == 1 ? 5 : 6;
                }
            }
        }
        newMolecule.moleculeGroup = moleculeGroup;

        // this changes the atomic number of a hydrogen atom to 100 + the atomic number of the closest heavy atom
        // giving the abbility to differentiate between hydrogen atoms depending on the atom they are bound to.
        if (param.sw_differentiateHydrogens == 1) {

            // the original atomic numbers are used, so that already relabeled hydrogens are still
            // recognized as hydrogens (e.g. H2 gets 101 for both hydrogen atoms)
            const Eigen::VectorXi originalAtomicNumbers = newMolecule.atomAtomicNumbers;

            for (int atomIndexI = 0; atomIndexI < numberOfAtoms; atomIndexI++) {
                double minimumDistanceAtomIAtomJSquared = 10e14;
                int closestAtomIndex = -1;
                if (originalAtomicNumbers(atomIndexI) != 1) {
                    continue;
                }

                for (int atomIndexJ = 0; atomIndexJ < numberOfAtoms; atomIndexJ++) {

                    if (atomIndexI == atomIndexJ) {
                        continue;
                    }

                    // search for next atom that is not a Hydrogen unless the molecule is H2
                    if (originalAtomicNumbers(atomIndexJ) == 1 && numberOfAtoms != 2) {
                        continue;
                    }

                    if (distanceAtomAtomSquared(atomIndexI, atomIndexJ) < minimumDistanceAtomIAtomJSquared) {
                        minimumDistanceAtomIAtomJSquared = distanceAtomAtomSquared(atomIndexI, atomIndexJ);
                        closestAtomIndex = atomIndexJ;
                    }
                }
                if (closestAtomIndex == -1) {
                    throw std::runtime_error("sw_SR_differentiateHydrogens = 1: no bonding partner was found for hydrogen atom " + std::to_string(atomIndexI) + " of molecule '" + newMolecule.name + "'.");
                }
                int newAtomicNumber = 100 + originalAtomicNumbers(closestAtomIndex);
                newMolecule.atomAtomicNumbers(atomIndexI) = newAtomicNumber;
                param.R_i_COSMO[newAtomicNumber] = param.R_i_COSMO[1];
            }
        }
        else if (param.sw_differentiateHydrogens != 0) {
            throw std::runtime_error("differentiateHydrogens accepts values [0, 1]");
        }

        int numberOfSegments = int(newMolecule.segmentAreas.size());
        newMolecule.segmentAtomicNumber = Eigen::VectorXi(numberOfSegments);

        for (int i = 0; i < numberOfSegments; i++) {

            int atomicNumber = newMolecule.atomAtomicNumbers(newMolecule.segmentAtomIndices(i));
            newMolecule.segmentAtomicNumber(i) = atomicNumber;
        }

        newMolecule.segmentHydrogenBondingType = Eigen::VectorXi(numberOfSegments);
        newMolecule.segmentAtomicPolariz = Eigen::VectorXd(numberOfSegments);

        projectSegments(param,newMolecule);
        averageAndClusterSegments(param, newMolecule);
    }

    newMolecule.clear_unneeded_matrices(param.sw_alwaysReloadSigmaProfiles);
    newMolecule.segments.shrink_to_fit();

#ifdef DEBUG_INFO
    newMolecule.segments.sort();
    WriteExtendedSigmaProfileToFile(newMolecule.name + ".extsp", newMolecule.segments);
#endif

    return newMolecule;
}

void reloadAllMolecules() {
    threadException e;
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int i = 0; i < molecules.size(); i++) {
        e.run([=] {
            // reference: the reloaded sigma profile has to be stored in the loaded molecule
            molecule& _molecule = *molecules[i];
            int previousNumberOfSegmentTypes = int(_molecule.segments.size());
            _molecule.segments.clear();
            // the polarizability projections of modes 6 to 8 depend on Rav as well
            projectSegments(param, _molecule);
            averageAndClusterSegments(param, _molecule, previousNumberOfSegmentTypes);
            });
    }

    // the segment types of all calculations have to be rebuilt from the reloaded sigma profiles,
    // also if the reload failed for some molecules
    for (calculation& _calculation : calculations) {
        _calculation.sigmaProfilesOutdated = true;
    }

    e.rethrow();
}

// sets area and volume of monoatomic cations from the given radii, needed once after loading the molecules
// or whenever the radii change (sw_alwaysCalculateSizeRelatedParameters == 1)
void resizeMonoatomicCations(parameters& param, std::vector<std::shared_ptr<molecule>> molecules);

void prepareMoleculesForCalculation(parameters& param) {
    if (param.sw_alwaysReloadSigmaProfiles == 1 && sigmaProfilesUsedInCalculation) {
        reloadAllMolecules();
    }
    if (param.sw_alwaysCalculateSizeRelatedParameters == 1 || !monoatomicCationsResized) {
        resizeMonoatomicCations(param, molecules);
        monoatomicCationsResized = true;
    }
}

void resizeMonoatomicCations(parameters& param, std::vector<std::shared_ptr<molecule>> molecules) {
#ifdef MEASURE_TIME
    std::chrono::high_resolution_clock::time_point rescaleSegments_last = std::chrono::high_resolution_clock::now();
#endif
    // scale A and V for monoatomic cations
    for (int i = 0; i < molecules.size(); i++) {

        if (molecules[i]->moleculeGroup == 3) {
            int AN = molecules[i]->atomAtomicNumbers(0);
            double R_i = param.R_i[AN];
            if (!(R_i > 0)) {
                throw std::runtime_error("No radius was given for the monoatomic cation '" + molecules[i]->name + "' (atomic number " + std::to_string(AN) +
                    "). With sw_SR_differentiateMoleculeGroups = 1 the radius of every monoatomic cation has to be given in parameters['radii'].");
            }
            molecules[i]->Area = (4 * PI * R_i * R_i);
            molecules[i]->Volume = (4.0 / 3 * PI * R_i * R_i * R_i);
        }
    }
#ifdef MEASURE_TIME
    rescaleSegments_total_ms += (const unsigned long)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - rescaleSegments_last).count());
#endif
}

void rescaleSegments(parameters& param, calculation& _calculation) {
    // rescale segments for monoatomic cations only if present
    if (_calculation.segments.numberOfSegmentsForGroup[3] != 0) {

        // monoatomic cation
        std::unordered_map<unsigned short, std::array<int, 2>> segmentIndicesBelongingToASpecificAtomicNumber;
        for (int j = _calculation.segments.lowerBoundIndexForGroup[3];
            j < _calculation.segments.upperBoundIndexForGroup[3]; j++) {

            unsigned short AN = _calculation.segments.SegmentTypeAtomicNumber[j];

            auto found = segmentIndicesBelongingToASpecificAtomicNumber.find(AN);
            if (found == segmentIndicesBelongingToASpecificAtomicNumber.end()) {
                segmentIndicesBelongingToASpecificAtomicNumber[AN] = { j, -1 };
            }
            else if (found->second[1] == -1) {
                found->second[1] = j;
            }
            else {
                throw std::runtime_error("More than 2 segment types were found for the monoatomic cation with atomic number " + std::to_string(AN) + ".");
            }

        }
        for (auto it = segmentIndicesBelongingToASpecificAtomicNumber.begin(); it != segmentIndicesBelongingToASpecificAtomicNumber.end(); ++it) {

            int AN = it->first;
            int ind_left = it->second[0];
            int ind_right = it->second[1];

            if (ind_right == -1) {
                throw std::runtime_error("Only one segment type was found for the monoatomic cation with atomic number " + std::to_string(AN) +
                    " (its screening charge density lies exactly on a charge raster point), rescaling its segments is not possible.");
            }

            if (ind_right - ind_left > 1) {
                throw std::runtime_error("More than 2 segment types were found for the monoatomic cation with atomic number " + std::to_string(AN) + ".");
            }

            int ind_molecule = -1;

            for (int k = 0; k < int(_calculation.segments.SegmentTypeAreas[ind_left].size()); k++) {
                if (_calculation.segments.SegmentTypeAreas[ind_left][k] > 0.0) {
                    ind_molecule = k;
                    break;
                }
            }

            if (ind_molecule == -1) {
                throw std::runtime_error("The molecule belonging to the monoatomic cation with atomic number " + std::to_string(AN) + " was not found.");
            }

            if (!(param.R_i[AN] > 0)) {
                throw std::runtime_error("No radius was given for the monoatomic cation with atomic number " + std::to_string(AN) + " in parameters['radii'].");
            }

            int screeningCharge = int(_calculation.components[ind_molecule]->moleculeCharge * -1);

            double newArea = (4 * PI * param.R_i[AN] * param.R_i[AN]);
            double newSigma = screeningCharge / newArea;

            int ind_Sigma_left = getLeftChargeRasterIndex(param, newSigma, "rescaled screening charge density (sigma)", _calculation.components[ind_molecule]->name);

            double sigmaLeft = param.ChargeRaster[ind_Sigma_left];
            double sigmaRight = param.ChargeRaster[ind_Sigma_left + 1];

            double AsigmaRight = newArea * (newSigma - sigmaLeft) / param.sigmaStep;
            double AsigmaLeft = newArea * (sigmaRight - newSigma) / param.sigmaStep;

            _calculation.segments.SegmentTypeAreas[ind_left][ind_molecule] = AsigmaLeft;
            _calculation.segments.SegmentTypeSigma[ind_left] = sigmaLeft;

            _calculation.segments.SegmentTypeAreas[ind_right][ind_molecule] = AsigmaRight;
            _calculation.segments.SegmentTypeSigma[ind_right] = sigmaRight;

        }
    }
}

void calculateSegmentConcentrations(calculation& _calculation) {
    // calculate the mole fraction of segments for each concentration
    for (int j = 0; j < _calculation.concentrations.size(); j++) {

        int firstNonZeroSegmentIndex = 0;
        int lastNonZeroSegmentIndex = int(_calculation.segments.size()) - 1;

        std::vector<double> segmentConcentration(_calculation.segments.size(), 0.0);

        double sumAreaSegmentsConcentrationj = 0.0;

        for (int k = 0; k < _calculation.segments.size(); k++) {

            double areaSegmentK = 0.0;

            for (int m = 0; m < _calculation.components.size(); m++) {
                double thisArea = _calculation.concentrations[j][m] * _calculation.segments.SegmentTypeAreas[k][m];
                areaSegmentK += thisArea;
                sumAreaSegmentsConcentrationj += thisArea;
            }

            segmentConcentration[k] = areaSegmentK;
        }

        double cumulativeSumOfAreasFromSegmentZeroOn = 0.0;

        bool firstNonZeroSegmentIndex_found = false;
        for (int k = 0; k < _calculation.segments.size(); k++) {

            cumulativeSumOfAreasFromSegmentZeroOn += segmentConcentration[k];

            if (cumulativeSumOfAreasFromSegmentZeroOn != 0 && !firstNonZeroSegmentIndex_found) {
                firstNonZeroSegmentIndex = k;
                firstNonZeroSegmentIndex_found = true;
            }
            if (segmentConcentration[k] != 0) {
                lastNonZeroSegmentIndex = k;
            }

            _calculation.segmentConcentrations(k, j) = calcType(segmentConcentration[k] / sumAreaSegmentsConcentrationj);
        }

        _calculation.lowerBoundIndexForCOSMOSPACECalculation[j] = RoundDownToNextMultipleOfEight(firstNonZeroSegmentIndex);
        _calculation.upperBoundIndexForCOSMOSPACECalculation[j] = RoundUpToNextMultipleOfEight(lastNonZeroSegmentIndex + 1);
    }
}

// checks all concentrations (including the ones of the reference states) and temperatures of a calculation
void checkConcentrationsAndTemperatures(const calculation& _calculation) {

    const std::string prefix = "For calculation number " + std::to_string(_calculation.number) + " the mixture number ";

    for (int j = 0; j < _calculation.concentrations.size(); j++) {

        const double temperature = _calculation.temperatures[j];
        if (!std::isfinite(temperature) || temperature <= 0) {
            throw std::runtime_error(prefix + std::to_string(j) + " has the invalid temperature " + std::to_string(temperature) + " K.");
        }

        double sumOfConcentrations = 0;
        double mix_chrg = 0;
        for (int k = 0; k < _calculation.components.size(); k++) {
            const double x = _calculation.concentrations[j][k];
            if (!std::isfinite(x) || x < 0) {
                throw std::runtime_error(prefix + std::to_string(j) + " has the invalid mole fraction " + std::to_string(x) + " for component " + std::to_string(k) + ".");
            }
            sumOfConcentrations += x;
            mix_chrg = mix_chrg + _calculation.components[k]->moleculeCharge * x;
        }

        if (std::abs(1.0 - sumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
            throw std::runtime_error(prefix + std::to_string(j) + ": the mole fractions do not add up to unity. residual concentration: " + std::to_string(std::abs(1.0 - sumOfConcentrations)));
        }

        if (std::abs(mix_chrg) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
            throw std::runtime_error(prefix + std::to_string(j) + " is not electroneutral. residual charge: " + std::to_string(std::abs(mix_chrg)));
        }
    }
}

void finishCalculationInitiation(calculation& _calculation) {

    if (_calculation.concentrations.size() > 65535) {
        throw std::runtime_error("Too many calculations, other datatype would be necessary for newCalculation.referenceStates to cope with this amount. (unsigned short used so far allowing for up to 65535)");
    }

    checkConcentrationsAndTemperatures(_calculation);

    // save sorting of concentrations according to conditions, this clusters the calculation of the interaction matrix
    // sorting by similar concentration does not make sense as the gammas from the previous iteration
    // are saved for every concentration already accelerating the COSMOSPACE convergence
    std::vector<int> sortingVector(_calculation.concentrations.size());
    std::iota(sortingVector.begin(), sortingVector.end(), 0);
    std::sort(sortingVector.begin(), sortingVector.end(),
        [&](int i, int j) {
            if (_calculation.temperatures[i] != _calculation.temperatures[j])
                return _calculation.temperatures[i] < _calculation.temperatures[j];
            else
                return i < j;
        });

    apply_vector_permutation_in_place(_calculation.concentrations, sortingVector);
    apply_vector_permutation_in_place(_calculation.temperatures, sortingVector);

    _calculation.actualConcentrationIndices = std::vector<int>(_calculation.concentrations.size());
    for (int j = 0; j < _calculation.concentrations.size(); j++) {
        _calculation.actualConcentrationIndices[sortingVector[j]] = j;
    }

    for (int j = 0; j < _calculation.concentrations.size(); j++) {
        int TauIndex = _calculation.addOrFindTauIndexForConditions(_calculation.temperatures[j]);
        _calculation.TauConcentrationIndices[TauIndex].push_back(j);
    }

    // initiate arrays for the calculation
    _calculation.segmentGammas = MatrixCalcType::Constant(RoundUpToNextMultipleOfEight(int(_calculation.segments.size())), int(_calculation.concentrations.size()), 1.0);
    _calculation.segmentConcentrations = MatrixCalcType::Zero(RoundUpToNextMultipleOfEight(int(_calculation.segments.size())), int(_calculation.concentrations.size()));

    _calculation.PhiDash_pxi = Eigen::MatrixXd::Zero(_calculation.concentrations.size(), _calculation.components.size());
    _calculation.ThetaDash_pxi = Eigen::MatrixXd::Zero(_calculation.concentrations.size(), _calculation.components.size());

    for (int j = 0; j < _calculation.concentrations.size(); j++) {
        _calculation.lowerBoundIndexForCOSMOSPACECalculation.push_back(0);
        _calculation.upperBoundIndexForCOSMOSPACECalculation.push_back(int(_calculation.segments.size()));
    }

    _calculation.shrink_to_fit();
}

void calculateLnGammaCombinatorial(parameters& param, calculation& _calculation) {

    std::vector<double> averageVolumes(_calculation.concentrations.size(), 0.0);
    std::vector<double> averageAreas(_calculation.concentrations.size(), 0.0);

    // Volume fraction to mole fraction ratio for all compositions
    // Area fraction to mole fraction ratio for all compositions
    for (int j = 0; j < _calculation.concentrations.size(); j++) {
        for (int k = 0; k < _calculation.components.size(); k++) {
            averageVolumes[j] += _calculation.components[k]->Volume * _calculation.concentrations[j][k];
            averageAreas[j] += _calculation.components[k]->Area * _calculation.concentrations[j][k];
        }

        for (int k = 0; k < _calculation.components.size(); k++) {
            _calculation.PhiDash_pxi(j, k) = _calculation.components[k]->Volume / averageVolumes[j];
            _calculation.ThetaDash_pxi(j, k) = _calculation.components[k]->Area / averageAreas[j];
        }
    }

    std::vector<double> qi_std(_calculation.components.size());

    /* Calculation of molecule specific parameters for combinatorial contribution */
    for (int i = 0; i < _calculation.components.size(); i++) {
        qi_std[i] = _calculation.components[i]->Area / param.comb_SG_A_std;
    }

    Eigen::MatrixXd lnGamaForCalculations = Eigen::MatrixXd::Zero(_calculation.concentrations.size(), _calculation.components.size());

    if (param.sw_combTerm == 1) { // tested
        /* Staverman-Guggenheim term
           cp. Kikic (1980) or Lin & Sandler (2002)
           Remark: square brackets in L&S paper for the li term are wrong as can be seen in monograph Prausnitz Lichtenthaler */
        for (int i = 0; i < _calculation.concentrations.size(); i++) {

            double buffdb1 = 0;
            for (int j = 0; j < _calculation.components.size(); j++) {
                buffdb1 = _calculation.PhiDash_pxi(i, j) / _calculation.ThetaDash_pxi(i, j);
                lnGamaForCalculations(i, j) = log(_calculation.PhiDash_pxi(i, j)) + 1 - _calculation.PhiDash_pxi(i, j) - \
                    param.comb_SG_z_coord * 0.5 * qi_std[j] * (log(buffdb1) + 1 - buffdb1);
            }
        }
    }
    else if (param.sw_combTerm == 2) { // tested
        /* Klamt (2003) */
        for (int i = 0; i < _calculation.concentrations.size(); i++) {
            double buffdb1 = 0;
            for (int j = 0; j < _calculation.components.size(); j++) {
                buffdb1 = buffdb1 + _calculation.concentrations[i][j] * log(_calculation.components[j]->Volume);
            }

            for (int j = 0; j < _calculation.components.size(); j++) {
                lnGamaForCalculations(i, j) = (param.comb_lambda0 * (log(_calculation.components[j]->Volume) - buffdb1)) \
                    - (param.comb_lambda1 * (_calculation.PhiDash_pxi(i, j) - 1)) \
                    - (param.comb_lambda2 * (_calculation.ThetaDash_pxi(i, j) - 1));
            }
        }
    }
    else if (param.sw_combTerm == 3) { // tested
        /*  mod. Staverman-Guggenheim with exponential scaling, compare e.g. Soares (2011), Kikic (1980), Donohue Prausnitz (1975) */
        for (int i = 0; i < _calculation.concentrations.size(); i++) {

            double buffdb1 = 0;
            for (int j = 0; j < _calculation.components.size(); j++) {
                buffdb1 += pow(_calculation.components[j]->Volume, param.comb_modSG_exp) * _calculation.concentrations[i][j];
            }

            double buffdb2 = 0;
            for (int j = 0; j < _calculation.components.size(); j++) {
                buffdb2 = _calculation.PhiDash_pxi(i, j) / _calculation.ThetaDash_pxi(i, j);

                double PhiDash_pxi_mod_i = pow(_calculation.components[j]->Volume, param.comb_modSG_exp) / buffdb1;

                lnGamaForCalculations(i, j) = log(PhiDash_pxi_mod_i) + 1 - PhiDash_pxi_mod_i - \
                    param.comb_SG_z_coord * 0.5 * qi_std[j] * (log(buffdb2) + 1 - buffdb2);
            }
        }
    }
    else if (param.sw_combTerm == 4) {
        /*  mod. Staverman-Guggenheim by Grensemann & Gmehling, Ind. Eng. Chem. Res. 44 (2005) 1610, Eqs. (7) to (15),
            especially developed for COSMO-RS.
            Deviation from the paper: r_i# and q_i# (Eqs. 12 and 15) are calculated as the mole fraction averaged volume
            and area of all other molecules, sum_{j!=i} x_j r_j / (1 - x_i). This is identical to the paper for binary
            mixtures, but stays finite for mole fractions of zero and is an average for more than two components.
            For x_i = 1 the values do not matter (V_i = F_i = 1), r_i and q_i are used. */
        const int numberOfComponents = int(_calculation.components.size());

        for (int i = 0; i < _calculation.concentrations.size(); i++) {

            const std::vector<double>& x = _calculation.concentrations[i];

            double sum_q_times_x = 0;
            double sum_r_times_x = 0;

            for (int j = 0; j < numberOfComponents; j++) {
                sum_q_times_x += _calculation.components[j]->Area * x[j];
                sum_r_times_x += _calculation.components[j]->Volume * x[j];
            }

            for (int j = 0; j < numberOfComponents; j++) {

                const double q_j = _calculation.components[j]->Area;
                const double r_j = _calculation.components[j]->Volume;

                // mole fraction averaged area and volume of the other molecules
                double q_hash = q_j;
                double r_hash = r_j;
                const double sumOfOtherMoleFractions = 1.0 - x[j];
                if (sumOfOtherMoleFractions > 0) {
                    double sum_q_others = 0;
                    double sum_r_others = 0;
                    for (int k = 0; k < numberOfComponents; k++) {
                        if (k == j) continue;
                        sum_q_others += _calculation.components[k]->Area * x[k];
                        sum_r_others += _calculation.components[k]->Volume * x[k];
                    }
                    q_hash = sum_q_others / sumOfOtherMoleFractions;
                    r_hash = sum_r_others / sumOfOtherMoleFractions;
                }

                const double q_min = std::min(q_j, q_hash);
                const double q_max = std::max(q_j, q_hash);

                const double r_min = std::min(r_j, r_hash);
                const double r_max = std::max(r_j, r_hash);

                const double F = pow(q_j / sum_q_times_x, param.comb_SGG_lambda * (1 - (q_min / q_max)));
                const double V = pow(r_j / sum_r_times_x, param.comb_SGG_beta * (1 - (r_min / r_max)));

                // Eq. (7)
                lnGamaForCalculations(i, j) = 1 - V + log(V) + 1 - (V / F) + log(V / F);
            }
        }
    }
    else if (param.sw_combTerm == 5) { // not tested yet
        /*  Franke & Hannebauer (2011) */
        for (int i = 0; i < _calculation.concentrations.size(); i++) {
            for (int j = 0; j < _calculation.components.size(); j++) {
                lnGamaForCalculations(i, j) = param.comb_lambda0 * log(_calculation.components[j]->Volume) \
                    + param.comb_lambda1 * (1 - (_calculation.components[j]->Volume / averageVolumes[i]) - log(averageVolumes[i])) \
                    + param.comb_lambda2 * (1 - (_calculation.components[j]->Area / averageAreas[i]) - log(averageAreas[i]));
            }
        }
    }
    else if (param.sw_combTerm != 0) { 
        throw std::runtime_error("Error: Invalid switch value for combinatorial term.\n");
    }


    /* Convert activity coefficients to correct reference states. */
    for (int h = 0; h < _calculation.originalNumberOfCalculations; h++) {

        int i = _calculation.actualConcentrationIndices[h];

        for (int j = 0; j < _calculation.components.size(); j++) {

            double multiplier = 1.0;
            // for an ion with reference state PureComponentsOnlyNeutral
            if (_calculation.referenceStateType[h] == 1 && _calculation.referenceStateCalculationIndices[h][j] == -1) {
                multiplier = NAN;
            }
            _calculation.lnGammaCombinatorial(h, j) = multiplier * lnGamaForCalculations(i, j);
#ifdef PRINT_DEBUG_INFO
            display("lnGammaCombinatorial_1_" + std::to_string(h) + ": " + std::to_string(_calculation.lnGammaCombinatorial(h, j)) + "\n");
#endif

            // only subtract reference state value if reference state type is not COSMO
            // and a reference state exists (not the case for ions with reference state PureComponentsOnlyNeutral)
            if (_calculation.referenceStateType[h] != 3 && _calculation.referenceStateType[h] != 4 && _calculation.referenceStateCalculationIndices[h][j] != -1) {
                int indexOfReferenceStateCalculation = _calculation.actualConcentrationIndices[_calculation.referenceStateCalculationIndices[h][j]];

                _calculation.lnGammaCombinatorial(h, j) -= lnGamaForCalculations(indexOfReferenceStateCalculation, j);
            }
        }
    }
}

// the COSMOSPACE equations are solved by successive substitution
// instead of calculating the complete matrix-vector product and then calculating the gammas for the next iteration
// each row is calculated and the gamma for the row just calculated is already used in the calculation of the next row.
// This was done firstly by accident, but it was found, that it accelerated the convergence by a factor of at least 4.
void calculateLnGammaResidual(parameters& param, calculation& _calculation) {

    const int numberOfSegments = int(_calculation.segments.size());
    const int nMultipleOfEight = RoundUpToNextMultipleOfEight(numberOfSegments);

    Eigen::MatrixXd temporary_lnGammaMolecule = Eigen::MatrixXd::Zero(_calculation.concentrations.size(), _calculation.components.size());

    // temporary arrays are indexed by the internal (sorted) concentration index i,
    // the output arrays by the original concentration index h
    Eigen::Tensor<double, 3, Eigen::RowMajor> temporary_contactStatistics;
    Eigen::Tensor<double, 4, Eigen::RowMajor> temporary_averageInteractionEnergies;
    Eigen::Tensor<double, 3, Eigen::RowMajor> temporary_partialMolarEnergies;

    if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {

        temporary_contactStatistics = Eigen::Tensor<double, 3, Eigen::RowMajor>(int(_calculation.concentrations.size()),
            int(_calculation.components.size()),
            int(_calculation.components.size()));

        temporary_contactStatistics.setZero();

        temporary_averageInteractionEnergies = Eigen::Tensor<double, 4, Eigen::RowMajor>(int(_calculation.concentrations.size()),
            param.numberOfPartialInteractionMatrices + 1, // +1 because A_int is the first one
            int(_calculation.components.size()),
            int(_calculation.components.size()));

        temporary_averageInteractionEnergies.setZero();

        if (param.sw_calculateContactStatisticsAndAdditionalProperties == 2) {
            temporary_partialMolarEnergies = Eigen::Tensor<double, 3, Eigen::RowMajor>(int(_calculation.concentrations.size()),
                param.numberOfPartialInteractionMatrices + 1, // +1 because A_int is the first one
                int(_calculation.components.size()));

            temporary_partialMolarEnergies.setZero();

        }

    }


    const double dampingFactor = 0.6;
    const double dampingFactorComplement = 1 - dampingFactor;
    const double convergenceThreshhold = 0.00001;   // originally this was 0.0000001;
                                                    // however when using AVX and vectorized code (especially with floats)
                                                    // the accuracy of the calculation drops
                                                    // without increasing this value the gammas do not converge
                                                    // however since the convergence criteria was changed from
                                                    // gamma to log(gamma) the larger convergence threshhold 
                                                    // still leads to numerically very similar results

    const int maximumNumberOfIterations = 50000;
    double newGamma = 0.0;


    MatrixCalcType A_int = MatrixCalcType::Zero(numberOfSegments, numberOfSegments);
    MatrixCalcType Tau = MatrixCalcType::Zero(nMultipleOfEight, numberOfSegments);
    MatrixCalcType TauX = MatrixCalcType::Zero(nMultipleOfEight, numberOfSegments);

    std::vector<Eigen::MatrixXd> partialInteractionMatrices;

    if (param.numberOfPartialInteractionMatrices > 0) {
        for (int h = 0; h < param.numberOfPartialInteractionMatrices; h++) {
            partialInteractionMatrices.push_back(Eigen::MatrixXd::Zero(numberOfSegments, numberOfSegments));
        }
    }

    for (int g = 0; g < _calculation.TauConcentrationIndices.size(); g++) {

        // conditions
        double temperature = _calculation.TauTemperatures[g];

#ifdef MEASURE_TIME
        std::chrono::high_resolution_clock::time_point calculateTau_last = std::chrono::high_resolution_clock::now();
#endif
        // calculate interaction matrices
        if (param.numberOfPartialInteractionMatrices > 0) {
            for (int h = 0; h < param.numberOfPartialInteractionMatrices; h++) {
                partialInteractionMatrices[h].setZero();
            }
        }

        // reset, as not all entries are rewritten (e.g. ion segments) and the segment reference state
        // would otherwise be applied again on the values of the previous temperature
        A_int.setZero();
        calculateInteractionMatrix(_calculation.segments, A_int, partialInteractionMatrices, param, temperature);

        // if the full interaction matrix is needed, fill upper right half
        if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0){
            for (int i = 0; i < numberOfSegments; i++) {
                for (int j = i + 1; j < numberOfSegments; j++) {
                    // always j >= i + 1
                    A_int(i, j) = A_int(j, i);
                }
            }
            for (int h = 0; h < param.numberOfPartialInteractionMatrices; h++) {
                for (int i = 0; i < numberOfSegments; i++) {
                    for (int j = i + 1; j < numberOfSegments; j++) {
                        // always j >= i + 1
                        partialInteractionMatrices[h](i, j) = partialInteractionMatrices[h](j, i);
                    }
                }
            }
        }

        // calculate Tau from the interaction matrix A_int
        int idx = -1;
        int columnSum = 0;
        const calcType minus_div_RT = calcType(-1.0 / (R_GAS_CONSTANT * temperature));

        calcType val;
        calcType* Tau_1D = &(Tau(0, 0));
        for (int i = 0; i < numberOfSegments; i++) {
            columnSum = i * nMultipleOfEight;
            for (int j = i; j < numberOfSegments; j++) {
                // always j >= i
                val = exp(A_int(j, i) * minus_div_RT);

                idx = columnSum + j;
                Tau_1D[idx] = val;
                Tau_1D[j * nMultipleOfEight + i] = val;
            }
        }

#ifdef MEASURE_TIME
        calculateTau_total_ms += (const unsigned long)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateTau_last).count());
#endif

#ifdef PRINT_DEBUG_INFO
        Eigen::MatrixXd Tau_d = Tau.cast<double>();
        WriteEigenMatrixToFile("Tau_" + std::to_string(_calculation.number) + "_" + std::to_string(g), Tau_d);
#endif

        for (int h = 0; h < _calculation.TauConcentrationIndices[g].size(); h++) {

            int i = _calculation.TauConcentrationIndices[g][h];

            calcType* gammas = &(_calculation.segmentGammas(0, i));
#ifdef MEASURE_TIME
            std::chrono::high_resolution_clock::time_point calculateCOSMOSPACE_last = std::chrono::high_resolution_clock::now();
#endif

            int lowerBoundIndexForCOSMOSPACECalculation = _calculation.lowerBoundIndexForCOSMOSPACECalculation[i];
            int upperBoundIndexForCOSMOSPACECalculation = _calculation.upperBoundIndexForCOSMOSPACECalculation[i];

            { // as local as possible scope for all variables used in your loops

                const calcType* vTau_1D = &(Tau(0, 0));
                const calcType* vX = &(_calculation.segmentConcentrations(0, i));
                calcType* vTauX_1D = &(TauX(0, 0));
                int idx = 0;
                int columnSum = 0;
                using SimdCalc = Simd<calcType>;

                for (int j = 0; j < numberOfSegments; j++) {

                    columnSum = j * nMultipleOfEight;

                    for (int k = lowerBoundIndexForCOSMOSPACECalculation; k < upperBoundIndexForCOSMOSPACECalculation; k += SimdCalc::lanes) {
                        idx = columnSum + k;
                        SimdCalc::store(
                            vTauX_1D + idx,
                            SimdCalc::mul(SimdCalc::load(vTau_1D + idx), SimdCalc::load(vX + k))
                        );
                    }

                }
            }

#ifdef PRINT_DEBUG_INFO
            Write1DArrayToFile<calcType>("X_" + std::to_string(_calculation.number) + "_" + std::to_string(g) + "_" + std::to_string(i), &(_calculation.segmentConcentrations(0, i)), 1, nMultipleOfEight);
            Eigen::MatrixXd TauX_d = TauX.cast<double>();
            WriteEigenMatrixToFile("TauX_" + std::to_string(_calculation.number) + "_" + std::to_string(g) + "_" + std::to_string(i), TauX_d);
#endif

            int numberOfConvergedGammas = 0;
            int numberOfIteration = 0;

            while (numberOfConvergedGammas != numberOfSegments) {

                numberOfIteration++;

                if (numberOfIteration > maximumNumberOfIterations) {

                    // no output here: this runs in OpenMP threads (in Python without the GIL),
                    // all information is passed on in the exception message instead
                    std::string information = "";

                    if (param.sw_skip_COSMOSPACE_errors == 0) {

                        information += "\nSYSTEM_COMPONENTS:\n";
                        for (int j = 0; j < _calculation.components.size(); j++) {
                            information += " -" + _calculation.components[j]->name + "\n";
                        }

                        calcType* Tau_1D = &(Tau(0, 0));
                        for (int j = 0; j < numberOfSegments * nMultipleOfEight; j++) {
                            if (std::isnan(Tau_1D[j])) {
                                information += " Some Tau entries are NaN.";
                                break;
                            }
                        }

                        for (int j = 0; j < numberOfSegments * nMultipleOfEight; j++) {
                            if (std::isinf(Tau_1D[j])) {
                                information += " Some Tau entries are inf.";
                                break;
                            }
                        }

                        for (int j = 0; j < numberOfSegments; j++) {
                            if (std::isnan(_calculation.segmentGammas(j, i))) {
                                information += " Some gammas are NaN.";
                                break;
                            }
                        }
                        for (int j = 0; j < numberOfSegments; j++) {
                            if (std::isinf(_calculation.segmentGammas(j, i))) {
                                information += " Some gammas are inf.";
                                break;
                            }
                        }
                    }
                    else {
                        // set all segment gammas to 1 as initial point for next execution
                        _calculation.segmentGammas = MatrixCalcType::Constant(RoundUpToNextMultipleOfEight(int(_calculation.segments.size())), int(_calculation.concentrations.size()), 1.0f);
                    }

                    throw std::runtime_error("COSMOSPACE did not converge for calculation " + std::to_string(_calculation.number) + " on concentration " + std::to_string(i) + ": maximum number of iterations reached. (" + std::to_string(numberOfConvergedGammas) + "/" + std::to_string(numberOfSegments) + ")" + information);
                }

                numberOfConvergedGammas = 0;

                int columnSum = 0;

                for (int j = 0; j < numberOfSegments; j++) {
                    newGamma = 0.0;

                    columnSum = j * nMultipleOfEight;

                    { // as local as possible scope for all variables used in your loops

                        const calcType* vTauX = &(TauX(0, 0));
                        calcType* vGammas = gammas;
                        using SimdCalc = Simd<calcType>;

                        SimdCalc::Vec tempVectorSum = SimdCalc::zero();

                        for (int k = lowerBoundIndexForCOSMOSPACECalculation; k < upperBoundIndexForCOSMOSPACECalculation; k += SimdCalc::lanes) {
                            tempVectorSum = SimdCalc::fmadd(
                                SimdCalc::load(vTauX + columnSum + k),
                                SimdCalc::load(vGammas + k),
                                tempVectorSum
                            );
                        }

                        newGamma = double(SimdCalc::horizontal_sum(tempVectorSum));
                }

                    newGamma = 1 / newGamma;
                    double oldGamma = double(gammas[j]);

                    // check convergence criteria on the undamped update, otherwise
                    // damping would loosen the criterion by the factor 1 / (1 - dampingFactor)
                    if (std::abs(log(newGamma) - log(oldGamma)) <= convergenceThreshhold) {
                        numberOfConvergedGammas++;
                    }

                    // apply damping
                    if (numberOfIteration > 200) {
                        newGamma = dampingFactor * newGamma + dampingFactorComplement * oldGamma;
                    }

                    gammas[j] = calcType(newGamma);
                }
            }

#ifdef PRINT_DEBUG_INFO
            display("COSMOSPACE converged with this number of iterations: " + std::to_string(numberOfIteration) + "\n");
            Write1DArrayToFile<calcType>("Gammas_" + std::to_string(_calculation.number) + "_" + std::to_string(g) + "_" + std::to_string(i), gammas, 1, numberOfSegments);
#endif

#ifdef MEASURE_TIME
            calculateCOSMOSPACE_total_ms += (const unsigned long)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateCOSMOSPACE_last).count());
            std::chrono::high_resolution_clock::time_point calculateContactStatistics_last = std::chrono::high_resolution_clock::now();
#endif


            if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {
                calculateContactStatistics(_calculation, A_int, partialInteractionMatrices, Tau, gammas, int(i), temporary_contactStatistics, temporary_averageInteractionEnergies, temporary_partialMolarEnergies, param);
            }

#ifdef MEASURE_TIME
            calculateContactStatistics_total_ms += (const unsigned long)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateContactStatistics_last).count());
            std::chrono::high_resolution_clock::time_point calculateGammasForMolecules_last = std::chrono::high_resolution_clock::now();
#endif
            double div_Aeff = 1 / param.Aeff;

            for (int j = 0; j < _calculation.components.size(); j++) {

                double lnGammaMolecule = 0;

                for (int k = 0; k < numberOfSegments; k++) {
                    lnGammaMolecule += _calculation.segments.SegmentTypeAreas[k][j] * div_Aeff * log(double(gammas[k])); // log(gammas[k]) is calculated more than once although not necessary
                }
#ifdef PRINT_DEBUG_INFO
                display("lnGammaResidual_1_" + std::to_string(h) + ": " + std::to_string(lnGammaMolecule) + "\n");
#endif
                
                temporary_lnGammaMolecule(i, j) = lnGammaMolecule;
            }

#ifdef MEASURE_TIME
            calculateGammasForMolecules_total_ms += (const unsigned long)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateGammasForMolecules_last).count());
#endif
        }
    }

#ifdef MEASURE_TIME
    std::chrono::high_resolution_clock::time_point calculateGammasForMolecules_last = std::chrono::high_resolution_clock::now();
#endif

    /* Convert to correct reference states. */
    for (int h = 0; h < _calculation.originalNumberOfCalculations; h++) {

        int i = _calculation.actualConcentrationIndices[h];

        for (int j = 0; j < _calculation.components.size(); j++) {


            double multiplier = 1.0;
            // for an ion with reference state PureComponentsOnlyNeutral
            if (_calculation.referenceStateType[h] == 1 && _calculation.referenceStateCalculationIndices[h][j] == -1) {
                multiplier = NAN;
            }
            _calculation.lnGammaResidual(h, j) = multiplier * temporary_lnGammaMolecule(i, j);

#ifdef PRINT_DEBUG_INFO
            display("lnGammaResidual_2_" + std::to_string(h) + ": " + std::to_string(_calculation.lnGammaResidual(h, j)) + "\n");
#endif

            if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {
                for (int k = 0; k < _calculation.components.size(); k++) {
                    // contact statistics do not depend on the reference state, hence they are copied as is
                    _calculation.contactStatistics(h, j, k) = temporary_contactStatistics(i, j, k);

                    _calculation.averageSurfaceEnergies(h, 0, j, k) = multiplier * temporary_averageInteractionEnergies(i, 0, j, k);

                    for (int l = 0; l < param.numberOfPartialInteractionMatrices; l++) {
                        _calculation.averageSurfaceEnergies(h, l + 1, j, k) = multiplier * temporary_averageInteractionEnergies(i, l + 1, j, k);
                    }
                }

                if (param.sw_calculateContactStatisticsAndAdditionalProperties == 2) {
                    _calculation.partialMolarEnergies(h, 0, j) = multiplier * temporary_partialMolarEnergies(i, 0, j);

                    for (int l = 0; l < param.numberOfPartialInteractionMatrices; l++) {
                        _calculation.partialMolarEnergies(h, l + 1, j) = multiplier * temporary_partialMolarEnergies(i, l + 1, j);
                    }
                }
            }

            // only subtract reference state value if reference state type is not COSMO
            // and a reference state exists (not the case for ions with reference state PureComponentsOnlyNeutral)
            if (_calculation.referenceStateType[h] != 3 && _calculation.referenceStateType[h] != 4 && _calculation.referenceStateCalculationIndices[h][j] != -1) {
                int indexOfReferenceStateCalculation = _calculation.actualConcentrationIndices[_calculation.referenceStateCalculationIndices[h][j]];

                _calculation.lnGammaResidual(h, j) -= temporary_lnGammaMolecule(indexOfReferenceStateCalculation, j);

                if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {
                    for (int k = 0; k < _calculation.components.size(); k++) {
                        _calculation.averageSurfaceEnergies(h, 0, j, k) -= temporary_averageInteractionEnergies(indexOfReferenceStateCalculation, 0, j, k);

                        for (int l = 0; l < param.numberOfPartialInteractionMatrices; l++) {
                            _calculation.averageSurfaceEnergies(h, l + 1, j, k) -= temporary_averageInteractionEnergies(indexOfReferenceStateCalculation, l + 1, j, k);
                        }
                    }

                    if (param.sw_calculateContactStatisticsAndAdditionalProperties == 2) {
                        _calculation.partialMolarEnergies(h, 0, j) -= temporary_partialMolarEnergies(indexOfReferenceStateCalculation, 0, j);

                        for (int l = 0; l < param.numberOfPartialInteractionMatrices; l++) {
                            _calculation.partialMolarEnergies(h, l + 1, j) -= temporary_partialMolarEnergies(indexOfReferenceStateCalculation, l + 1, j);
                        }
                    }
                }
            }
        }
    }

#ifdef MEASURE_TIME
    calculateGammasForMolecules_total_ms += (const unsigned long)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateGammasForMolecules_last).count());
#endif
}

void calculate(std::vector<int>& calculationIndices) {

    // this is needed to catch exceptions in the OPENMP threads and rethrow them after the parallel section ends
    threadException e;
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (int i = 0; i < calculationIndices.size(); i++) {

        e.run([=] {

            int calculationIndex = calculationIndices[i];

#ifdef MEASURE_TIME
            std::chrono::high_resolution_clock::time_point rescaleSegments_last = std::chrono::high_resolution_clock::now();
#endif
            if (calculations[calculationIndex].sigmaProfilesOutdated) {

                // the segment types changed, all size related properties have to be calculated again
                calculations[calculationIndex].sizeRelatedPropertiesInitialized = false;
                calculations[calculationIndex].segments.clear();

                for (int j = 0; j < calculations[calculationIndex].components.size(); j++) {

                    std::shared_ptr<molecule> thisMolecule = calculations[calculationIndex].components[j];

                    for (int k = 0; k < thisMolecule->segments.size(); k++) {
                        calculations[calculationIndex].segments.add((unsigned short)j, thisMolecule->segments.SegmentTypeGroup[k],
                            thisMolecule->segments.SegmentTypeSigma[k],
                            thisMolecule->segments.SegmentTypeSigmaCorr[k],
                            thisMolecule->segments.SegmentTypeHBtype[k],
                            thisMolecule->segments.SegmentTypeAtomicNumber[k],
                            thisMolecule->segments.SegmentTypeAtomicPolariz[k],
                            thisMolecule->segments.SegmentTypeAreas[k][0]);
                    }
                }
                calculations[calculationIndex].segments.sort();
                calculations[calculationIndex].segmentGammas =MatrixCalcType::Constant(RoundUpToNextMultipleOfEight(int(calculations[calculationIndex].segments.size())), int(calculations[calculationIndex].concentrations.size()), 1.0f);
                calculations[calculationIndex].segmentConcentrations =MatrixCalcType::Zero(RoundUpToNextMultipleOfEight(int(calculations[calculationIndex].segments.size())), int(calculations[calculationIndex].concentrations.size()));
                calculations[calculationIndex].sigmaProfilesOutdated = false;
            }


            // segment concentrations and the combinatorial contribution are calculated once per calculation
            // (tracked per calculation, so that calculations not run before or that failed are initialized later)
            // and again whenever concentrations, sigma profiles or size related parameters change
            const bool calculateSizeRelatedProperties = param.sw_alwaysCalculateSizeRelatedParameters == 1
                || !calculations[calculationIndex].sizeRelatedPropertiesInitialized
                || param.sw_reloadConcentrations == 1 || param.sw_reloadReferenceConcentrations == 1;

            if (calculateSizeRelatedProperties) {

                rescaleSegments(param, calculations[calculationIndex]);
                calculateSegmentConcentrations(calculations[calculationIndex]);
            }

#ifdef PRINT_DEBUG_INFO
            display("number of components: " + std::to_string(calculations[calculationIndex].components.size()) + "\n");
            display("number of segments: " + std::to_string(calculations[calculationIndex].segments.size()) + "\n");
            display("number of concentrations: " + std::to_string(calculations[calculationIndex].concentrations.size()) + "\n");
#endif

#ifdef MEASURE_TIME
            rescaleSegments_total_ms += (const unsigned long)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - rescaleSegments_last).count());
            std::chrono::high_resolution_clock::time_point calculateCombinatorial_last = std::chrono::high_resolution_clock::now();
#endif
            // recalculate combinatorial term if needed
            if (calculateSizeRelatedProperties) {
                calculateLnGammaCombinatorial(param, calculations[calculationIndex]);
                calculations[calculationIndex].sizeRelatedPropertiesInitialized = true;
            }

#ifdef MEASURE_TIME
            calculateCombinatorial_total_ms += (const unsigned long)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateCombinatorial_last).count());
            std::chrono::high_resolution_clock::time_point calculateResidual_last = std::chrono::high_resolution_clock::now();
#endif

            // calculate residual part
            calculateLnGammaResidual(param, calculations[calculationIndex]);

#ifdef MEASURE_TIME
            calculateResidual_total_ms += (const unsigned long)(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateResidual_last).count());
#endif

            calculations[calculationIndex].lnGammaTotal = calculations[calculationIndex].lnGammaCombinatorial + calculations[calculationIndex].lnGammaResidual;
#ifdef PRINT_DEBUG_INFO
            for (int i_concentration = 0; i_concentration < calculations[calculationIndex].originalNumberOfCalculations; i_concentration++) {
                for (int i_component = 0; i_component < calculations[calculationIndex].components.size(); i_component++) {
                    display("lnGammaTotal_" + std::to_string(calculationIndex) + "_" + std::to_string(i_concentration) + "_" + std::to_string(i_component) + \
                        ": " + std::to_string(calculations[calculationIndex].lnGammaTotal(i_concentration, i_component)) + "\n");
                }
            }
#endif

            bool calculateSolvationEnergies = param.dGsolv_E_gas.size() > 0;
            if (calculateSolvationEnergies){
                double kcalPerMol_per_Hartree = 2625.499639479/4.184;
                double reference_pressure = 101325; // Pa = 1 atm;
                std::vector<int> atomicNumbersWithout_dGsolv_tau = std::vector<int>();
                double approximate_dGsolv_tau = 0.0262; // median of other values
                calculation& thisCalculation = calculations[calculationIndex];
                for (int i_concentration = 0; i_concentration < thisCalculation.originalNumberOfCalculations; i_concentration++) {
                    if (thisCalculation.referenceStateType[i_concentration] == 4) {

                        if (thisCalculation.dGsolv.data() == nullptr) {
                            throw std::runtime_error("Calculation number " + std::to_string(thisCalculation.number) + " uses reference state type 4 (solvation energy), but no 'dGsolv' output array was given.");
                        }

                        // concentrations and temperatures are stored sorted internally
                        const int i_internal = thisCalculation.actualConcentrationIndices[i_concentration];
                        const std::vector<double>& thisConcentration = thisCalculation.concentrations[i_internal];

                        int i_solvent_component = -1;
                        for (int i_component = 0; i_component < thisCalculation.components.size(); i_component++) {
                            if (thisConcentration[i_component] == 1.0) {
                                i_solvent_component = i_component;
                                break;
                            }
                        }
                        if (i_solvent_component == -1) {
                            throw std::runtime_error("Calculation number " + std::to_string(thisCalculation.number) + ", concentration " + std::to_string(i_concentration) +
                                ": reference state type 4 (solvation energy) requires one component with a mole fraction of 1 (the solvent) and all other components with a mole fraction of 0.");
                        }
                        for (int i_component = 0; i_component < thisCalculation.components.size(); i_component++) {
                            double dGsolv = 0.0;
                            if (thisConcentration[i_component] == 0.0) {

                                const int i_molecule = thisCalculation.componentIndices[i_component];
                                const molecule& solute = *thisCalculation.components[i_component];

                                if (std::isnan(solute.epsilonInfinityTotalEnergy)) {
                                    throw std::runtime_error("The COSMO file of molecule '" + solute.name + "' does not contain the total energy, which is needed to calculate solvation energies.");
                                }

                                double RT = R_GAS_CONSTANT * thisCalculation.temperatures[i_internal];
                                double molar_volume_ideal_gas = RT / reference_pressure;
                                double RT_kcalPerMol = RT / (1000 * 4.184);

                                // all energies calculated below this line are in kcal/mol
                                double E_vdw = 0.0;
                                const segmentTypeCollection& segments = solute.segments;
                                std::unordered_map<int, double> areasByAtomicNumber;
                                for (int i_segment = 0; i_segment < int(segments.SegmentTypeAtomicNumber.size()); i_segment++) {
                                    int AN = getElementForParameterLookup(segments.SegmentTypeAtomicNumber[i_segment]);

                                    if (areasByAtomicNumber.find(AN) == areasByAtomicNumber.end())
                                        areasByAtomicNumber[AN] = 0.0;

                                    areasByAtomicNumber[AN] += segments.SegmentTypeAreas[i_segment][0];
                                }

                                for (auto& it : areasByAtomicNumber) {
                                    double this_atom_dGsolv_tau = std::abs(param.dGsolv_tau[it.first]);
                                    if (this_atom_dGsolv_tau == 0.0) {
                                        this_atom_dGsolv_tau = approximate_dGsolv_tau;
                                        atomicNumbersWithout_dGsolv_tau.push_back(it.first);
                                    }
                                    
                                    E_vdw += this_atom_dGsolv_tau * it.second;
                                }
                                double referenceStateCorrection = RT_kcalPerMol * log(molar_volume_ideal_gas / (thisCalculation.components[i_solvent_component]->molarVolumeAt25C / 1E6));

                                // dGsolv_E_gas and dGsolv_numberOfAtomsInRing are given per loaded molecule
                                double E_diel = (solute.epsilonInfinityTotalEnergy - param.dGsolv_E_gas[i_molecule]) * kcalPerMol_per_Hartree;
                                double mu_liquid = RT_kcalPerMol * thisCalculation.lnGammaTotal(i_concentration, i_component);
                                double E_ring = param.dGsolv_omega_ring * param.dGsolv_numberOfAtomsInRing[i_molecule];
                                dGsolv = E_diel + mu_liquid - E_vdw - E_ring - referenceStateCorrection - param.dGsolv_eta;
#ifdef PRINT_DEBUG_INFO
                                display("dGsolv_1_" + std::to_string(calculationIndex) + ": " + std::to_string(dGsolv) + "\n");
                                display("E_diel: " + std::to_string(E_diel) + "\n");
                                display("mu_liquid: " + std::to_string(mu_liquid) + "\n");
                                display("E_vdw: " + std::to_string(E_vdw) + "\n");
                                display("E_ring: " + std::to_string(E_ring) + "\n");
                                display("referenceStateCorrection: " + std::to_string(referenceStateCorrection) + "\n");
                                display("dGsolv_eta: " + std::to_string(param.dGsolv_eta) + "\n");
#endif

                            }

                            thisCalculation.dGsolv(i_concentration, i_component) = dGsolv;
                        }
                    }
                }
                if (atomicNumbersWithout_dGsolv_tau.size() > 0) {
                    sort(atomicNumbersWithout_dGsolv_tau.begin(), atomicNumbersWithout_dGsolv_tau.end());
                    atomicNumbersWithout_dGsolv_tau.erase(unique(atomicNumbersWithout_dGsolv_tau.begin(), atomicNumbersWithout_dGsolv_tau.end()), atomicNumbersWithout_dGsolv_tau.end());
                    std::string ANs = "";
                    for (int i_AN = 0; i_AN < atomicNumbersWithout_dGsolv_tau.size(); i_AN++) {
                        ANs += std::to_string(atomicNumbersWithout_dGsolv_tau[i_AN]) + ";";
                    }
                    if (param.sw_dGsolv_calculation_strict == 1) {
                        throw std::runtime_error("For the following atomic numbers not all parameters are available for the calculation of solvation energies: " + ANs);
                    }
                    else {
                        addWarning(" - For the following atomic numbers not all parameters are available for the calculation of solvation energies, an estimate was used: " + ANs);
                    }

                }
            }
        });
    }
    sigmaProfilesUsedInCalculation = true;
    e.rethrow();
}
