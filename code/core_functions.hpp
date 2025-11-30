/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/


#pragma once

#include "interaction_matrix.hpp"
#include "contact_statistics.hpp"
#include "COSMOfile_functions.hpp"
#include <stdexcept>
#include <functional>
#include <iostream>


#if defined(MEASURE_TIME) 
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

    oneIteration_total_ms += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - oneIteration_last).count();

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

float hsum_ps_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);
    __m128 sums = _mm_add_ps(v, shuf);
    shuf = _mm_movehl_ps(shuf, sums);
    sums = _mm_add_ss(sums, shuf);
    return _mm_cvtss_f32(sums);
}

#if defined(__AVX__) || defined(__FMA__)
float hsum256_ps_avx(__m256 v) {
    __m128 vlow = _mm256_castps256_ps128(v);
    __m128 vhigh = _mm256_extractf128_ps(v, 1);
    vlow = _mm_add_ps(vlow, vhigh);
    return hsum_ps_sse3(vlow);
}
#endif

void initialize(parameters& param, bool showBinarySpecs = true) {

    n_ex = 0;
    molecules.clear();
    calculations.clear();

    param.ChargeRaster.clear();
    param.exp_param.clear();
    param.R_i = std::vector<double>(118, 0.0);
    param.R_i_COSMO = std::vector<double>(118, 0.0);
    param.HBClassElmnt = std::vector<int>(300, 0);

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

    if (showBinarySpecs)
        display("\nBINARY SPECS\n-------------------------\n" + compilation_mode + "\n" + OPENMP_parallelization + "\n" + vectorization_level + "\n-------------------------\n\n");
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
    if (calculateSolvationEnergies) {
        if (_molecule.qmMethod != "DFT_CPCM_BP86_def2-TZVP+def2-TZVPD_SP" && _molecule.qmMethod != "DFT_BP86_def2-TZVPD_SP") {
            if (param.sw_dGsolv_calculation_strict == 1) {
                throw std::runtime_error("The QSPR model for the molar volume only works for the quantum chemistry method DFT_BP86_def2-TZVPD_SP");
            }
            else {
                warnings.push_back(" - The QSPR model for the molar volume was parametrized using a different quantum chemistry method than the one you are using. Recommended method: DFT_BP86_def2-TZVPD_SP");
            }
        }

        int numberOfSiAtoms = 0;
        int numberOfHAtoms = 0;
        int numberOfOAtoms = 0;
        for (int i = 0; i < numberOfAtoms; i++) {
            if (_molecule.atomAtomicNumbers[i] == 14)
                numberOfSiAtoms += 1;
            else if (_molecule.atomAtomicNumbers[i] == 1 || _molecule.atomAtomicNumbers[i] > 100)
                numberOfHAtoms += 1;
            else if (_molecule.atomAtomicNumbers[i] == 8)
                numberOfOAtoms += 1;
        }

        if (numberOfAtoms == 3 && numberOfHAtoms == 2 && numberOfOAtoms == 1) {
            _molecule.molarVolumeAt25C = 18.06863632;
        }
        else {
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
    float sigmaLeft = -1;
    float sigmaRight = -1;
    double AsigmaLeft = -1;
    double AsigmaRight = -1;

    float sigmaCorrLeft = -1;
    float sigmaCorrRight = -1;

    double AsigmaLeftSigmaCorrLeft = -1;
    double AsigmaLeftSigmaCorrRight = -1;
    double AsigmaRightSigmaCorrLeft = -1;
    double AsigmaRightSigmaCorrRight = -1;

    unsigned short ind_SigmaCorr_left = 0;

    for (int j = 0; j < numberOfSegments; j++) {

        unsigned short ind_Sigma_left = int((averagedSigmas(j) - param.sigmaMin) / param.sigmaStep);

        sigmaLeft = (float)param.ChargeRaster[ind_Sigma_left];
        sigmaRight = (float)param.ChargeRaster[ind_Sigma_left + 1];

        AsigmaRight = _molecule.segmentAreas(j) * (averagedSigmas(j) - sigmaLeft) / param.sigmaStep;
        AsigmaLeft = _molecule.segmentAreas(j) * (sigmaRight - averagedSigmas(j)) / param.sigmaStep;


        if (calculateMisfitCorrelation == true) {
            ind_SigmaCorr_left = int((averagedSigmaCorrs(j) - param.sigmaMin) / param.sigmaStep);

            sigmaCorrLeft = (float)param.ChargeRaster[ind_SigmaCorr_left];
            sigmaCorrRight = (float)param.ChargeRaster[ind_SigmaCorr_left + 1];

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

    int numberOfAtoms = int(_molecule.atomAtomicNumbers.size());
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
        for (int segmentIndexI = 0; segmentIndexI < numberOfSegments; segmentIndexI++) {
            int atomIndex = segmentAtomIndex(segmentIndexI);
            areaOfAtom(atomIndex) = areaOfAtom(atomIndex) + surfaceAreasInAu(segmentIndexI);
            segmentsperAtom(atomIndex) = segmentsperAtom(atomIndex) + 1;
        }


        for (int segmentIndexI = 0; segmentIndexI < numberOfSegments; segmentIndexI++) {
            int atomIndex = segmentAtomIndex(segmentIndexI);
            for (int k = 0; k < 3; k++) {
                distanceSegmentAtom(k) = abs(surfaceCoordinatesInAu(segmentIndexI, k) - atomCoordinatesInAu(atomIndex, k));
            }
            double distanceSegmentAtom_norm = distanceSegmentAtom(0) + distanceSegmentAtom(1) + distanceSegmentAtom(2);
            distanceSegmentAtom = distanceSegmentAtom / distanceSegmentAtom_norm;
            polarizabilityTensor = extractAtomicPolarizabilityTensor(_molecule, atomIndex);

            Eigen::RowVector3d distanceSegmentAtom_row = distanceSegmentAtom.transpose();
            Eigen::RowVector3d tmp_vec = distanceSegmentAtom_row * polarizabilityTensor;

            polarizability_projection = tmp_vec.norm();

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

    int numberOfAtoms = int(newMolecule.atomAtomicNumbers.size());

    float sumOfScreeningCharge = float((newMolecule.segmentAreas.array() * newMolecule.segmentSigmas.array()).matrix().sum());

    newMolecule.moleculeCharge = (signed char)(std::round(-1.0f * sumOfScreeningCharge));

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

        for (int atomIndexI = 0; atomIndexI < numberOfAtoms; atomIndexI++) {
            double minimumDistanceAtomIAtomJSquared = 10e14;
            int closestAtomIndex = -1;
            if (newMolecule.atomAtomicNumbers(atomIndexI) != 1) {
                continue;
            }

            for (int atomIndexJ = 0; atomIndexJ < numberOfAtoms; atomIndexJ++) {

                if (atomIndexI == atomIndexJ) {
                    continue;
                }

                // search for next atom that is not a Hydrogen
                if (newMolecule.atomAtomicNumbers(atomIndexJ) == 1 && newMolecule.atomAtomicNumbers(atomIndexJ) < 100) {
                    continue;
                }

                if (distanceAtomAtomSquared(atomIndexI, atomIndexJ) < minimumDistanceAtomIAtomJSquared) {
                    minimumDistanceAtomIAtomJSquared = distanceAtomAtomSquared(atomIndexI, atomIndexJ);
                    closestAtomIndex = atomIndexJ;
                }
            }
            newMolecule.atomAtomicNumbers(atomIndexI) = 100 + newMolecule.atomAtomicNumbers(closestAtomIndex);
            param.R_i_COSMO[100 + newMolecule.atomAtomicNumbers(closestAtomIndex)] = param.R_i_COSMO[1];
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

    newMolecule.clear_unneeded_matrices(param.sw_alwaysReloadSigmaProfiles);
    newMolecule.segments.shrink_to_fit();

#ifdef DEBUG_INFO
    newMolecule.segments.sort();
    WriteExtendedSigmaProfiletoFile(newMolecule.name + ".extsp", newMolecule.segments);
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
            molecule _molecule = *molecules[i];
            int previousNumberOfSegmentTypes = int(_molecule.segments.size());
            _molecule.segments.clear();
            averageAndClusterSegments(param, _molecule, previousNumberOfSegmentTypes);
            });
    }
    e.rethrow();
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
            molecules[i]->Area = (4 * PI * R_i * R_i);
            molecules[i]->Volume = (4.0 / 3 * PI * R_i * R_i * R_i);
        }
    }
#ifdef MEASURE_TIME
    rescaleSegments_total_ms += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - rescaleSegments_last).count();
#endif
}

void rescaleSegments(parameters& param, calculation& _calculation) {
    // rescale segments for monoatomic cations only if present
    if (_calculation.segments.numberOfSegmentsForGroup[3] != 0) {

        // monoatomic cation
        std::unordered_map<unsigned short, int[2]> segmentIndicesBelongingToASpecificAtomicNumber;
        for (int j = _calculation.segments.lowerBoundIndexForGroup[3];
            j < _calculation.segments.upperBoundIndexForGroup[3]; j++) {

            unsigned short AN = _calculation.segments.SegmentTypeAtomicNumber[j];

            if (segmentIndicesBelongingToASpecificAtomicNumber.find(AN) == segmentIndicesBelongingToASpecificAtomicNumber.end()) {
                segmentIndicesBelongingToASpecificAtomicNumber[AN][0] = (unsigned short)j;
            }
            else {
                segmentIndicesBelongingToASpecificAtomicNumber[AN][1] = (unsigned short)j;
            }

        }
        for (auto it = segmentIndicesBelongingToASpecificAtomicNumber.begin(); it != segmentIndicesBelongingToASpecificAtomicNumber.end(); ++it) {

            int AN = it->first;
            int ind_left = it->second[0];
            int ind_right = it->second[1];

            if (ind_right - ind_left > 1) {
                throw std::runtime_error("More than 2 segmentTypes was found for monoatomic cation: " + std::to_string(AN) + ".");
            }

            int ind_molecule = -1;

            for (int k = 0; _calculation.segments.SegmentTypeAreas[ind_left].size(); k++) {
                if (_calculation.segments.SegmentTypeAreas[ind_left][k] > 0.0) {
                    ind_molecule = k;
                    break;
                }
            }

            int screeningCharge = int(_calculation.components[ind_molecule]->moleculeCharge * -1);

            double newArea = (4 * PI * param.R_i[AN] * param.R_i[AN]);
            double newSigma = screeningCharge / newArea;

            unsigned short ind_Sigma_left = int((newSigma - param.sigmaMin) / param.sigmaStep);

            double sigmaLeft = param.ChargeRaster[ind_Sigma_left];
            double sigmaRight = param.ChargeRaster[ind_Sigma_left + 1];

            double AsigmaRight = newArea * (newSigma - sigmaLeft) / param.sigmaStep;
            double AsigmaLeft = newArea * (sigmaRight - newSigma) / param.sigmaStep;

            _calculation.segments.SegmentTypeAreas[ind_left][ind_molecule] = AsigmaLeft;
            _calculation.segments.SegmentTypeSigma[ind_left] = (float)sigmaLeft;

            _calculation.segments.SegmentTypeAreas[ind_right][ind_molecule] = AsigmaRight;
            _calculation.segments.SegmentTypeSigma[ind_right] = (float)sigmaRight;

        }
    }
}


void calculateSegmentConcentrations(calculation& _calculation) {
    // calculate the mole fraction of segments for each concentration
    for (int j = 0; j < _calculation.concentrations.size(); j++) {

        int firstNonZeroSegmentIndex = 0;
        int lastNonZeroSegmentIndex = int(_calculation.concentrations.size()) - 1;

        std::vector<float> segmentConcentration(_calculation.segments.size(), 0.0f);

        double sumAreaSegmentsConcentrationj = 0.0;

        for (int k = 0; k < _calculation.segments.size(); k++) {

            double areaSegmentK = 0.0;

            for (int m = 0; m < _calculation.components.size(); m++) {
                double thisArea = _calculation.concentrations[j][m] * _calculation.segments.SegmentTypeAreas[k][m];
                areaSegmentK += thisArea;
                sumAreaSegmentsConcentrationj += thisArea;
            }

            segmentConcentration[k] = float(areaSegmentK);
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

            _calculation.segmentConcentrations(k, j) = float(segmentConcentration[k] / sumAreaSegmentsConcentrationj);
        }

        _calculation.lowerBoundIndexForCOSMOSPACECalculation[j] = RoundDownToNextMultipleOfEight(firstNonZeroSegmentIndex);
        _calculation.upperBoundIndexForCOSMOSPACECalculation[j] = RoundUpToNextMultipleOfEight(lastNonZeroSegmentIndex + 1);
    }
}

void finishCalculationInitiation(calculation& _calculation) {

    if (_calculation.concentrations.size() > 65535) {
        throw std::runtime_error("Too many calculations, other datatype would be necessary for newCalculation.referenceStates to cope with this amount. (unsigned short used so far allowing for up to 65535)");
    }

    // check charge of all concentrations
    for (int j = 0; j < _calculation.concentrations.size(); j++) {

        double mix_chrg = 0;
        for (int k = 0; k < _calculation.components.size(); k++) {
            mix_chrg = mix_chrg + _calculation.components[k]->moleculeCharge * _calculation.concentrations[j][k];
        }

        if (abs(mix_chrg) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
            throw std::runtime_error("For calculation number " + std::to_string(_calculation.number) + " the mixture number " + std::to_string(j) + " is not electroneutral. residual charge: " + std::to_string(abs(mix_chrg)));
        }
    }

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
    _calculation.segmentGammas = Eigen::MatrixXf::Constant(RoundUpToNextMultipleOfEight(int(_calculation.segments.size())), int(_calculation.concentrations.size()), 1.0f);
    _calculation.segmentConcentrations = Eigen::MatrixXf::Zero(RoundUpToNextMultipleOfEight(int(_calculation.segments.size())), int(_calculation.concentrations.size()));

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

    std::vector<float> qi_std(_calculation.components.size());

    /* Calculation of molecule specific parameters for combinatorial contribution */
    for (int i = 0; i < _calculation.components.size(); i++) {
        qi_std[i] = float(_calculation.components[i]->Area / param.comb_SG_A_std);
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
    else if (param.sw_combTerm == 4) { // not tested yet
        /*  mod. Staverman-Guggenheim by Grensemann published in Grensemann & Gmehling (2005) especially developed for COSMO-RS */
        for (int i = 0; i < _calculation.concentrations.size(); i++) {

            double sum_qi_div_xi = 0;
            double sum_ri_div_xi = 0;
            double sum_qi_times_xi = 0;
            double sum_ri_times_xi = 0;

            for (int j = 0; j < _calculation.components.size(); j++) {
                sum_qi_div_xi = sum_qi_div_xi + (_calculation.components[j]->Area / _calculation.concentrations[i][j]);
                sum_ri_div_xi = sum_ri_div_xi + (_calculation.components[j]->Volume / _calculation.concentrations[i][j]);
                sum_qi_times_xi = sum_qi_times_xi + (_calculation.components[j]->Area * _calculation.concentrations[i][j]);
                sum_ri_times_xi = sum_ri_times_xi + (_calculation.components[j]->Volume * _calculation.concentrations[i][j]);
            }

            std::vector<double> lnGamaForThisCalculation;
            for (int j = 0; j < _calculation.components.size(); j++) {

                double qi_hash = (1 - _calculation.concentrations[i][j]) * sum_qi_div_xi;
                double ri_hash = (1 - _calculation.concentrations[i][j]) * sum_ri_div_xi;

                double qi_min = std::min(_calculation.components[j]->Area, qi_hash);
                double qi_max = std::max(_calculation.components[j]->Area, qi_hash);

                double ri_min = std::min(_calculation.components[j]->Volume, ri_hash);
                double ri_max = std::max(_calculation.components[j]->Volume, ri_hash);

                double Fi = pow(_calculation.components[j]->Area / sum_qi_times_xi, param.comb_SGG_lambda * (1 - (qi_min / qi_max)));
                double Vi = pow(_calculation.components[j]->Volume / sum_ri_times_xi, param.comb_SGG_beta * (1 - (ri_min / ri_max)));

                lnGamaForCalculations(i, j) = 1 - Vi - log(Vi) + 1 - (Vi / Fi) - log(Vi / Fi);
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

            float multiplier = 1.0f;
            // for an ion with reference state PureComponentsOnlyNeutral
            if (_calculation.referenceStateType[h] == 1 && _calculation.referenceStateCalculationIndices[h][j] == -1) {
                multiplier = NAN;
            }
            _calculation.lnGammaCombinatorial(h, j) = float(multiplier * lnGamaForCalculations(i, j));

            // only subtract reference state value if reference state type is not COSMO
            if (_calculation.referenceStateType[h] != 3 && _calculation.referenceStateType[h] != 4) {
                int indexOfReferenceStateCalculation = _calculation.actualConcentrationIndices[_calculation.referenceStateCalculationIndices[h][j]];

                _calculation.lnGammaCombinatorial(h, j) -= float(lnGamaForCalculations(indexOfReferenceStateCalculation, j));
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

    Eigen::MatrixXf temporary_lnGammaMolecule = Eigen::MatrixXf::Zero(_calculation.concentrations.size(), _calculation.components.size());

    Eigen::Tensor<float, 4, Eigen::RowMajor> temporary_averageInteractionEnergies;
    Eigen::Tensor<float, 3, Eigen::RowMajor> temporary_partialMolarEnergies;

    if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {

        temporary_averageInteractionEnergies = Eigen::Tensor<float, 4, Eigen::RowMajor>(int(_calculation.concentrations.size()),
            param.numberOfPartialInteractionMatrices + 1, // +1 because A_int is the first one
            int(_calculation.components.size()),
            int(_calculation.components.size()));

        temporary_averageInteractionEnergies.setZero();

        if (param.sw_calculateContactStatisticsAndAdditionalProperties == 2) {
            temporary_partialMolarEnergies = Eigen::Tensor<float, 3, Eigen::RowMajor>(int(_calculation.concentrations.size()),
                param.numberOfPartialInteractionMatrices + 1, // +1 because A_int is the first one
                int(_calculation.components.size()));

            temporary_partialMolarEnergies.setZero();

        }

    }


    const double dampingFactor = 0.6;
    const double dampingFactorComplement = 1 - dampingFactor;
    const double convergenceThreshhold = 0.00001;   // originally this was 0.0000001;
                                                    // however when using AVX and vectorized code with floats
                                                    // the accuracy of the calculation drops
                                                    // without increasing this value the gammas do not converge
                                                    // however since the convergence criteria was changed from
                                                    // gamma to log(gamma) the larger convergence threshhold 
                                                    // still leads to numerically very similar results

    const int maximumNumberOfIterations = 50000;
    double tempSum = 0.0;


    Eigen::MatrixXf A_int = Eigen::MatrixXf::Zero(numberOfSegments, numberOfSegments);
    Eigen::MatrixXf Tau = Eigen::MatrixXf::Zero(nMultipleOfEight, numberOfSegments);
    Eigen::MatrixXf TauX = Eigen::MatrixXf::Zero(nMultipleOfEight, numberOfSegments);

    std::vector<Eigen::MatrixXd> partialInteractionMatrices;

    if (param.numberOfPartialInteractionMatrices > 0) {
        for (int h = 0; h < param.numberOfPartialInteractionMatrices; h++) {
            partialInteractionMatrices.push_back(Eigen::MatrixXd::Zero(numberOfSegments, numberOfSegments));
        }
    }

    for (int g = 0; g < _calculation.TauConcentrationIndices.size(); g++) {

        // conditions
        float temperature = _calculation.TauTemperatures[g];

#ifdef MEASURE_TIME
        std::chrono::high_resolution_clock::time_point calculateTau_last = std::chrono::high_resolution_clock::now();
#endif
        // calculate interaction matrices
        if (param.numberOfPartialInteractionMatrices > 0) {
            for (int h = 0; h < param.numberOfPartialInteractionMatrices; h++) {
                partialInteractionMatrices[h].setZero();
            }
        }

        calculateInteractionMatrix(_calculation.segments, A_int, partialInteractionMatrices, param, temperature);

        // if the full interaction matrix is needed, fill upper right half
        if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {
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
        const float minus_div_RT = float(-1.0f / (R_GAS_CONSTANT * temperature));

        float val;
        float* Tau_1D = &(Tau(0, 0));
        for (int i = 0; i < numberOfSegments; i++) {
            columnSum = i * nMultipleOfEight;
            for (int j = i; j < numberOfSegments; j++) {
                // always j >= i
                val = expf(A_int(j, i) * minus_div_RT);

                idx = columnSum + j;
                Tau_1D[idx] = val;
                Tau_1D[j * nMultipleOfEight + i] = val;
            }
        }

#ifdef MEASURE_TIME
        calculateTau_total_ms += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateTau_last).count();
#endif

#ifdef DEBUG_INFO
        Eigen::MatrixXd Tau_d = Tau.cast<double>();
        WriteEigenMatrixtoFile("Tau_" + std::to_string(_calculation.number) + "_" + std::to_string(g), Tau_d);
#endif

        for (int h = 0; h < _calculation.TauConcentrationIndices[g].size(); h++) {

            int i = _calculation.TauConcentrationIndices[g][h];

            float* gammas = &(_calculation.segmentGammas(0, i));
#ifdef MEASURE_TIME
            std::chrono::high_resolution_clock::time_point calculateCOSMOSPACE_last = std::chrono::high_resolution_clock::now();
#endif

            int lowerBoundIndexForCOSMOSPACECalculation = _calculation.lowerBoundIndexForCOSMOSPACECalculation[i];
            int upperBoundIndexForCOSMOSPACECalculation = _calculation.upperBoundIndexForCOSMOSPACECalculation[i];

            { // as local as possible scope for all variables used in your loops

                const float* vTau_1D = &(Tau(0, 0));
                const float* vX = &(_calculation.segmentConcentrations(0, i));
                float* vTauX_1D = &(TauX(0, 0));
                int idx = 0;
                int columnSum = 0;


                for (int j = 0; j < numberOfSegments; j++) {

                    columnSum = j * nMultipleOfEight;

#if defined(__AVX__) || defined(__FMA__)//AVX
                    for (int k = lowerBoundIndexForCOSMOSPACECalculation; k < upperBoundIndexForCOSMOSPACECalculation; k += 8) {
                        idx = columnSum + k;

                        _mm256_store_ps(vTauX_1D + idx, _mm256_mul_ps(_mm256_load_ps(vTau_1D + idx), _mm256_load_ps(vX + k)));
                    }
#else //SSE3
                    for (int k = lowerBoundIndexForCOSMOSPACECalculation; k < upperBoundIndexForCOSMOSPACECalculation; k += 4) {
                        idx = columnSum + k;

                        _mm_store_ps(vTauX_1D + idx, _mm_mul_ps(_mm_load_ps(vTau_1D + idx), _mm_load_ps(vX + k)));
                    }
#endif

                }
            }

#ifdef DEBUG_INFO
            Write1DArraytoFile<float>("X_" + std::to_string(_calculation.number) + "_" + std::to_string(g) + "_" + std::to_string(i), &(_calculation.segmentConcentrations(0, i)), 1, nMultipleOfEight);
            Eigen::MatrixXd TauX_d = TauX.cast<double>();
            WriteEigenMatrixtoFile("TauX_" + std::to_string(_calculation.number) + "_" + std::to_string(g) + "_" + std::to_string(i), TauX_d);
#endif

            int numberOfConvergedGammas = 0;
            int numberOfIteration = 0;

            while (numberOfConvergedGammas != numberOfSegments) {

                numberOfIteration++;

                if (numberOfIteration > maximumNumberOfIterations) {

                    if (param.sw_skip_COSMOSPACE_errors == 0) {

                        std::string information = "";
                        information += "\nSYSTEM_COMPONENTS:\n";
                        for (int j = 0; j < _calculation.components.size(); j++) {
                            information += " -" + _calculation.components[j]->name + "\n";
                        }

                        float* Tau_1D = &(Tau(0, 0));
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
                        display(information + "\nCOSMOSPACE did not converge for calculation " + std::to_string(_calculation.number) + " on concentration " + std::to_string(i) + ": maximum number of iterations reached. (" + std::to_string(numberOfConvergedGammas) + "/" + std::to_string(numberOfSegments) + ")");
                    }
                    else {
                        // set all segment gammas to 1 as initial point for next execution
                        _calculation.segmentGammas = Eigen::MatrixXf::Constant(RoundUpToNextMultipleOfEight(int(_calculation.segments.size())), int(_calculation.concentrations.size()), 1.0f);
                    }

                    throw std::runtime_error("COSMOSPACE did not converge for calculation " + std::to_string(_calculation.number) + " on concentration " + std::to_string(i) + ": maximum number of iterations reached. (" + std::to_string(numberOfConvergedGammas) + "/" + std::to_string(numberOfSegments) + ")");
                }

                numberOfConvergedGammas = 0;

                int columnSum = 0;

                for (int j = 0; j < numberOfSegments; j++) {
                    tempSum = 0.0;

                    columnSum = j * nMultipleOfEight;

                    { // as local as possible scope for all variables used in your loops

                        const float* vTauX = &(TauX(0, 0));
                        float* vGammas = gammas;

#if defined(__FMA__)
                        __m256 tempVectorSum = _mm256_set_ps(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

                        for (int k = lowerBoundIndexForCOSMOSPACECalculation; k < upperBoundIndexForCOSMOSPACECalculation; k += 8) {
                            tempVectorSum = _mm256_fmadd_ps(_mm256_load_ps(vTauX + columnSum + k), _mm256_load_ps(vGammas + k), tempVectorSum); //FMA
                        }

                        tempSum = hsum256_ps_avx(tempVectorSum);

#elif defined(__AVX__)
                        __m256 tempVectorSum = _mm256_set_ps(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

                        for (int k = lowerBoundIndexForCOSMOSPACECalculation; k < upperBoundIndexForCOSMOSPACECalculation; k += 8) {
                            //tempSum += hsum256_ps_avx(_mm256_mul_ps(_mm256_load_ps(vTauX + columnSum + k), _mm256_load_ps(vGammas + k))); // less accurate
                            tempVectorSum = _mm256_add_ps(_mm256_mul_ps(_mm256_load_ps(vTauX + columnSum + k), _mm256_load_ps(vGammas + k)), tempVectorSum); //AVX
                        }

                        tempSum = hsum256_ps_avx(tempVectorSum);
#else//SSE3
                        __m128 tempVectorSum = _mm_setr_ps(0.0, 0.0, 0.0, 0.0);

                        for (int k = lowerBoundIndexForCOSMOSPACECalculation; k < upperBoundIndexForCOSMOSPACECalculation; k += 4) {
                            tempVectorSum = _mm_add_ps(_mm_mul_ps(_mm_load_ps(vTauX + columnSum + k), _mm_load_ps(vGammas + k)), tempVectorSum); //SSE
                        }

                        tempSum = hsum_ps_sse3(tempVectorSum);
#endif
                    }

                    // newGamma
                    tempSum = 1 / tempSum;


                    // apply damping
                    if (numberOfIteration > 200) {
                        tempSum = dampingFactor * tempSum + dampingFactorComplement * gammas[j];
                    }

                    // check convergence criteria
                    if (abs(log(tempSum) - log(gammas[j])) <= convergenceThreshhold) {
                        numberOfConvergedGammas++;
                    }

                    gammas[j] = float(tempSum);
                }
            }

#ifdef DEBUG_INFO
            display("COSMOSPACE converged with this number of iterations: " + std::to_string(numberOfIteration) + "\n");
            Write1DArraytoFile<float>("Gammas_" + std::to_string(_calculation.number) + "_" + std::to_string(g) + "_" + std::to_string(i), gammas, 1, numberOfSegments);
#endif

#ifdef MEASURE_TIME
            calculateCOSMOSPACE_total_ms += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateCOSMOSPACE_last).count();
            std::chrono::high_resolution_clock::time_point calculateContactStatistics_last = std::chrono::high_resolution_clock::now();
#endif


            if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {
                calculateContactStatistics(_calculation, A_int, partialInteractionMatrices, Tau, gammas, int(i), temporary_averageInteractionEnergies, temporary_partialMolarEnergies, param);
            }

#ifdef MEASURE_TIME
            calculateContactStatistics_total_ms += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateContactStatistics_last).count();
            std::chrono::high_resolution_clock::time_point calculateGammasForMolecules_last = std::chrono::high_resolution_clock::now();
#endif
            double div_Aeff = 1 / param.Aeff;

            for (int j = 0; j < _calculation.components.size(); j++) {

                double lnGammaMolecule = 0;

                for (int k = 0; k < numberOfSegments; k++) {
                    lnGammaMolecule += _calculation.segments.SegmentTypeAreas[k][j] * div_Aeff * log(gammas[k]); // log(gammas[k]) is calculated more than once although not necessary
                }

                temporary_lnGammaMolecule(i, j) = float(lnGammaMolecule);
            }

#ifdef MEASURE_TIME
            calculateGammasForMolecules_total_ms += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateGammasForMolecules_last).count();
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


            float multiplier = 1.0f;
            // for an ion with reference state PureComponentsOnlyNeutral
            if (_calculation.referenceStateType[h] == 1 && _calculation.referenceStateCalculationIndices[h][j] == -1) {
                multiplier = NAN;
            }
            _calculation.lnGammaResidual(h, j) = multiplier * temporary_lnGammaMolecule(i, j);

            if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {
                for (int k = 0; k < _calculation.components.size(); k++) {
                    _calculation.averageSurfaceEnergies(h, 0, j, k) = multiplier * temporary_averageInteractionEnergies(i, 0, j, k);

                    for (int l = 0; l < param.numberOfPartialInteractionMatrices; l++) {
                        _calculation.averageSurfaceEnergies(h, l + 1, j, k) = multiplier * temporary_averageInteractionEnergies(i, l + 1, j, k);
                    }
                }

                if (param.sw_calculateContactStatisticsAndAdditionalProperties == 2) {
                    _calculation.partialMolarEnergies(h, 0, j) = multiplier * temporary_partialMolarEnergies(i, 0, j);

                    for (int l = 0; l < param.numberOfPartialInteractionMatrices; l++) {
                        _calculation.partialMolarEnergies(h, l + 1, j) = multiplier * temporary_partialMolarEnergies(h, l + 1, j);
                    }
                }
            }

            // only subtract reference state value if reference state type is not COSMO
            if (_calculation.referenceStateType[h] != 3 && _calculation.referenceStateType[h] != 4) {
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
    calculateGammasForMolecules_total_ms += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateGammasForMolecules_last).count();
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
            if (param.sw_alwaysReloadSigmaProfiles == 1 && n_ex > 3) {

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
                calculations[calculationIndex].segmentGammas = Eigen::MatrixXf::Constant(RoundUpToNextMultipleOfEight(int(calculations[calculationIndex].segments.size())), int(calculations[calculationIndex].concentrations.size()), 1.0f);
                calculations[calculationIndex].segmentConcentrations = Eigen::MatrixXf::Zero(RoundUpToNextMultipleOfEight(int(calculations[calculationIndex].segments.size())), int(calculations[calculationIndex].concentrations.size()));
            }


            if (param.sw_alwaysCalculateSizeRelatedParameters == 1 || (param.sw_alwaysCalculateSizeRelatedParameters == 0 && n_ex == 3) \
                || (param.sw_alwaysReloadSigmaProfiles == 1 && n_ex > 3) || param.sw_reloadConcentrations == 1 || param.sw_reloadReferenceConcentrations == 1) {

                rescaleSegments(param, calculations[calculationIndex]);
                calculateSegmentConcentrations(calculations[calculationIndex]);
            }

#ifdef DEBUG_INFO
            display("number of components: " + std::to_string(calculations[calculationIndex].components.size()) + "\n");
            display("number of segments: " + std::to_string(calculations[calculationIndex].segments.size()) + "\n");
            display("number of concentrations: " + std::to_string(calculations[calculationIndex].concentrations.size()) + "\n");
#endif

#ifdef MEASURE_TIME
            rescaleSegments_total_ms += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - rescaleSegments_last).count();
            std::chrono::high_resolution_clock::time_point calculateCombinatorial_last = std::chrono::high_resolution_clock::now();
#endif
            // recalculate combinatorial term if needed
            if ((param.sw_alwaysCalculateSizeRelatedParameters == 0 && n_ex == 3) || param.sw_alwaysCalculateSizeRelatedParameters == 1 || param.sw_reloadConcentrations == 1 || param.sw_reloadReferenceConcentrations == 1) {
                calculateLnGammaCombinatorial(param, calculations[calculationIndex]);;
            }

#ifdef MEASURE_TIME
            calculateCombinatorial_total_ms += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateCombinatorial_last).count();
            std::chrono::high_resolution_clock::time_point calculateResidual_last = std::chrono::high_resolution_clock::now();
#endif

            // calculate residual part
            calculateLnGammaResidual(param, calculations[calculationIndex]);
            //  using std::cout;
            //  using std::endl;
            //  std::cout << "residual is calculated";

#ifdef MEASURE_TIME
            calculateResidual_total_ms += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - calculateResidual_last).count();
#endif

            calculations[calculationIndex].lnGammaTotal = calculations[calculationIndex].lnGammaCombinatorial + calculations[calculationIndex].lnGammaResidual;

            bool calculateSolvationEnergies = param.dGsolv_E_gas.size() > 0;
            if (calculateSolvationEnergies) {
                double kcalPerMol_per_Hartree = 2625.499639479 / 4.184;
                double reference_pressure = 101325; // Pa = 1 atm;
                std::vector<int> atomicNumbersWithout_dGsolv_tau = std::vector<int>();
                double approximate_dGsolv_tau = 0.0262; // median of other values
                for (int i_concentration = 0; i_concentration < calculations[calculationIndex].originalNumberOfCalculations; i_concentration++) {
                    if (calculations[calculationIndex].referenceStateType[i_concentration] == 4) {
                        int i_solvent_component = -1;
                        for (int i_component = 0; i_component < calculations[calculationIndex].components.size(); i_component++) {
                            if (calculations[calculationIndex].concentrations[i_concentration][i_component] == 1.0f) {
                                i_solvent_component = i_component;
                                break;
                            }
                        }
                        for (int i_component = 0; i_component < calculations[calculationIndex].components.size(); i_component++) {
                            double dGsolv = 0.0;
                            if (calculations[calculationIndex].concentrations[i_concentration][i_component] == 0.0f) {

                                double RT = R_GAS_CONSTANT * calculations[calculationIndex].temperatures[i_concentration];
                                double molar_volume_ideal_gas = RT / reference_pressure;
                                double RT_kcalPerMol = RT / (1000 * 4.184);

                                // all energies calculated below this line are in kcal/mol
                                double E_vdw = 0.0;
                                segmentTypeCollection segments = calculations[calculationIndex].components[i_component]->segments;
                                std::unordered_map<int, double> areasByAtomicNumber;
                                for (int i_segment = 0; i_segment < segments.size(); i_segment++) {
                                    int AN = segments.SegmentTypeAtomicNumber[i_segment];

                                    if (areasByAtomicNumber.find(AN) == areasByAtomicNumber.end())
                                        areasByAtomicNumber[AN] = 0.0;

                                    areasByAtomicNumber[AN] += segments.SegmentTypeAreas[i_segment][0];
                                }

                                for (auto& it : areasByAtomicNumber) {
                                    double this_atom_dGsolv_tau = abs(param.dGsolv_tau[it.first]);
                                    if (this_atom_dGsolv_tau == 0.0) {
                                        this_atom_dGsolv_tau = approximate_dGsolv_tau;
                                        atomicNumbersWithout_dGsolv_tau.push_back(it.first);
                                    }

                                    E_vdw += this_atom_dGsolv_tau * it.second;
                                }
                                double referenceStateCorrection = RT_kcalPerMol * log(molar_volume_ideal_gas / (calculations[calculationIndex].components[i_solvent_component]->molarVolumeAt25C / 1E6));

                                double E_diel = (calculations[calculationIndex].components[i_component]->epsilonInfinityTotalEnergy - param.dGsolv_E_gas[i_component]) * kcalPerMol_per_Hartree;
                                double mu_liquid = RT_kcalPerMol * calculations[calculationIndex].lnGammaTotal(i_concentration, i_component);
                                double E_ring = param.dGsolv_omega_ring * param.dGsolv_numberOfAtomsInRing[i_component];
                                dGsolv = E_diel + mu_liquid - E_vdw - E_ring - referenceStateCorrection - param.dGsolv_eta;
                            }

                            calculations[calculationIndex].dGsolv(i_concentration, i_component) = float(dGsolv);
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
                        warnings.push_back(" - For the following atomic numbers not all parameters are available for the calculation of solvation energies, an estimate was used: " + ANs);
                    }

                }
            }
            });
    }
    e.rethrow();
}

