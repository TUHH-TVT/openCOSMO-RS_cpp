/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/


#pragma once

#include "helper_functions.hpp"
#include <stdexcept>
#include <cmath>
#include <map>

// Partial interaction matrices that can be requested via sw_SR_partialInteractionMatrices.
// Their contact-weighted values are written as additional entries to average_surface_energies
// and partial_molar_energies (index 0 is always A_int, index k + 1 the k-th requested matrix).
//   E_mf  : misfit energy, including the fluorine correction E_F_corr if polarizabilities are used
//   G_hb  : hydrogen bonding free energy (temperature dependent, as used in A_int)
//   E_hb  : hydrogen bonding enthalpy, from G_hb via Gibbs-Helmholtz
//   E_vdW : dispersion energy from the polarizabilities
// A_int = E_mf + G_hb + E_vdW
void setPartialInteractionMatrices(parameters& param, const std::vector<std::string>& names) {

    param.partialInteractionMatrixIndex_E_mf = -1;
    param.partialInteractionMatrixIndex_G_hb = -1;
    param.partialInteractionMatrixIndex_E_hb = -1;
    param.partialInteractionMatrixIndex_E_vdW = -1;

    for (int i = 0; i < int(names.size()); i++) {
        int* index = nullptr;
        if (names[i] == "E_mf") index = &param.partialInteractionMatrixIndex_E_mf;
        else if (names[i] == "G_hb") index = &param.partialInteractionMatrixIndex_G_hb;
        else if (names[i] == "E_hb") index = &param.partialInteractionMatrixIndex_E_hb;
        else if (names[i] == "E_vdW") index = &param.partialInteractionMatrixIndex_E_vdW;
        else throw std::runtime_error("Unknown partial interaction matrix '" + names[i] + "'. Possible values are: E_mf, G_hb, E_hb, E_vdW.");

        if (*index != -1) throw std::runtime_error("The partial interaction matrix '" + names[i] + "' was requested more than once.");
        *index = i;
    }

    param.numberOfPartialInteractionMatrices = int(names.size());
}

// returns lower left triangular matrices, this is because this way the matrix is accessed in sequential order in memory
void calculateInteractionMatrix(segmentTypeCollection& segments, MatrixCalcType& A_int, std::vector<Eigen::MatrixXd>& partialInteractionMatrices, parameters& param, double temperature) {

    const int numberOfSegments = int(segments.size());

    double sigmai = 0.0;
    double sigmaj = 0.0;
    double alphai = 0.0;
    double alphaj = 0.0;
    double sigmaCorri = 0.0;
    double sigmaCorrj = 0.0;
    double sigmaTransi = 0.0;
    double sigmaTransj = 0.0;
    double sigmaMFij = 0.0;

    double const sigmaHB = param.SigmaHB;
    double const minus_sigmaHB = -1 * sigmaHB;

    std::vector<double> ChargeRaster = param.ChargeRaster;

    // neutral - neutral interactions -------------------------------------------------------------------------------------------------
    // misfit 
    double const misfit_prefactor = param.Aeff * param.alpha * 0.5;

    double ionizationPotentialCorrection = 0;
    const bool hydrogensAreDifferentiated = param.sw_differentiateHydrogens == 1;

    // hb
    double CHB_T = 0;
    double buffdb1 = 1.0 - param.CHBT + param.CHBT * (298.15 / (temperature));
    if (buffdb1 > 0) CHB_T = param.CHB * buffdb1;

    double const hb_prefactor = param.Aeff * CHB_T;

    // factor to convert the hydrogen bonding free energy to the enthalpy: d(G/T)/d(1/T) / G
    double hb_enthalpyFactor = 0;
    if (buffdb1 > 0) hb_enthalpyFactor = (1.0 - param.CHBT + 2 * param.CHBT * (298.15 / temperature)) / buffdb1;

    const bool fillPartialInteractionMatrices = param.sw_calculateContactStatisticsAndAdditionalProperties > 0 && param.numberOfPartialInteractionMatrices > 0;

    double val = 0;
    double mf_val = 0;
    double hb_val = 0;
    double vdw_val = 0;
 
    int UpperBoundIndexForNeutralComponents = std::max(segments.upperBoundIndexForGroup[0], segments.upperBoundIndexForGroup[1]);
    UpperBoundIndexForNeutralComponents = std::max(UpperBoundIndexForNeutralComponents, segments.upperBoundIndexForGroup[2]);

    for (int i = segments.lowerBoundIndexForGroup[0]; i < UpperBoundIndexForNeutralComponents; i++) {

        sigmai = segments.SegmentTypeSigma[i];
        alphai = segments.SegmentTypeAtomicPolariz[i];

        if (param.sw_misfit > 0) {
            sigmaCorri = segments.SegmentTypeSigmaCorr[i];
            sigmaTransi = sigmaCorri - 0.816 * sigmai;
        }

        for (int j = i; j < UpperBoundIndexForNeutralComponents; j++) {

            sigmaj = segments.SegmentTypeSigma[j];
            alphaj = segments.SegmentTypeAtomicPolariz[j];
            sigmaMFij = sigmai + sigmaj;

            if (param.sw_misfit > 0) {
                sigmaCorrj = segments.SegmentTypeSigmaCorr[j];
                sigmaTransj = sigmaCorrj - 0.816 * sigmaj;

                mf_val = misfit_prefactor * sigmaMFij * (sigmaMFij + param.fCorr * (sigmaTransi + sigmaTransj));
            }
            else {

                mf_val = misfit_prefactor * sigmaMFij * sigmaMFij;
            }

            hb_val = 0.0;

            if (sigmai < minus_sigmaHB && sigmaj > sigmaHB) {

                if (segments.SegmentTypeHBtype[i] == 1 && segments.SegmentTypeHBtype[j] == 2) {
                    hb_val = hb_prefactor * (sigmaj - sigmaHB) * (sigmai + sigmaHB);
                }
                else if (segments.SegmentTypeHBtype[i] == 2 && segments.SegmentTypeHBtype[j] == 1) {
                    throw std::runtime_error("This should not happen. Wrong assumption calculating interaction matrix. 1");
                }

            }
            else if (sigmaj < minus_sigmaHB && sigmai > sigmaHB) {

                if (segments.SegmentTypeHBtype[i] == 2 && segments.SegmentTypeHBtype[j] == 1) {
                    hb_val = hb_prefactor * (sigmai - sigmaHB) * (sigmaj + sigmaHB);
                }
                else if (segments.SegmentTypeHBtype[i] == 1 && segments.SegmentTypeHBtype[j] == 2) {
                    throw std::runtime_error("This should not happen. Wrong assumption calculating interaction matrix. 2");
                }

            }
            vdw_val = 0.0;

           
            if (param.sw_SR_polarizabilities > 0) {

                if (param.sw_SR_polarizabilities == 7) {
                    ionizationPotentialCorrection = 1 / getIonizationPotentialEv(segments.SegmentTypeAtomicNumber[i], hydrogensAreDifferentiated) / getIonizationPotentialEv(segments.SegmentTypeAtomicNumber[j], hydrogensAreDifferentiated);
                    vdw_val = param.Aeff * param.m_vdW * sqrt(alphai * alphaj * ionizationPotentialCorrection);
                }
                else if (param.sw_SR_polarizabilities == 8) {
                    ionizationPotentialCorrection = 1 / (getIonizationPotentialEv(segments.SegmentTypeAtomicNumber[i], hydrogensAreDifferentiated) + getIonizationPotentialEv(segments.SegmentTypeAtomicNumber[j], hydrogensAreDifferentiated));
                    vdw_val = param.Aeff * param.m_vdW * sqrt(alphai * alphaj) * ionizationPotentialCorrection;
                }
                else {
                    vdw_val = param.Aeff * param.m_vdW * sqrt(alphai * alphaj);
                }

                // the dispersion contribution is attractive
                vdw_val = -vdw_val;

                // repulsive fluorine correction of the misfit energy, Eq. (14) in Grigorash et al., Chem. Eng. Sci. 319 (2026) 122170,
                // part of the parameterization with polarizabilities
                if (segments.SegmentTypeAtomicNumber[i] == 9 || segments.SegmentTypeAtomicNumber[j] == 9) {
                    mf_val += param.E_F_corr;
                }
            }

            val = mf_val + hb_val + vdw_val;

            // always j >= i
            A_int(j, i) = static_cast<calcType>(val);

            if (fillPartialInteractionMatrices) {
                if (param.partialInteractionMatrixIndex_E_mf != -1) partialInteractionMatrices[param.partialInteractionMatrixIndex_E_mf](j, i) = mf_val;
                if (param.partialInteractionMatrixIndex_G_hb != -1) partialInteractionMatrices[param.partialInteractionMatrixIndex_G_hb](j, i) = hb_val;
                if (param.partialInteractionMatrixIndex_E_hb != -1) partialInteractionMatrices[param.partialInteractionMatrixIndex_E_hb](j, i) = hb_enthalpyFactor * hb_val;
                if (param.partialInteractionMatrixIndex_E_vdW != -1) partialInteractionMatrices[param.partialInteractionMatrixIndex_E_vdW](j, i) = vdw_val;
            }
        }
    }

    if (param.sw_useSegmentReferenceStateForInteractionMatrix == 1) {

        // apply the required reference state for the interaction matrix
        for (int i = 0; i < numberOfSegments; i++) {
            for (int j = i + 1; j < numberOfSegments; j++) {
                // always j >= i + 1
                A_int(j, i) = A_int(j, i) - 0.5 * (A_int(i, i) + A_int(j, j));
            }
        }

        // assign zero to the diagonal
        for (int i = 0; i < numberOfSegments; i++) {
            A_int(i, i) = 0.0f;
        }

        if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {
            for (int h = 0; h < param.numberOfPartialInteractionMatrices; h++) {
                for (int i = 0; i < numberOfSegments; i++) {
                    for (int j = i + 1; j < numberOfSegments; j++) {
                        // always j >= i + 1
                        partialInteractionMatrices[h](j, i) = partialInteractionMatrices[h](j, i) - 0.5 * (partialInteractionMatrices[h](i, i) + partialInteractionMatrices[h](j, j));
                    }
                }
            }

            for (int h = 0; h < param.numberOfPartialInteractionMatrices; h++) {
                for (int i = 0; i < numberOfSegments; i++) {
                    partialInteractionMatrices[h](i, i) = 0.0;
                }
            }
        }
    }
}
