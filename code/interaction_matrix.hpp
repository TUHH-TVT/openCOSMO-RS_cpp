/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/


#pragma once

#include "helper_functions.hpp"
#include <stdexcept>
#include <immintrin.h>// always include this as at least SSE3 is required
#include <cmath>
#include <map>
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

    // hb
    double CHB_T = 0;
    double buffdb1 = 1.0 - param.CHBT + param.CHBT * (298.15 / (temperature));
    if (buffdb1 > 0) CHB_T = param.CHB * buffdb1;

    double const hb_prefactor = param.Aeff * CHB_T;
    double val = 0;
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

                val =misfit_prefactor * sigmaMFij * (sigmaMFij + param.fCorr * (sigmaTransi + sigmaTransj));

                if (segments.SegmentTypeAtomicNumber[i] == 9 || segments.SegmentTypeAtomicNumber[j] == 9 ) {
                    val = misfit_prefactor * sigmaMFij * (sigmaMFij + param.fCorr * (sigmaTransi + sigmaTransj)) + param.E_F_corr;
                }
            }
            else {

                val = misfit_prefactor * sigmaMFij * sigmaMFij;
            }


            if (sigmai < minus_sigmaHB && sigmaj > sigmaHB) {

                if (segments.SegmentTypeHBtype[i] == 1 && segments.SegmentTypeHBtype[j] == 2) {
                    val += hb_prefactor * (sigmaj - sigmaHB) * (sigmai + sigmaHB);
                }
                else if (segments.SegmentTypeHBtype[i] == 2 && segments.SegmentTypeHBtype[j] == 1) {
                    throw std::runtime_error("This should not happen. Wrong assumption calculating interaction matrix. 1");
                }

            }
            else if (sigmaj < minus_sigmaHB && sigmai > sigmaHB) {

                if (segments.SegmentTypeHBtype[i] == 2 && segments.SegmentTypeHBtype[j] == 1) {
                    val += hb_prefactor * (sigmai - sigmaHB) * (sigmaj + sigmaHB);
                }
                else if (segments.SegmentTypeHBtype[i] == 1 && segments.SegmentTypeHBtype[j] == 2) {
                    throw std::runtime_error("This should not happen. Wrong assumption calculating interaction matrix. 2");
                }

            }
            vdw_val = 0.0;

           
            if (param.sw_SR_polarizabilities > 0) {

                if (param.sw_SR_polarizabilities == 7) {
                    ionizationPotentialCorrection = 1/ ioniz_potential_ev[segments.SegmentTypeAtomicNumber[i]]/ ioniz_potential_ev[segments.SegmentTypeAtomicNumber[j]];
                    vdw_val = param.Aeff * param.m_vdW * sqrt(alphai * alphaj * ionizationPotentialCorrection);
                }
                else if (param.sw_SR_polarizabilities == 8) {
                    ionizationPotentialCorrection = 1/(ioniz_potential_ev[segments.SegmentTypeAtomicNumber[i]] + ioniz_potential_ev[segments.SegmentTypeAtomicNumber[j]]);
                    vdw_val = param.Aeff * param.m_vdW * sqrt(alphai * alphaj) * ionizationPotentialCorrection;
                }
                else {
                    vdw_val = param.Aeff * param.m_vdW * sqrt(alphai * alphaj);
                }
            }

            // always j >= i
            val -= vdw_val;
            A_int(j, i) = static_cast<calcType>(val);
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
