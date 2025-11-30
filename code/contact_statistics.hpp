/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/

#pragma once

#include "types.hpp"

// this code is not optimized for fast calculation yet
void calculateContactStatistics(calculation& _calculation, Eigen::MatrixXf& A_int_float, std::vector<Eigen::MatrixXd>& partialInteractionMatrices, Eigen::MatrixXf& Tau_float, float* Gamma, int i_concentration, Eigen::Tensor<float, 4, Eigen::RowMajor>& temporary_averageInteractionEnergies, Eigen::Tensor<float, 3, Eigen::RowMajor>& temporary_partialMolarEnergies, parameters& param) {

	/* Calculate contact statistics for all compositions.*/
	const size_t numberOfComponents = _calculation.components.size();
	const size_t numberOfSegments = _calculation.segments.size();

	float* x = &(_calculation.concentrations[i_concentration][0]);

	Eigen::VectorXd X_vector = _calculation.segmentConcentrations(Eigen::seqN(0, numberOfSegments), i_concentration).cast<double>();
	Eigen::VectorXd x_vector = Eigen::Map<Eigen::VectorXf>(x, numberOfComponents).cast<double>();
	Eigen::VectorXd gamma_vector = Eigen::Map<Eigen::VectorXf>(Gamma, numberOfSegments).cast<double>();

	Eigen::ArrayXXd Tau = Tau_float(Eigen::seqN(0, numberOfSegments), Eigen::seqN(0, numberOfSegments)).cast<double>();
	Eigen::ArrayXXd A_int = A_int_float(Eigen::seqN(0, numberOfSegments), Eigen::seqN(0, numberOfSegments)).cast<double>();

	double div_Aeff = 1 / param.Aeff;

	Eigen::MatrixXd N_i_I(numberOfComponents, numberOfSegments); // number of segments of type I on molecule i

	for (int j = 0; j < numberOfComponents; j++) {
		for (int k = 0; k < numberOfSegments; k++) {
			N_i_I(j, k) = _calculation.segments.SegmentTypeAreas[k][j] * div_Aeff;
		}
	}

	Eigen::VectorXd N_mol = N_i_I.rowwise().sum();
	double sum_xN = (N_mol.array() * x_vector.array()).sum(); // total number of segments per total number of molecules in mixture

	// set some temporary variables
	Eigen::ArrayXXd gammai_gammajT = gamma_vector * gamma_vector.transpose();
	Eigen::ArrayXXd gammai_gammajT_Tau = gammai_gammajT * Tau;
	Eigen::ArrayXXd X_gammai_gammajT_Tau = gammai_gammajT_Tau.rowwise() * X_vector.transpose().array();

	// molecular contact probabilities and average energies

	for (int mol_i = 0; mol_i < numberOfComponents; mol_i++) {
		for (int mol_j = 0; mol_j < numberOfComponents; mol_j++) {

			double factor = x[mol_j] / sum_xN;

			Eigen::ArrayXXd temp = (gamma_vector.array() * N_i_I(mol_i, Eigen::indexing::all).transpose().array()).matrix() * (gamma_vector.array() * N_i_I(mol_j, Eigen::indexing::all).transpose().array()).matrix().transpose();
			Eigen::ArrayXXd weighting = temp * Tau;

			// contact statistics are directly written to the output
			if (i_concentration < _calculation.originalNumberOfCalculations) {
				_calculation.contactStatistics(i_concentration, mol_i, mol_j) = float(weighting.sum() * factor / N_mol(mol_i));
			}
			temporary_averageInteractionEnergies(i_concentration, 0, mol_i, mol_j) = float(N_AVOGADRO * 0.5 * factor * (weighting * A_int.array()).sum());

			for (int k = 0; k < param.numberOfPartialInteractionMatrices; k++) {
				for (int mol_j = 0; mol_j < numberOfComponents; mol_j++) {
					temporary_averageInteractionEnergies(i_concentration, k + 1, mol_i, mol_j) = float(N_AVOGADRO * 0.5 * factor * (weighting * partialInteractionMatrices[k].array()).sum());
				}
			}
		}
	}

	// Calculate partial molar properties
	if (param.sw_calculateContactStatisticsAndAdditionalProperties == 2) {

		// Calculate gradients of segment activity coefficients towards mole fraction and store in newly allocated array
		// Fill RHSs of LES for the derivative towards each molecule if molecular mole fraction is unequal unity
		Eigen::MatrixXd b(numberOfSegments, numberOfComponents);
		for (int mol_i = 0; mol_i < numberOfComponents; mol_i++) {

			Eigen::ArrayXXd temp = gamma_vector.array() * (Tau.matrix() * (N_i_I(mol_i, Eigen::indexing::all).transpose().array() * gamma_vector.array()).matrix()).array();
			b(Eigen::indexing::all, mol_i) = (gamma_vector.array() / sum_xN) * (N_mol[mol_i] - temp);
		}

		// Fill array of coefficients for the LES
		Eigen::MatrixXd A = (((gamma_vector.array() * gamma_vector.array()).matrix() * X_vector.transpose()).array() * Tau).matrix() + Eigen::MatrixXd::Identity(numberOfSegments, numberOfSegments);

		// calculate decomposition and solve for several RHSs
		Eigen::HouseholderQR<Eigen::MatrixXd> decomposition = A.householderQr();
		Eigen::MatrixXd ndGamma_dni(numberOfSegments, numberOfComponents);

#if defined(DEBUG_INFO)
		double relativeError = 0;
#endif
		for (int mol_i = 0; mol_i < numberOfComponents; mol_i++) {
			ndGamma_dni(Eigen::indexing::all, mol_i) = decomposition.solve(b(Eigen::indexing::all, mol_i));

#if defined(DEBUG_INFO)
			relativeError = std::max(relativeError, (A * ndGamma_dni(Eigen::indexing::all, 0) - b(Eigen::indexing::all, 0)).norm() / b(Eigen::indexing::all, 0).norm());
#endif
		}

#if defined(DEBUG_INFO)
		if (relativeError > 10e-14) {
			throw std::runtime_error("When calculating ndGamma_dni the equation system was not solved correctly.");
		}
#endif

		Eigen::ArrayXXd ndX_dni = (N_i_I.transpose() - X_vector * N_mol.transpose()) / sum_xN;
		Eigen::ArrayXXd temp1 = Tau.rowwise() * (X_vector.array() * gamma_vector.array()).transpose();
		Eigen::ArrayXXd temp2 = Tau * (gamma_vector * X_vector.transpose()).array();

		for (int mol_i = 0; mol_i < numberOfComponents; mol_i++) {

			Eigen::ArrayXXd ndpss_dni_mol_i = gammai_gammajT_Tau.rowwise() * ndX_dni(Eigen::indexing::all, mol_i).transpose();
			ndpss_dni_mol_i += temp1.colwise() * ndGamma_dni(Eigen::indexing::all, mol_i).array();
			ndpss_dni_mol_i += temp2.rowwise() * ndGamma_dni(Eigen::indexing::all, mol_i).array().transpose();

			Eigen::ArrayXXd temp3 = sum_xN * ((ndpss_dni_mol_i.colwise() * X_vector.array()) + (X_gammai_gammajT_Tau.colwise() * ndX_dni(Eigen::indexing::all, mol_i)));
			Eigen::ArrayXXd weighting = (0.5 * (X_gammai_gammajT_Tau.colwise() * X_vector.array() * N_mol(mol_i) + temp3));

			temporary_partialMolarEnergies(i_concentration, 0, mol_i) = float((weighting * (A_int.cast<double>()).array()).sum());

			for (int k = 0; k < param.numberOfPartialInteractionMatrices; k++) {
				temporary_partialMolarEnergies(i_concentration, k + 1, mol_i) = float((weighting * partialInteractionMatrices[k].array()).sum());
			}
		}

	}

}