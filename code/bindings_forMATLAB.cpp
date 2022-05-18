/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @authors: Simon Mueller & Andrés González de Castilla, 2022
*/



//#define MEASURE_TIME
//#define DEBUG_INFO

// this can be uncommented to get Intellisense
#include "C:\Program Files\MATLAB\R2019a\extern\include\mex.hpp"
#include "C:\Program Files\MATLAB\R2019a\extern\include\mexAdapter.hpp"

#include "general.hpp"
#include "core_functions.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
private:
	std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
	int mode;


	// Helper function to print output string on MATLAB command prompt.
	void displayOnMATLAB(std::string message) {
		ArrayFactory factory;
		matlabPtr->feval(u"fprintf", 0, std::vector<Array>
			({ factory.createScalar(message) }));
	}

	void displayTimeOnMATLAB(std::string message, unsigned long durationInMicroseconds) {
		displayOnMATLAB(message + ": " + std::to_string(durationInMicroseconds) + " microseconds\n");
	}

	template<typename T>
	T getNumericFieldValue(StructArray sArray, int index, std::string fieldName) {
		TypedArray<double> ta = sArray[index][fieldName];

		double val = ta[0];

		if (std::isnan(val)) {
			throw std::runtime_error("While reading a structure a NaN value was found.");
		}

		return (T)val;
	}

public:

    MexFunction()
    {
        matlabPtr = getEngine();
        n_ex = 0;
		
		display = [&](std::string message) { displayOnMATLAB(message); };
		displayTime = [&](std::string message, unsigned long durationInMicroseconds) { displayTimeOnMATLAB(message, durationInMicroseconds); };

        initialize(param);
    }
    
	void operator()(ArgumentList outputs, ArgumentList inputs) {
		try {
			calculateOnMatlab(outputs, inputs);
		}
		catch (const std::exception& e) {
			// this is not so nice as all exceptions are of the same type
			// but it serves the purpose and relays the error to MATLAB

			ArrayFactory factory;
			matlabPtr->feval(u"error", 0, std::vector<Array>
				({ factory.createScalar(e.what()) }));
		}
	}

	void calculateOnMatlab(ArgumentList outputs, ArgumentList inputs) {

#ifdef MEASURE_TIME
		startCalculationMeasurement();
#endif

		bool reloadConcentrations = false;
		bool reloadReferenceConcentrations = false;

		// check arguments
		n_ex += 1;
		if (inputs.size() > 0) {
			if (inputs[0].getType() != ArrayType::STRUCT) {
				throw std::runtime_error("Argument 1 is of wrong type: struct expected.");
			}
		}

		if (inputs.size() > 1) {
			if (inputs[1].getType() != ArrayType::STRUCT) {
				throw std::runtime_error("Argument 2 is of wrong type: struct expected.");
			}
		}

		if (inputs.size() == 3) {

			if (inputs[2].getType() == ArrayType::CELL && n_ex == 1) {
				mode = 0; // load components
			}
			else if (inputs[2].getType() == ArrayType::STRUCT && n_ex == 2) {
				mode = 1; // load calculations
			}
			else {
				if (inputs[2].getType() == ArrayType::STRUCT && n_ex > 2) {
					mode = 2; // perform calculations
					
					if (outputs.size() != 1) {
						throw std::runtime_error("Wrong number of outputs specified: 1 expected.");
					}
				}
				else {
					throw std::runtime_error("Argument 3 is of wrong type: cellarray[string] or array[struct] expected.");
				}
			}
		}
		else {
			if (inputs.size() == 5 && inputs[2].getType() == ArrayType::STRUCT && n_ex > 2 \
				&& inputs[3].getType() == ArrayType::LOGICAL && inputs[4].getType() == ArrayType::LOGICAL) {

				TypedArray<bool> temparray = TypedArray<bool>(inputs[3]);
				reloadConcentrations = temparray[0];

				temparray = TypedArray<bool>(inputs[4]);
				reloadReferenceConcentrations = temparray[0];
			}
			else {
				throw std::runtime_error("Wrong size of inputs.");
			}
		}

		
		// read the options
		StructArray const matlabStructArrayOpt = inputs[0];
		// load options
		param.sw_combTerm = getNumericFieldValue<int>(matlabStructArrayOpt, 0, "sw_SR_combTerm");
		param.sw_alwaysCalculateSizeRelatedParameters = getNumericFieldValue<int>(matlabStructArrayOpt, 0, "sw_SR_alwaysCalculateSizeRelatedParameters");
		param.sw_alwaysReloadSigmaProfiles = getNumericFieldValue<int>(matlabStructArrayOpt, 0, "sw_SR_alwaysReloadSigmaProfiles");
		param.sw_useSegmentReferenceStateForInteractionMatrix = getNumericFieldValue<int>(matlabStructArrayOpt, 0, "sw_SR_useSegmentReferenceStateForInteractionMatrix");
		param.sw_calculateContactStatisticsAndAdditionalProperties = getNumericFieldValue<int>(matlabStructArrayOpt, 0, "sw_SR_calculateContactStatisticsAndAdditionalProperties");
		param.sw_differentiateHydrogens = getNumericFieldValue<int>(matlabStructArrayOpt, 0, "sw_SR_differentiateHydrogens");
		TypedArray<MATLABString> sw_COSMOfiles_type = matlabStructArrayOpt[0]["sw_SR_COSMOfiles_type"];
		param.sw_COSMOfiles_type = sw_COSMOfiles_type[0];

		
		if (mode == 1) {
			if (param.sw_calculateContactStatisticsAndAdditionalProperties != 0) {
				CellArray partialInteractionMatrices = matlabStructArrayOpt[0]["sw_SR_partialInteractionMatrices"]; //untested
				param.numberOfPartialInteractionMatrices = partialInteractionMatrices.getNumberOfElements(); //untested
			}
			else {
				param.numberOfPartialInteractionMatrices = 0;
			}
		}
		
		param.sw_atomicNumber = getNumericFieldValue<int>(matlabStructArrayOpt, 0, "sw_SR_atomicNumber");
		param.sw_misfit = getNumericFieldValue<int>(matlabStructArrayOpt, 0, "sw_SR_misfit");

		if (param.sw_misfit < 0 && param.sw_misfit > 2) {
			throw std::runtime_error("sw_SR_misfit should have one of the following values: [0, 1, 2].");
		}
		
		
		// load parameters
		StructArray const matlabStructArrayPar = inputs[1];
		loadParametersOnMATLAB(matlabStructArrayPar);
		if (mode == 0) { // load the molecules
			CellArray cellArrayComponents = inputs[2];
			loadMoleculesOnMATLAB(cellArrayComponents);

		}
		else if (mode == 1) { // load the calculations
			StructArray const matlabStructArrayCalc = inputs[2];
			loadCalculationsOnMATLAB(matlabStructArrayCalc);
		}
		else if (mode == 2) {
			if (n_ex < 3) {
				throw std::runtime_error("Before trying to run a calculation please first load the molecules, then load the calculations and then execute the calculations.");
			}

			StructArray matlabStructArrayCalc(inputs[2]);

			if (param.sw_alwaysReloadSigmaProfiles == 1 && n_ex > 3) {
				reloadAllMolecules();
			}
			if (param.sw_alwaysCalculateSizeRelatedParameters == 1 || (param.sw_alwaysCalculateSizeRelatedParameters == 0 && n_ex == 3)) {
				resizeMonoatomicCations(param, molecules);
			}
			const size_t numCalcs = matlabStructArrayCalc.getNumberOfElements();
			std::vector<int> calculationIndices(numCalcs);
			param.sw_reloadConcentrations = 0;
			for (int i = 0; i < numCalcs; i++) {
				calculationIndices[i] = getNumericFieldValue<double>(matlabStructArrayCalc, i, "index");;
				if (reloadConcentrations == true) {
					param.sw_reloadConcentrations = 1;
					TypedArray<double> matlabArrayConcentrations = matlabStructArrayCalc[i]["concentrations"];
					for (int h = 0; h < calculations[calculationIndices[i]].originalNumberOfCalculations; h++) {
						int j = calculations[calculationIndices[i]].actualConcentrationIndices[h];
						std::vector<float> rowConcentration;

						float tempSumOfConcentrations = 0;
						for (int k = 0; k < calculations[calculationIndices[i]].components.size(); k++) {
							float val = (float)(matlabArrayConcentrations[j][k]);
							tempSumOfConcentrations += val;
							calculations[calculationIndices[i]].concentrations[j][k] = val;
						}

						if (abs(1.0f - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
							throw std::runtime_error("A concentration does not add up to unity.");
						}
					}
				}
				if (reloadReferenceConcentrations == true) {
					throw std::runtime_error("reloadReferenceConcentrations has not been implemented yet.");
				}
			}

			calculate(calculationIndices);
			
			// copy results to matlab matrices
			for (int i = 0; i < calculationIndices.size(); i++) {

				int calculationIndex = calculationIndices[i];

				TypedArrayRef<float> matlabArraylnGammaCombinatorial = std::move(matlabStructArrayCalc[calculationIndex]["ln_gamma_x_SR_combinatorial_calc"]);
				TypedArrayRef<float> matlabArraylnGammaResidual = std::move(matlabStructArrayCalc[calculationIndex]["ln_gamma_x_SR_residual_calc"]);
				TypedArrayRef<float> matlabArraylnGammaTotal = std::move(matlabStructArrayCalc[calculationIndex]["ln_gamma_x_SR_calc"]);

				for (int j = 0; j < calculations[calculationIndex].originalNumberOfCalculations; j++) {

					for (int k = 0; k < calculations[calculationIndex].components.size(); k++) {
						matlabArraylnGammaCombinatorial[j][k] = calculations[calculationIndex].lnGammaCombinatorial(j, k);
						matlabArraylnGammaResidual[j][k] = calculations[calculationIndex].lnGammaResidual(j, k);
						matlabArraylnGammaTotal[j][k] = calculations[calculationIndex].lnGammaTotal(j, k);
					}
				}

				if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {
					TypedArrayRef<float> matlabArrayContactStatistics = std::move(matlabStructArrayCalc[calculationIndex]["contact_statistics"]);
					TypedArrayRef<float> matlabArrayAverageSurfaceEnergies = std::move(matlabStructArrayCalc[calculationIndex]["average_surface_energies"]);

					for (int j = 0; j < calculations[calculationIndex].originalNumberOfCalculations; j++) {

						for (int k = 0; k < calculations[calculationIndex].components.size(); k++) {
							for (int l = 0; l < calculations[calculationIndex].components.size(); l++) {
								matlabArrayContactStatistics[j][k][l] = calculations[calculationIndex].contactStatistics(j, k, l);

								for (int p = 0; p < param.numberOfPartialInteractionMatrices + 1; p++) {
									matlabArrayAverageSurfaceEnergies[j][p][k][l] = calculations[calculationIndex].averageSurfaceEnergies(j, p, k, l);
								}
							}
						}
					}

					if (param.sw_calculateContactStatisticsAndAdditionalProperties == 2) {
						TypedArrayRef<float> matlabArrayPartialMolarEnergies = std::move(matlabStructArrayCalc[calculationIndex]["partial_molar_energies"]);
						for (int j = 0; j < calculations[calculationIndex].originalNumberOfCalculations; j++) {
							for (int p = 0; p < param.numberOfPartialInteractionMatrices + 1; p++) {
								for (int k = 0; k < calculations[calculationIndex].components.size(); k++) {
									matlabArrayPartialMolarEnergies[j][p][k] = calculations[calculationIndex].partialMolarEnergies(j, p, k);
								}
							}
						}		
					}
				}
			}
			outputs[0] = matlabStructArrayCalc;

#ifdef MEASURE_TIME
			stopCalculationMeasurement();
#endif
		}
	}

	void loadParametersOnMATLAB(StructArray const matlabStructArrayPar) {

		param.Aeff = getNumericFieldValue<double>(matlabStructArrayPar, 0, "Aeff");
		param.alpha = exp(getNumericFieldValue<double>(matlabStructArrayPar, 0, "ln_alpha"));
		
		param.CHB = exp(getNumericFieldValue<double>(matlabStructArrayPar, 0, "ln_CHB"));
		param.CHBT = getNumericFieldValue<double>(matlabStructArrayPar, 0, "CHBT");
		param.SigmaHB = getNumericFieldValue<double>(matlabStructArrayPar, 0, "SigmaHB");

		param.Rav = getNumericFieldValue<double>(matlabStructArrayPar, 0, "Rav");
		param.RavCorr = getNumericFieldValue<double>(matlabStructArrayPar, 0, "RavCorr");
		param.fCorr = getNumericFieldValue<double>(matlabStructArrayPar, 0, "fCorr");

		param.comb_SG_A_std = getNumericFieldValue<double>(matlabStructArrayPar, 0, "comb_SG_A_std");
		param.comb_SG_z_coord = getNumericFieldValue<double>(matlabStructArrayPar, 0, "comb_SG_z_coord");
		param.comb_modSG_exp = getNumericFieldValue<double>(matlabStructArrayPar, 0, "comb_modSG_exp");


		param.comb_lambda0 = getNumericFieldValue<double>(matlabStructArrayPar, 0, "comb_lambda0");
		param.comb_lambda1 = getNumericFieldValue<double>(matlabStructArrayPar, 0, "comb_lambda1");
		param.comb_lambda2 = getNumericFieldValue<double>(matlabStructArrayPar, 0, "comb_lambda2");

		param.comb_SGG_lambda = getNumericFieldValue<double>(matlabStructArrayPar, 0, "comb_SGG_lambda");
		param.comb_SGG_beta = getNumericFieldValue<double>(matlabStructArrayPar, 0, "comb_SGG_beta");

		TypedArray<double> matlabArrayRadii = matlabStructArrayPar[0]["radii"];
		for (int i = 0; i < param.R_i.size(); i++) {

			double val = matlabArrayRadii[i];

			if (!std::isnan(val)) {
				param.R_i[i + 1] = val; // change from zero based index to AN
			}
		}

		StructArray const matlabStructArrayParExp = matlabStructArrayPar[0]["exp"];
		auto fieldsExp = matlabStructArrayParExp.getFieldNames();
		std::vector<std::string> fieldNamesExp(fieldsExp.begin(), fieldsExp.end());

		for (int i = 0; i < fieldNamesExp.size(); i++) {
			param.exp_param[fieldNamesExp[i]] = getNumericFieldValue<double>(matlabStructArrayParExp, 0, fieldNamesExp[i]);
		}
	}

	void loadMoleculesOnMATLAB(CellArray cellArrayComponents) {

		if (n_ex != 1) {
			throw std::runtime_error("loadMolecules should only be executed once.");
		}

		const size_t numComponents = cellArrayComponents.getNumberOfElements();

		for (int i = 0; i < numComponents; i++) {
			TypedArray<MATLABString> cPath = cellArrayComponents[i];
			std::string componentPath = cPath[0];

			molecule newMolecule = loadNewMolecule(param, componentPath);

			molecules.push_back(std::move(std::make_shared<molecule>(newMolecule)));
		}
	}

	void loadCalculationsOnMATLAB(StructArray const matlabStructArrayCalc) {

		if (n_ex != 2) {
			throw std::runtime_error("The loading of the calculations can only be executed once.");
		}

		// load the calculations
		// array of structs = calculations(i)
		// calculations(i).index = 0
		// calculations(i).components = [1 2];
		// calculations(i).temperatures = [298.15, 298.15, 298.15];
		// calculations(i).concentrations = [0.2 0.8; 0.5 0.5; 0.1 0.9];
		// calculations(i).referenceStateTypes = [0 2 3];  0: PureComponents | 1: PureComponentsOnlyNeutral | 2: ReferenceMixture (also for InfiniteDilution) | 3: SegmentReference
		// calculations(i).referenceStateConcentrations = {[] [] [0 1]};  # cell array of concentrations where referenceStateType == 3 is an array

		const size_t numCalcs = matlabStructArrayCalc.getNumberOfElements();

		// important as otherwise the data behind Eigen::Maps
		// is lost when the array is resized.
		calculations.reserve(numCalcs);

		if (numCalcs == 0) {
			throw std::runtime_error("Please specify at least one calculation.");
		}

		for (int i = 0; i < numCalcs; i++) {
			
			// array of component indices
			TypedArray<double> matlabArrayComponents = matlabStructArrayCalc[i]["components"];
			int numberOfComponents = int(matlabArrayComponents.getNumberOfElements());
			calculation newCalculation(numberOfComponents);
			for (int j = 0; j < numberOfComponents; j++) {
				std::shared_ptr<molecule> thisMolecule = molecules[int(matlabArrayComponents[j]) - 1]; // minus one because matlab indexing is one-based
				for (int k = 0; k < thisMolecule->segments.size(); k++) {
					newCalculation.segments.add(j, thisMolecule->segments.SegmentTypeGroup[k],
						thisMolecule->segments.SegmentTypeSigma[k],
						thisMolecule->segments.SegmentTypeSigmaCorr[k],
						thisMolecule->segments.SegmentTypeHBtype[k],
						thisMolecule->segments.SegmentTypeAtomicNumber[k],
						thisMolecule->segments.SegmentTypeAreas[k][0]);
				}
				newCalculation.components.push_back(thisMolecule);
			}
			newCalculation.segments.sort();
			newCalculation.segments.shrink_to_fit();
			// temperature
			TypedArray<double> matlabArrayTemperatures = matlabStructArrayCalc[i]["temperatures"];
			// array of compositions
			TypedArray<double> matlabArrayConcentrations = matlabStructArrayCalc[i]["concentrations"];
			ArrayDimensions ad = matlabArrayConcentrations.getDimensions();
			TypedArray<double> matlabArrayReferenceStateTypes = matlabStructArrayCalc[i]["reference_state_types"];
			CellArray matlabCellArrayReferenceStateConcentrations = matlabStructArrayCalc[i]["reference_state_concentrations"];
			const size_t numReferenceStateConcentrations = matlabCellArrayReferenceStateConcentrations.getNumberOfElements();

			// concentrations
			for (int j = 0; j < ad[0]; j++) {
				std::vector<float> rowConcentration;

				float tempSumOfConcentrations = 0;
				for (int k = 0; k < newCalculation.components.size(); k++) {
					float val = (float)(matlabArrayConcentrations[j][k]);
					tempSumOfConcentrations += val;
					rowConcentration.push_back(val);
				}

				if (abs(1.0f - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
					throw std::runtime_error("For calculation number " + std::to_string(i) + ", the concentrations do not add up to unity. residual concentration: " + std::to_string(abs(1.0f - tempSumOfConcentrations)));
				}

				float temperature = (float)matlabArrayTemperatures[j];

				newCalculation.temperatures.push_back(temperature);
				newCalculation.concentrations.push_back(rowConcentration);
			}

			if (numReferenceStateConcentrations != newCalculation.concentrations.size()) {
				throw std::runtime_error("concentrations and referenceStateConcentrations of calculation number " + std::to_string(i) + " have different sizes.\n");
			}

			newCalculation.originalNumberOfCalculations = (unsigned short)newCalculation.concentrations.size();

			// reference states

			for (int j = 0; j < ad[0]; j++) {

				TypedArray<double> matlabArrayReferenceStateConcentrations = matlabCellArrayReferenceStateConcentrations[j];

				int referenceStateType = int(matlabArrayReferenceStateTypes[j]);
				newCalculation.referenceStateType.push_back((unsigned short)referenceStateType);

				float tempSumOfConcentrations = 0;
				for (int k = 0; k < (size_t)matlabArrayReferenceStateConcentrations.getNumberOfElements(); k++) {
					tempSumOfConcentrations += (float)matlabArrayReferenceStateConcentrations[k];
				}

				if (referenceStateType == 0) { // PureComponents: untested
					if (matlabArrayReferenceStateConcentrations.getNumberOfElements() != 0) {
						throw std::runtime_error("A reference state concentration was specified for a calculation with reference state PureComponents, this does not make sense.");
					}

					std::vector<int> thisReferenceStateCalculationIndices;
					for (int k = 0; k < newCalculation.components.size(); k++) {

						std::vector<float> referenceStateConcentration;
						for (int m = 0; m < newCalculation.components.size(); m++) {
							referenceStateConcentration.push_back(k == m ? 1.0f : 0.0f);
						}

						float temperature = newCalculation.temperatures[j];
						int referenceStateCalculationIndex = newCalculation.addOrFindArrayIndexForConcentration(referenceStateConcentration, temperature);

						thisReferenceStateCalculationIndices.push_back(referenceStateCalculationIndex);

					}
					newCalculation.referenceStateCalculationIndices.push_back(thisReferenceStateCalculationIndices);
				}
				else if (referenceStateType == 1) { // PureComponentsOnlyNeutral: untested
					if (matlabArrayReferenceStateConcentrations.getNumberOfElements() != 0) {
						throw std::runtime_error("A reference state concentration was specified for a calculation with reference state PureComponentsOnlyNeutral, this does not make sense.");
					}

					std::vector<int> thisReferenceStateCalculationIndices;
					for (int k = 0; k < newCalculation.components.size(); k++) {

						if (newCalculation.components[k]->moleculeCharge == 0) {

							std::vector<float> referenceStateConcentration;
							for (int m = 0; m < newCalculation.components.size(); m++) {
								referenceStateConcentration.push_back(k == m ? 1.0f : 0.0f);
							}

							float temperature = newCalculation.temperatures[0];
							int referenceStateCalculationIndex = newCalculation.addOrFindArrayIndexForConcentration(referenceStateConcentration, temperature);

							thisReferenceStateCalculationIndices.push_back(referenceStateCalculationIndex);
						}
						else {
							thisReferenceStateCalculationIndices.push_back(-1);
						}

					}
					newCalculation.referenceStateCalculationIndices.push_back(thisReferenceStateCalculationIndices);
				}
				else if (referenceStateType == 2) { // ReferenceMixture
					if (matlabArrayReferenceStateConcentrations.getNumberOfElements() != newCalculation.components.size()) {
						std::runtime_error("A reference state concentration was specified with the wrong amount of concentrations.");
					}

					std::vector<float> referenceStateConcentration;
					for (int m = 0; m < newCalculation.components.size(); m++) {
						referenceStateConcentration.push_back((float)matlabArrayReferenceStateConcentrations[m]);
					}

					if (abs(1.0f - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
						throw std::runtime_error("For calculation number " + std::to_string(i) + ", the reference concentrations do not add up to unity. residual concentration: " + std::to_string(abs(1.0f - tempSumOfConcentrations)));
					}

					float temperature = newCalculation.temperatures[j];
					int referenceStateCalculationIndex = newCalculation.addOrFindArrayIndexForConcentration(referenceStateConcentration, temperature);

					std::vector<int> thisReferenceStateCalculationIndices;
					for (int m = 0; m < newCalculation.components.size(); m++) {
						thisReferenceStateCalculationIndices.push_back(referenceStateCalculationIndex);
					}
					newCalculation.referenceStateCalculationIndices.push_back(thisReferenceStateCalculationIndices);
				}
				else if (referenceStateType == 3) { // COSMO

					if (matlabArrayReferenceStateConcentrations.getNumberOfElements() != 0) {
						std::runtime_error("A reference state concentration was specified for a calculation with reference state COSMO, this does not make sense.");
					}

					std::vector<int> thisReferenceStateCalculationIndices(numberOfComponents, -1);
					newCalculation.referenceStateCalculationIndices.push_back(thisReferenceStateCalculationIndices);
				}
				else {
					std::runtime_error("An unknown reference state type was given.");
				}

			}

			// I tried lots of things to directly bind to MATLAB matrices like was done for python
			// to avoid the step of copying the values explicitly, I however did not get it to work.

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
			finishCalculationInitiation(newCalculation, param);
			// the following std::move statement is very important as otherwise the Eigen::Map
			// are pointing to matrices deleted after exiting this function
			calculations.push_back(std::move(newCalculation));
		}
	} 
};