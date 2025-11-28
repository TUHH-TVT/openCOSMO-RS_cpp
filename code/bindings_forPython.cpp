/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/


//#define MEASURE_TIME
//#define DEBUG_INFO

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
using namespace pybind11::literals;

#include "general.hpp"
#include "core_functions.hpp"

void displayOnPython(std::string message) {
	py::print(message, "end"_a = "");
}

void displayTimeOnPython(std::string message, unsigned long durationInMicroseconds) {
	displayOnPython(message + ": " + std::to_string(durationInMicroseconds) + " microseconds\n");
}

void initializeOnPython() {
	display = displayOnPython;
	displayTime = displayTimeOnPython;

    initialize(param);
}

void loadParametersOnPython(py::dict parameters) {

	param.Aeff = parameters["Aeff"].cast<double>();
	param.alpha = exp(parameters["ln_alpha"].cast<double>());

	param.CHB = exp(parameters["ln_CHB"].cast<double>());
	param.CHBT = parameters["CHBT"].cast<double>();
	param.SigmaHB = parameters["SigmaHB"].cast<double>();

	param.Rav = parameters["Rav"].cast<double>();

	param.E_F_corr = parameters["E_F_corr"].cast<double>();
	param.m_vdW = parameters["m_vdW"].cast<double>();

	if (param.sw_misfit > 0) {
		param.fCorr = parameters["fCorr"].cast<double>();
		param.RavCorr = parameters["RavCorr"].cast<double>();
	}

	if (param.sw_combTerm == 1 || param.sw_combTerm == 3) {
		param.comb_SG_A_std = parameters["comb_SG_A_std"].cast<double>();
		param.comb_SG_z_coord = parameters["comb_SG_z_coord"].cast<double>();

		if(param.sw_combTerm == 3)
			param.comb_modSG_exp = parameters["comb_modSG_exp"].cast<double>();
	}
	
	if (param.sw_combTerm == 2 || param.sw_combTerm == 5) {
		param.comb_lambda0 = parameters["comb_lambda0"].cast<double>();
		param.comb_lambda1 = parameters["comb_lambda1"].cast<double>();
		param.comb_lambda2 = parameters["comb_lambda2"].cast<double>();
	}

	if (param.sw_combTerm == 4) {
		param.comb_SGG_lambda = parameters["comb_SGG_lambda"].cast<double>();
		param.comb_SGG_beta = parameters["comb_SGG_beta"].cast<double>();
	}


	if (parameters.contains("dGsolv_eta")) {
		param.dGsolv_eta = parameters["dGsolv_eta"].cast<double>();
		param.dGsolv_omega_ring = parameters["dGsolv_omega_ring"].cast<double>();

		py::dict dGsolv_tau = parameters["dGsolv_tau"];

		for (auto item : dGsolv_tau) {
			std::string key = item.first.cast<std::string>();
			param.dGsolv_tau[std::stoi(key)] = item.second.cast<double>();
		}
	}

	if (parameters.contains("radii")) {
		py::dict radii = parameters["radii"];

		for (auto item : radii) {
			std::string key = item.first.cast<std::string>();
			param.R_i[std::stoi(key)] = item.second.cast<double>();
		}
	}
	// experimental parameters for prototyping
	if (parameters.contains("exp")) {
		py::dict exp = parameters["exp"];
		for (auto item : exp)
			param.exp_param[item.first.cast<std::string>()] = item.second.cast<double>();
	}

}

void loadMoleculesOnPython(py::dict options, py::dict parameters, py::list componentPaths) {

	// if uninitialized
	if (n_ex == -1) {
		initializeOnPython();
	}

	n_ex += 1;

	if (n_ex != 1) {
		throw std::runtime_error("loadMolecules should only be executed once after calling initiate.");
	}

	// options

	param.sw_combTerm = options["sw_SR_combTerm"].cast<int>();
	param.sw_alwaysCalculateSizeRelatedParameters = options["sw_SR_alwaysCalculateSizeRelatedParameters"].cast<int>();
	param.sw_alwaysReloadSigmaProfiles = options["sw_SR_alwaysReloadSigmaProfiles"].cast<int>();
	param.sw_useSegmentReferenceStateForInteractionMatrix = options["sw_SR_useSegmentReferenceStateForInteractionMatrix"].cast<int>();

	param.sw_calculateContactStatisticsAndAdditionalProperties = options["sw_SR_calculateContactStatisticsAndAdditionalProperties"].cast<int>();

	param.sw_differentiateHydrogens = options["sw_SR_differentiateHydrogens"].cast<int>();
	param.sw_differentiateMoleculeGroups = options["sw_SR_differentiateMoleculeGroups"].cast<int>();
	param.sw_COSMOfiles_type = options["sw_SR_COSMOfiles_type"].cast<std::string>();
	param.sw_SR_polarizabilities = options["sw_SR_polarizabilities"].cast<int>();

	if (param.sw_calculateContactStatisticsAndAdditionalProperties != 0) {
		py::list partialInteractionMatrices = options["sw_SR_partialInteractionMatrices"];
		param.numberOfPartialInteractionMatrices = int(partialInteractionMatrices.size());
	} else {
		param.numberOfPartialInteractionMatrices = 0;
	}

	param.sw_atomicNumber = options["sw_SR_atomicNumber"].cast<int>();
	param.sw_misfit = options["sw_SR_misfit"].cast<int>();

	if (param.sw_misfit < 0 && param.sw_misfit > 2) {
		throw std::runtime_error("sw_SR_misfit should have one of the following values: [0, 1, 2].");
	}
	param.sw_skip_COSMOSPACE_errors = options["sw_skip_COSMOSPACE_errors"].cast<int>();

	// parameters
	loadParametersOnPython(parameters);

	for (auto componentPath : componentPaths) {
		molecule newMolecule = loadNewMolecule(param, componentPath.cast<std::string>());
		molecules.push_back(std::make_shared<molecule>(newMolecule));
	}

	if (molecules.size() == 0) {
		throw std::runtime_error("Please load at least one molecule.");
	}

}

void loadCalculationsOnPython(py::list calculationsOnPython) {

	n_ex += 1;

	if (n_ex != 2) {
		throw std::runtime_error("loadCalculations should only be executed once after calling loadMolecules.");
	}

	// load the calculations
	// calculationsOnPython list of dictionaries
	// calculationsList[i]["components"] = np.array([0 2])
	// calculationsList[i]["temperatures"] = np.array([298.15])
	// calculationsList[i]["concentrations"] = np.array([0.2 0.8; 0.5 0.5; 0.1 0.9])
	// calculationsList[i]["reference_state_types"] = np.array([0 2 3])  0: PureComponents | 1: PureComponentsOnlyNeutral | 2: ReferenceMixture (also for InfiniteDilution) | 3: COSMO
	// calculationsList[i]["reference_state_concentrations"] = np.array([[] [] [0 1]])  # np.array of concentrations where referenceStateType == 2

	const size_t numCalcs = calculationsOnPython.size();

	if (numCalcs == 0) {
		throw std::runtime_error("Please specify at least one calculation.");
	}

	for (int i = 0; i < numCalcs; i++) {

		py::dict calculationDict = calculationsOnPython[i];

		// array of component indices
 		py::list componentList = calculationDict["component_indices"];
		int numberOfComponents = int(componentList.size());

		calculation newCalculation(numberOfComponents);

		for (int j = 0; j < numberOfComponents; j++) {
			std::shared_ptr<molecule> thisMolecule = molecules[componentList[j].cast<int>()];

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
		auto temperatures = py::array_t<double>(calculationDict["temperatures"]).unchecked<1>();
		auto concentrations = py::array_t<double>(calculationDict["concentrations"]).unchecked<2>();

		for (int j = 0; j < (size_t)concentrations.shape(0); j++) {

			std::vector<float> rowConcentration;

			float tempSumOfConcentrations = 0;
			for (int k = 0; k < numberOfComponents; k++) {
				float val = (float)concentrations(j, k);
				tempSumOfConcentrations += val;
				rowConcentration.push_back(val);
			}

			if (abs(1.0f - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
				throw std::runtime_error("For calculation number " + std::to_string(i) + ", the concentrations do not add up to unity. residual concentration: " + std::to_string(abs(1.0f - tempSumOfConcentrations)));
			}

			float temperature = (float)temperatures(j);
			
			newCalculation.temperatures.push_back(temperature);
			newCalculation.concentrations.push_back(rowConcentration);
		}

		auto referenceStateConcentrations = py::array_t<double>(calculationDict["reference_state_concentrations"]).unchecked<2>();

		if (referenceStateConcentrations.shape(0) != newCalculation.concentrations.size()) {
			throw std::runtime_error("concentrations and referenceStateConcentrations of calculation number " + std::to_string(i) + " have different sizes.\n");
		}

		newCalculation.originalNumberOfCalculations = (unsigned short)newCalculation.concentrations.size();

		// reference states
		auto referenceStateTypes = py::array_t<int>(calculationDict["reference_state_types"]).unchecked<1>();
		for (int j = 0; j < (size_t)referenceStateTypes.shape(0); j++) {


			int referenceStateType = referenceStateTypes(j);

			newCalculation.referenceStateType.push_back((unsigned short)referenceStateType);
			
			float tempSumOfConcentrations = 0;
			if (calculationDict.contains("reference_state_concentrations")) {
				for (int k = 0; k < (size_t)referenceStateConcentrations.shape(1); k++) {
					tempSumOfConcentrations += (float)referenceStateConcentrations(j, k);
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

					float temperature = (float)temperatures(j);
					int referenceStateCalculationIndex = (int)newCalculation.addOrFindArrayIndexForConcentration(referenceStateConcentration, temperature);
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

						float temperature = (float)temperatures(j);
						int referenceStateCalculationIndex = (int)newCalculation.addOrFindArrayIndexForConcentration(referenceStateConcentration, temperature);
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
					referenceStateConcentration.push_back((float)referenceStateConcentrations(j, m));
				}

				if (referenceStateConcentrations.shape(1) != newCalculation.components.size()) {
					throw std::runtime_error("A reference state concentration was specified with the wrong amount of concentrations.");
				}

				if (abs(1.0f - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
					throw std::runtime_error("For calculation number " + std::to_string(i) + ", the reference concentrations do not add up to unity. residual concentration: " + std::to_string(abs(1.0f - tempSumOfConcentrations)));
				}

				float temperature = (float)temperatures(j);
				int referenceStateCalculationIndex = (int)newCalculation.addOrFindArrayIndexForConcentration(referenceStateConcentration, temperature);

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

		// directly bind to python numpy arrays
		// for this to work correctly the sizes of the n-dimensional numpy arrays and the type
		// must be the same in python and c++ (float32/float) with same storage order: row major.

		py::array_t<float> tempArray = py::array_t<float>(calculationDict["ln_gamma_x_SR_combinatorial_calc"]);
		new (&newCalculation.lnGammaCombinatorial) Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(tempArray.mutable_data(),
			int(newCalculation.originalNumberOfCalculations),
			int(newCalculation.components.size()));

		tempArray = py::array_t<float>(calculationDict["ln_gamma_x_SR_residual_calc"]);
		new (&newCalculation.lnGammaResidual) Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(tempArray.mutable_data(),
			int(newCalculation.originalNumberOfCalculations),
			int(newCalculation.components.size()));

		tempArray = py::array_t<float>(calculationDict["ln_gamma_x_SR_calc"]);
		new (&newCalculation.lnGammaTotal) Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(tempArray.mutable_data(),
			int(newCalculation.originalNumberOfCalculations),
			int(newCalculation.components.size()));

		if (calculationDict.contains("dGsolv")) {
			tempArray = py::array_t<float>(calculationDict["dGsolv"]);
			new (&newCalculation.dGsolv) Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(tempArray.mutable_data(),
				int(newCalculation.originalNumberOfCalculations),
				1);
		}

		if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {

			tempArray = py::array_t<float>(calculationDict["contact_statistics"]);
			new (&newCalculation.contactStatistics) Eigen::TensorMap<Eigen::Tensor<float, 3, Eigen::RowMajor>>(tempArray.mutable_data(),
				int(newCalculation.originalNumberOfCalculations),
				int(newCalculation.components.size()),
				int(newCalculation.components.size()));

			tempArray = py::array_t<float>(calculationDict["average_surface_energies"]);
			new (&newCalculation.averageSurfaceEnergies) Eigen::TensorMap<Eigen::Tensor<float, 4, Eigen::RowMajor>>(tempArray.mutable_data(),
				int(newCalculation.originalNumberOfCalculations),
				int(param.numberOfPartialInteractionMatrices) + 1, // +1 because A_int is the first one
				int(newCalculation.components.size()),
				int(newCalculation.components.size()));

			if (param.sw_calculateContactStatisticsAndAdditionalProperties == 2) {

				tempArray = py::array_t<float>(calculationDict["partial_molar_energies"]);
				new (&newCalculation.partialMolarEnergies) Eigen::TensorMap<Eigen::Tensor<float, 3, Eigen::RowMajor>>(tempArray.mutable_data(),
				int(newCalculation.originalNumberOfCalculations),
				int(param.numberOfPartialInteractionMatrices) + 1, // +1 because A_int is the first one
				int(newCalculation.components.size()));
			}
		}
		
		newCalculation.number = (int)i;
		finishCalculationInitiation(newCalculation);
		calculations.push_back(newCalculation);
	}

}

py::list calculateOnPython(py::dict parameters, py::list calculationsOnPython, bool reloadConcentrations = false, bool reloadReferenceConcentrations = false) {

	n_ex += 1;

	if (n_ex < 3) {
		throw std::runtime_error("Before trying to run a calculation please first execute initiate and loadMolecules");
	}
#ifdef MEASURE_TIME
	startCalculationMeasurement();
#endif

	loadParametersOnPython(parameters);

	if (param.sw_alwaysReloadSigmaProfiles == 1 && n_ex > 3) {
		reloadAllMolecules();
	}

	if (param.sw_alwaysCalculateSizeRelatedParameters == 1 || (param.sw_alwaysCalculateSizeRelatedParameters == 0 && n_ex == 3)) {
		resizeMonoatomicCations(param, molecules);
	}

	const size_t numCalcs = calculationsOnPython.size();
	std::vector<int> calculationIndices(numCalcs);

	for (int i = 0; i < numCalcs; i++) {

		calculationIndices[i] = calculationsOnPython[i]["index"].cast<int>();

		param.sw_reloadConcentrations = 0;
		param.sw_reloadReferenceConcentrations = 0;

		if (reloadConcentrations == true) {
			param.sw_reloadConcentrations = 1;
			auto concentrations = py::array_t<double>(calculationsOnPython[i]["concentrations"]).unchecked<2>();
			for (int h = 0; h < calculations[calculationIndices[i]].originalNumberOfCalculations; h++) {

				int j = calculations[calculationIndices[i]].actualConcentrationIndices[h];

				std::vector<float> rowConcentration;

				float tempSumOfConcentrations = 0;
				for (int k = 0; k < calculations[calculationIndices[i]].components.size(); k++) {
					float val = (float)concentrations(h, k);
					tempSumOfConcentrations += val;
					calculations[calculationIndices[i]].concentrations[j][k] = val;
				}

				if (abs(1.0f - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
					throw std::runtime_error("A concentration does not add up to unity.");
				}
			}
		}
		if (reloadReferenceConcentrations == true) {
			param.sw_reloadReferenceConcentrations = 1;
			auto referenceStateTypes = py::array_t<int>(calculationsOnPython[i]["reference_state_types"]).unchecked<1>();
			auto referenceStateConcentrations = py::array_t<double>(calculationsOnPython[i]["reference_state_concentrations"]).unchecked<2>();


			for (int j = 0; j < (size_t)referenceStateConcentrations.shape(0); j++) {
				int referenceStateType = referenceStateTypes(j);
				if (referenceStateType != 2) {
					throw std::runtime_error("reloadReferenceConcentrations only makes sense if the referenceStateTypes == 2.\n");
				}
			}
			if (calculations[calculationIndices[i]].concentrations.size() != calculations[calculationIndices[i]].originalNumberOfCalculations * 2) {
				throw std::runtime_error("The implementation currently assumes that every concentratoin has a unique reference concentration, this could and should be changed in the future.\n");
			}

			for (int h = 0; h < calculations[calculationIndices[i]].originalNumberOfCalculations; h++) {

				std::vector<int> referenceStateCalculationIndices = calculations[calculationIndices[i]].referenceStateCalculationIndices[h];

				float tempSumOfConcentrations = 0;
				for (int k = 0; k < calculations[calculationIndices[i]].components.size(); k++) {
					int referenceStateCalculationIndex = referenceStateCalculationIndices[k];
					float val = (float)referenceStateConcentrations(h, k);
					tempSumOfConcentrations += val;
					calculations[calculationIndices[i]].concentrations[referenceStateCalculationIndex][k] = val;
				}

				if (abs(1.0f - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
					throw std::runtime_error("A concentration does not add up to unity.");
				}
			}
		}
	}

	calculate(calculationIndices);

#ifdef MEASURE_TIME
	stopCalculationMeasurement();
#endif

	return calculationsOnPython;

}

PYBIND11_MODULE(openCOSMORS, m) {
	m.doc() = R"pbdoc(
        openCOSMO-RS
        -----------------------

        .. currentmodule:: openCOSMORS

        .. autosummary::
           :toctree: _generate

           initialize
    )pbdoc";

	m.def("initialize", &initializeOnPython, R"pbdoc(
        Sets the stage to start running the module again.
    )pbdoc");

	m.def("loadMolecules", &loadMoleculesOnPython, R"pbdoc(
        Loads all the sigma profiles of the molecules.
		This needs to be called before calling loadCalculations.
    )pbdoc");

	m.def("loadCalculations", &loadCalculationsOnPython, R"pbdoc(
        Loads all calculations.
		This needs to be called before calling calculate.
    )pbdoc");
	m.def("calculate", &calculateOnPython, py::arg("parameters"), py::arg("calculationsOnPython"), py::arg("reloadConcentrations") = false, py::arg("reloadReferenceConcentrations") = false, py::return_value_policy::reference, R"pbdoc(
        Calculates the complete list of calculations with the provided set of parameters.
		These calculations should have been loaded with loadCalculations prior to executing calculate, otherwise this will produce an error.
    )pbdoc");


#ifdef VERSION_INFO
	m.attr("__version__") = VERSION_INFO;
#else
	m.attr("__version__") = "dev";
#endif
}