/*
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @author: Simon Mueller, 2022
*/


//#define MEASURE_TIME
//#define DEBUG_INFO

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
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
	param.alpha = parameters["alpha"].cast<double>();

	param.CHB = parameters["CHB"].cast<double>();
	param.CHBT = parameters["CHBT"].cast<double>();
	param.SigmaHB = parameters["SigmaHB"].cast<double>();

	param.Rav = parameters["Rav"].cast<double>();

	if (param.sw_SR_polarizabilities > 0) {
		param.E_F_corr = parameters["E_F_corr"].cast<double>();
		param.m_vdW = parameters["m_vdW"].cast<double>();
	}

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
			std::string key = py::str(item.first).cast<std::string>();
			param.dGsolv_tau[parseAtomicNumberKey(key, "dGsolv_tau")] = item.second.cast<double>();
		}

		// one entry per loaded molecule (order of the component paths)
		if (!parameters.contains("dGsolv_E_gas") || !parameters.contains("dGsolv_numberOfAtomsInRing")) {
			throw std::runtime_error("To calculate solvation energies the parameters 'dGsolv_E_gas' and 'dGsolv_numberOfAtomsInRing' (one entry per loaded molecule) are needed in addition to 'dGsolv_eta'.");
		}
		param.dGsolv_E_gas = parameters["dGsolv_E_gas"].cast<std::vector<double>>();
		param.dGsolv_numberOfAtomsInRing = parameters["dGsolv_numberOfAtomsInRing"].cast<std::vector<int>>();
	}

	if (parameters.contains("radii")) {
		py::dict radii = parameters["radii"];

		for (auto item : radii) {
			std::string key = py::str(item.first).cast<std::string>();
			param.R_i[parseAtomicNumberKey(key, "radii")] = item.second.cast<double>();
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

	std::vector<std::string> partialInteractionMatrixNames;
	if (param.sw_calculateContactStatisticsAndAdditionalProperties != 0) {
		partialInteractionMatrixNames = options["sw_SR_partialInteractionMatrices"].cast<std::vector<std::string>>();
	}
	setPartialInteractionMatrices(param, partialInteractionMatrixNames);

	param.sw_atomicNumber = options["sw_SR_atomicNumber"].cast<int>();
	param.sw_misfit = options["sw_SR_misfit"].cast<int>();

	if (param.sw_misfit < 0 || param.sw_misfit > 2) {
		throw std::runtime_error("sw_SR_misfit should have one of the following values: [0, 1, 2].");
	}

	if (param.sw_SR_polarizabilities > 0 && param.sw_atomicNumber == 0) {
		throw std::runtime_error("sw_SR_polarizabilities > 0 requires sw_SR_atomicNumber = 1, as the fluorine correction and the ionization potentials depend on the atomic number.");
	}
	param.sw_skip_COSMOSPACE_errors = options["sw_skip_COSMOSPACE_errors"].cast<int>();

	// parameters
	loadParametersOnPython(parameters);
	validateSolvationEnergyParameters(param, int(componentPaths.size()));

	for (auto componentPath : componentPaths) {
		molecule newMolecule = loadNewMolecule(param, componentPath.cast<std::string>());
		molecules.push_back(std::make_shared<molecule>(newMolecule));
	}

	if (molecules.size() == 0) {
		throw std::runtime_error("Please load at least one molecule.");
	}

}

std::string shapeToString(const std::vector<long>& shape) {
	std::string s = "(";
	for (size_t k = 0; k < shape.size(); k++) s += (k ? ", " : "") + std::to_string(shape[k]);
	return s + ")";
}

// returns the data pointer of an output array the results are written into directly,
// throws if the array cannot be used without a copy
double* getOutputArrayData(py::dict calculationDict, const std::string& key, const std::vector<long>& expectedShape, int calculationNumber) {

	const std::string prefix = "Calculation number " + std::to_string(calculationNumber) + ": the output array '" + key + "' ";
	const std::string expected = " Expected a writeable numpy array with dtype float64, C order and shape " + shapeToString(expectedShape) + ".";

	if (!calculationDict.contains(key.c_str())) {
		throw std::runtime_error(prefix + "is missing." + expected);
	}
	py::object obj = calculationDict[key.c_str()];
	if (!py::isinstance<py::array>(obj)) {
		throw std::runtime_error(prefix + "is not a numpy array." + expected);
	}
	py::array arr = py::reinterpret_borrow<py::array>(obj);

	std::vector<long> shape;
	for (py::ssize_t k = 0; k < arr.ndim(); k++) shape.push_back(long(arr.shape(k)));

	if (!arr.dtype().is(py::dtype::of<double>())) {
		throw std::runtime_error(prefix + "has dtype " + py::str(arr.dtype()).cast<std::string>() + "." + expected);
	}
	if (shape != expectedShape) {
		throw std::runtime_error(prefix + "has shape " + shapeToString(shape) + "." + expected);
	}
	if (!(arr.flags() & py::array::c_style)) {
		throw std::runtime_error(prefix + "is not C-contiguous (e.g. Fortran order, transposed or sliced)." + expected);
	}
	if (!arr.writeable()) {
		throw std::runtime_error(prefix + "is read-only." + expected);
	}

	return static_cast<double*>(arr.mutable_data());
}

// The results are written directly into the numpy arrays of the calculation dictionary, therefore no copies
// are allowed: the arrays have to exist with dtype float64, C order and exactly the expected shape.
// This is done in loadCalculations and again in every calculate, so that the results are always written into
// the arrays currently stored in the dictionary (the dictionaries passed to calculate keep them alive).
void bindOutputArrays(calculation& _calculation, py::dict calculationDict, int calculationNumber) {

	const long nRows = long(_calculation.originalNumberOfCalculations);
	const long nComp = long(_calculation.components.size());
	const long nMatrices = long(param.numberOfPartialInteractionMatrices) + 1; // +1 because A_int is the first one

	using MapType = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;

	double* data = getOutputArrayData(calculationDict, "ln_gamma_x_SR_combinatorial_calc", { nRows, nComp }, calculationNumber);
	new (&_calculation.lnGammaCombinatorial) MapType(data, nRows, nComp);

	data = getOutputArrayData(calculationDict, "ln_gamma_x_SR_residual_calc", { nRows, nComp }, calculationNumber);
	new (&_calculation.lnGammaResidual) MapType(data, nRows, nComp);

	data = getOutputArrayData(calculationDict, "ln_gamma_x_SR_calc", { nRows, nComp }, calculationNumber);
	new (&_calculation.lnGammaTotal) MapType(data, nRows, nComp);

	if (calculationDict.contains("dGsolv")) {
		data = getOutputArrayData(calculationDict, "dGsolv", { nRows, nComp }, calculationNumber);
		new (&_calculation.dGsolv) MapType(data, nRows, nComp);
	}
	else {
		new (&_calculation.dGsolv) MapType(nullptr, 0, 0);
	}

	if (param.sw_calculateContactStatisticsAndAdditionalProperties > 0) {

		data = getOutputArrayData(calculationDict, "contact_statistics", { nRows, nComp, nComp }, calculationNumber);
		new (&_calculation.contactStatistics) Eigen::TensorMap<Eigen::Tensor<double, 3, Eigen::RowMajor>>(data, nRows, nComp, nComp);

		data = getOutputArrayData(calculationDict, "average_surface_energies", { nRows, nMatrices, nComp, nComp }, calculationNumber);
		new (&_calculation.averageSurfaceEnergies) Eigen::TensorMap<Eigen::Tensor<double, 4, Eigen::RowMajor>>(data, nRows, nMatrices, nComp, nComp);

		if (param.sw_calculateContactStatisticsAndAdditionalProperties == 2) {

			data = getOutputArrayData(calculationDict, "partial_molar_energies", { nRows, nMatrices, nComp }, calculationNumber);
			new (&_calculation.partialMolarEnergies) Eigen::TensorMap<Eigen::Tensor<double, 3, Eigen::RowMajor>>(data, nRows, nMatrices, nComp);
		}
	}
}

void loadCalculationsOnPython(py::list calculationsOnPython, bool reload = false) {

	n_ex += 1;

	if (reload) {
		initialize(param, false, false, true, false);
		// state after loadMolecules and loadCalculations (n_ex only checks the order of the calls,
		// the new calculations are initialized in the next calculation as they are new)
		n_ex = 2;
	}
	else {
		if (n_ex != 2) {
			throw std::runtime_error("loadCalculations should only be executed once after calling loadMolecules.");
		}
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

	// important as otherwise the data behind Eigen::Maps
	// can be lost when the vector is resized.
	calculations.reserve(numCalcs);

	for (int i = 0; i < numCalcs; i++) {

		py::dict calculationDict = calculationsOnPython[i];

		// array of component indices
		std::vector<int> componentIndices = calculationDict["component_indices"].cast<std::vector<int>>();
		int numberOfComponents = int(componentIndices.size());

		// input arrays, the converted arrays are kept alive as long as their unchecked views are used
		py::array_t<double, py::array::c_style | py::array::forcecast> temperaturesArray(calculationDict["temperatures"]);
		py::array_t<double, py::array::c_style | py::array::forcecast> concentrationsArray(calculationDict["concentrations"]);
		py::array_t<int, py::array::c_style | py::array::forcecast> referenceStateTypesArray(calculationDict["reference_state_types"]);
		py::array_t<double, py::array::c_style | py::array::forcecast> referenceStateConcentrationsArray(calculationDict["reference_state_concentrations"]);

		if (temperaturesArray.ndim() != 1 || concentrationsArray.ndim() != 2 || referenceStateTypesArray.ndim() != 1 || referenceStateConcentrationsArray.ndim() != 2) {
			throw std::runtime_error("Calculation number " + std::to_string(i) + ": 'temperatures' and 'reference_state_types' have to be one-dimensional, 'concentrations' and 'reference_state_concentrations' two-dimensional.");
		}

		validateCalculationSizes(i, numberOfComponents, componentIndices, int(molecules.size()),
			int(concentrationsArray.shape(0)), int(concentrationsArray.shape(1)), int(temperaturesArray.shape(0)),
			int(referenceStateTypesArray.shape(0)), int(referenceStateConcentrationsArray.shape(0)));

		calculation newCalculation(numberOfComponents);
		newCalculation.componentIndices = componentIndices;

		for (int j = 0; j < numberOfComponents; j++) {
			std::shared_ptr<molecule> thisMolecule = molecules[componentIndices[j]];

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
		auto temperatures = temperaturesArray.unchecked<1>();
		auto concentrations = concentrationsArray.unchecked<2>();

		for (int j = 0; j < (size_t)concentrations.shape(0); j++) {

			std::vector<double> rowConcentration;

			double tempSumOfConcentrations = 0;
			for (int k = 0; k < numberOfComponents; k++) {
				double val = (double)concentrations(j, k);
				tempSumOfConcentrations += val;
				rowConcentration.push_back(val);
			}

			if (std::abs(1.0 - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
				throw std::runtime_error("For calculation number " + std::to_string(i) + ", the concentrations do not add up to unity. residual concentration: " + std::to_string(std::abs(1.0 - tempSumOfConcentrations)));
			}

			double temperature = (double)temperatures(j);
			
			newCalculation.temperatures.push_back(temperature);
			newCalculation.concentrations.push_back(rowConcentration);
		}

		auto referenceStateConcentrations = referenceStateConcentrationsArray.unchecked<2>();

		newCalculation.originalNumberOfCalculations = (unsigned short)newCalculation.concentrations.size();
		const int startIndexOfReferenceStates = int(newCalculation.originalNumberOfCalculations);

		// reference states
		auto referenceStateTypes = referenceStateTypesArray.unchecked<1>();
		for (int j = 0; j < (size_t)referenceStateTypes.shape(0); j++) {


			int referenceStateType = referenceStateTypes(j);

			newCalculation.referenceStateType.push_back((unsigned short)referenceStateType);
			
			double tempSumOfConcentrations = 0;
			if (calculationDict.contains("reference_state_concentrations")) {
				for (int k = 0; k < (size_t)referenceStateConcentrations.shape(1); k++) {
					tempSumOfConcentrations += (double)referenceStateConcentrations(j, k);
				}
			}


			if (referenceStateType == 0) { // Pure component
				
				if (tempSumOfConcentrations != 0) {
					throw std::runtime_error("A reference state concentration was specified for a calculation with reference state PureComponents, this does not make sense.");
				}

				std::vector<int> thisReferenceStateCalculationIndices;
				for (int k = 0; k < numberOfComponents; k++) {

					std::vector<double> referenceStateConcentration;
					for (int m = 0; m < numberOfComponents; m++) {
						referenceStateConcentration.push_back(k == m ? 1.0f : 0.0f);
					}

					double temperature = (double)temperatures(j);
					int referenceStateCalculationIndex = (int)newCalculation.addOrFindArrayIndexForConcentration(referenceStateConcentration, temperature, startIndexOfReferenceStates);
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

						std::vector<double> referenceStateConcentration;
						for (int m = 0; m < numberOfComponents; m++) {
							referenceStateConcentration.push_back(k == m ? 1.0f : 0.0f);
						}

						double temperature = (double)temperatures(j);
						int referenceStateCalculationIndex = (int)newCalculation.addOrFindArrayIndexForConcentration(referenceStateConcentration, temperature, startIndexOfReferenceStates);
						thisReferenceStateCalculationIndices.push_back(referenceStateCalculationIndex);
					}
					else {
						thisReferenceStateCalculationIndices.push_back(-1);
					}

				}
				newCalculation.referenceStateCalculationIndices.push_back(thisReferenceStateCalculationIndices);
			}
			else if (referenceStateType == 2) { // Reference mixture

				if (referenceStateConcentrations.shape(1) != newCalculation.components.size()) {
					throw std::runtime_error("Calculation number " + std::to_string(i) + ", concentration " + std::to_string(j) + ": reference state type 2 needs a reference state concentration with " +
						std::to_string(newCalculation.components.size()) + " entries, but 'reference_state_concentrations' has " + std::to_string(referenceStateConcentrations.shape(1)) + " columns.");
				}

				std::vector<double> referenceStateConcentration;
				for (int m = 0; m < numberOfComponents; m++) {
					referenceStateConcentration.push_back((double)referenceStateConcentrations(j, m));
				}

				if (std::abs(1.0 - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
					throw std::runtime_error("For calculation number " + std::to_string(i) + ", the reference concentrations do not add up to unity. residual concentration: " + std::to_string(std::abs(1.0 - tempSumOfConcentrations)));
				}

				double temperature = (double)temperatures(j);
				int referenceStateCalculationIndex = (int)newCalculation.addOrFindArrayIndexForConcentration(referenceStateConcentration, temperature, startIndexOfReferenceStates);

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

		bindOutputArrays(newCalculation, calculationDict, i);

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

	// the solvation energy parameters can be changed here, the molar volume of the solvent is however
	// only calculated when loading the molecules, hence they have to be given there already
	validateSolvationEnergyParameters(param, int(molecules.size()));
	if (param.dGsolv_E_gas.size() > 0) {
		for (const std::shared_ptr<molecule>& _molecule : molecules) {
			if (std::isnan(_molecule->molarVolumeAt25C)) {
				throw std::runtime_error("To calculate solvation energies the parameters dGsolv_eta, dGsolv_E_gas, ... have to be given already when calling loadMolecules.");
			}
		}
	}

	prepareMoleculesForCalculation(param);

	const size_t numCalcs = calculationsOnPython.size();
	std::vector<int> calculationIndices(numCalcs);

	for (int i = 0; i < numCalcs; i++) {

		calculationIndices[i] = calculationsOnPython[i]["index"].cast<int>();

		if (calculationIndices[i] < 0 || calculationIndices[i] >= int(calculations.size())) {
			throw std::runtime_error("The calculation at position " + std::to_string(i) + " has the index " + std::to_string(calculationIndices[i]) +
				", but only " + std::to_string(calculations.size()) + " calculations were loaded.");
		}

		calculation& thisCalculation = calculations[calculationIndices[i]];
		const int nRows = int(thisCalculation.originalNumberOfCalculations);
		const int nComp = int(thisCalculation.components.size());

		// write the results into the arrays currently stored in the dictionary
		bindOutputArrays(thisCalculation, calculationsOnPython[i].cast<py::dict>(), calculationIndices[i]);

		param.sw_reloadConcentrations = 0;
		param.sw_reloadReferenceConcentrations = 0;

		if (reloadConcentrations == true) {
			param.sw_reloadConcentrations = 1;
			// kept alive as long as its unchecked view is used
			py::array_t<double, py::array::c_style | py::array::forcecast> concentrationsArray(calculationsOnPython[i]["concentrations"]);
			if (concentrationsArray.ndim() != 2 || concentrationsArray.shape(0) != nRows || concentrationsArray.shape(1) != nComp) {
				throw std::runtime_error("Calculation number " + std::to_string(calculationIndices[i]) + ": reloaded 'concentrations' must keep the shape (" +
					std::to_string(nRows) + ", " + std::to_string(nComp) + ").");
			}
			auto concentrations = concentrationsArray.unchecked<2>();
			for (int h = 0; h < nRows; h++) {

				// concentrations are stored sorted internally
				int j = thisCalculation.actualConcentrationIndices[h];

				double tempSumOfConcentrations = 0;
				for (int k = 0; k < nComp; k++) {
					double val = (double)concentrations(h, k);
					tempSumOfConcentrations += val;
					thisCalculation.concentrations[j][k] = val;
				}

				if (std::abs(1.0 - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
					throw std::runtime_error("Calculation number " + std::to_string(calculationIndices[i]) + ", concentration " + std::to_string(h) +
						": the reloaded concentrations do not add up to unity. residual concentration: " + std::to_string(std::abs(1.0 - tempSumOfConcentrations)));
				}
			}
		}
		if (reloadReferenceConcentrations == true) {
			param.sw_reloadReferenceConcentrations = 1;
			py::array_t<int, py::array::c_style | py::array::forcecast> referenceStateTypesArray(calculationsOnPython[i]["reference_state_types"]);
			py::array_t<double, py::array::c_style | py::array::forcecast> referenceStateConcentrationsArray(calculationsOnPython[i]["reference_state_concentrations"]);

			if (referenceStateTypesArray.ndim() != 1 || referenceStateTypesArray.shape(0) != nRows ||
				referenceStateConcentrationsArray.ndim() != 2 || referenceStateConcentrationsArray.shape(0) != nRows || referenceStateConcentrationsArray.shape(1) != nComp) {
				throw std::runtime_error("Calculation number " + std::to_string(calculationIndices[i]) + ": reloaded 'reference_state_types' must have " + std::to_string(nRows) +
					" entries and 'reference_state_concentrations' the shape (" + std::to_string(nRows) + ", " + std::to_string(nComp) + ").");
			}
			auto referenceStateTypes = referenceStateTypesArray.unchecked<1>();
			auto referenceStateConcentrations = referenceStateConcentrationsArray.unchecked<2>();

			for (int j = 0; j < nRows; j++) {
				if (referenceStateTypes(j) != 2 || thisCalculation.referenceStateType[j] != 2) {
					throw std::runtime_error("Calculation number " + std::to_string(calculationIndices[i]) + ", concentration " + std::to_string(j) +
						": reloadReferenceConcentrations is only possible if all reference state types are 2.");
				}
			}
			if (thisCalculation.concentrations.size() != size_t(nRows) * 2) {
				throw std::runtime_error("Calculation number " + std::to_string(calculationIndices[i]) +
					": reloadReferenceConcentrations requires every concentration to have its own unique reference concentration.");
			}

			for (int h = 0; h < nRows; h++) {

				const std::vector<int>& referenceStateCalculationIndices = thisCalculation.referenceStateCalculationIndices[h];

				double tempSumOfConcentrations = 0;
				for (int k = 0; k < nComp; k++) {
					// the stored index refers to the unsorted concentrations, they are stored sorted internally
					int referenceStateCalculationIndex = thisCalculation.actualConcentrationIndices[referenceStateCalculationIndices[k]];
					double val = (double)referenceStateConcentrations(h, k);
					tempSumOfConcentrations += val;
					thisCalculation.concentrations[referenceStateCalculationIndex][k] = val;
				}

				if (std::abs(1.0 - tempSumOfConcentrations) > MAX_CONCENTRATION_DIFF_FROM_ZERO) {
					throw std::runtime_error("Calculation number " + std::to_string(calculationIndices[i]) + ", concentration " + std::to_string(h) +
						": the reloaded reference concentrations do not add up to unity. residual concentration: " + std::to_string(std::abs(1.0 - tempSumOfConcentrations)));
				}
			}
		}
	}

	if (reloadConcentrations || reloadReferenceConcentrations) {
		for (int calculationIndex : calculationIndices) {
			checkConcentrationsAndTemperatures(calculations[calculationIndex]);
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

	m.def("loadCalculations", &loadCalculationsOnPython, py::arg("calculationsOnPython"), py::arg("reload") = false, R"pbdoc(
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