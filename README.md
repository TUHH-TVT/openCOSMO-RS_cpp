# openCOSMO-RS C++ with Python bindings

This repository contains a C++ implementation of the COSMO-RS model with a
pybind11 Python interface. It is much faster than the
[pure Python implementation](https://github.com/TUHH-TVT/openCOSMO-RS_py), and
newer model features and parameterizations are added here first.

The standalone binary is useful for simple JSON-driven batch calculations and
is currently being used by ORCA. MATLAB users should call the Python
module from MATLAB.

## Contents

- [Quick start: Python module](#quick-start-python-module)
- [Build options](#build-options)
- [Standalone binary and JSON input](#standalone-binary-and-json-input)
- [Python API reference](#python-api-reference)
- [MATLAB usage through Python](#matlab-usage-through-python)
- [Related projects](#related-projects)

## Quick start: Python module

### 1. Clone with submodules

```bash
git clone --recursive --shallow-submodules https://github.com/TUHH-TVT/openCOSMO-RS_cpp.git
cd openCOSMO-RS_cpp
```

If the repository was cloned without submodules:

```bash
git submodule update --init --recursive
```

### 2. Install basic Python dependencies

```bash
python -m pip install --upgrade pip numpy
```

### 3. Build the Python extension

The default CMake target builds the `openCOSMORS` Python extension.

```bash
cmake -S . -B build -DOPENCOMSORS_SIMD=AUTO
cmake --build build --config Release
```

### 4. Run the bundled example

Add the build output directory to `PYTHONPATH`, then run
[`bindings/run_example.py`](bindings/run_example.py).

Linux/macOS single-config generators usually place the extension directly in
`build`:

```bash
PYTHONPATH="$PWD/build:$PYTHONPATH" python bindings/run_example.py
```

Windows and other multi-config generators usually place the extension in
`build/Release`:

```powershell
$env:PYTHONPATH = "$PWD\build\Release;$env:PYTHONPATH"
python bindings\run_example.py
```

## Build options

### SIMD backend

The SIMD backend is selected with `OPENCOMSORS_SIMD`:

| Value | Meaning |
| --- | --- |
| `AUTO` | Default. Chooses `FMA` on x86/x86_64, `NEON` on ARM64, and `SCALAR` on unknown architectures. |
| `SCALAR` | Portable fallback for debugging, reference checks, or old CPUs. |
| `SSE3` | x86 SSE3 backend. |
| `AVX` | x86 AVX backend. |
| `FMA` | x86 AVX/FMA backend. |
| `NEON` | AArch64/ARM64 NEON backend, including Apple Silicon. |

Examples:

```bash
cmake -S . -B build -DOPENCOMSORS_SIMD=SCALAR
cmake -S . -B build -DOPENCOMSORS_SIMD=FMA
cmake -S . -B build -DOPENCOMSORS_SIMD=NEON
```

### Platform notes

**macOS**

Install command-line tools first:

```bash
xcode-select --install
```

On Apple Silicon, `AUTO` selects `NEON`. On Intel Macs, `AUTO` selects `FMA`.

**Windows**

Use Visual Studio with CMake. Activate your conda or virtual environment before
running CMake if you want pybind11 to bind against that Python installation.

```powershell
cmake -S . -B build -DOPENCOMSORS_SIMD=AUTO
cmake --build build --config Release
```

**Linux**

Install a compiler, CMake, Python development headers, and NumPy. If the target
machine does not support FMA, use `SCALAR`, `SSE3`, or `AVX` explicitly.

```bash
cmake -S . -B build -DOPENCOMSORS_SIMD=AUTO
cmake --build build --config Release
```

OpenMP is used for release builds when CMake can find it. If OpenMP is not
available, the project still builds and runs without OpenMP parallelization.

## Standalone binary and JSON input

### Build the binary

Pass `-DBINARY=1` to build the command-line executable instead of the Python
module:

```bash
cmake -S . -B build-cli -DBINARY=1 -DOPENCOMSORS_SIMD=AUTO
cmake --build build-cli --config Release
```

### Run the binary

The executable accepts an input JSON file and an optional output JSON file:

```bash
./build-cli/openCOSMORS input.json output.json
```

Windows multi-config builds usually place the executable under `Release`:

```powershell
.\build-cli\Release\openCOSMORS.exe input.json output.json
```

If the output path is omitted, the binary writes next to the input file using
the suffix `_out.json`, for example `input_out.json`.

### CLI output format

The CLI JSON format is not the same as the Python in-memory API. The binary
allocates result arrays internally and currently writes a compact output JSON:

```json
{
  "dGsolv": [
    [
      [0.0, 0.0, 0.0]
    ]
  ],
  "warnings": []
}
```

`dGsolv` shape:

- outer list: one entry per calculation.
- middle list: one entry per concentration/temperature row.
- inner list: one value per selected component.

If solvation-energy parameters such as `dGsolv_E_gas`, `dGsolv_tau`,
`dGsolv_eta`, and `dGsolv_numberOfAtomsInRing` are not provided, `dGsolv`
values remain zero. Use the Python API for activity-coefficient outputs and
detailed tensor outputs.

### Minimal CLI input JSON

This example uses the bundled COSMO files:

```json
{
  "sw_COSMOfiles_type": "ORCA_COSMO_TZVPD",
  "sw_SR_partialInteractionMatrices": [],

  "Aeff": 4.90825,
  "ln_alpha": 15.879012208938363,
  "ln_CHB": 17.71334077178702,
  "CHBT": 1.5,
  "SigmaHB": 0.009953,
  "Rav": 0.5,
  "RavCorr": 1.0,
  "fCorr": 2.4,
  "comb_SG_z_coord": 0.0,
  "comb_SG_A_std": 1.0,
  "m_vdW": 29.567,
  "E_F_corr": 346.82,
  "radii": {},
  "exp": {},

  "componentPaths": [
    "bindings/1112tetrachloroethane.orcacosmo",
    "bindings/methanol.orcacosmo",
    "bindings/water.orcacosmo"
  ],

  "calculations": [
    {
      "component_indices": [0, 1, 2],
      "concentrations": [[0.2, 0.5, 0.3]],
      "temperatures": [298.15],
      "reference_state_types": [3],
      "reference_state_concentrations": [[]]
    }
  ]
}
```

CLI input notes:

- `componentPaths` are relative to the process working directory unless
  absolute paths are used.
- `component_indices` select entries from `componentPaths`.
- `concentrations` has shape `(n_points, n_components)` and each row must sum
  to one.
- `temperatures` has shape `(n_points,)`.
- `reference_state_types` has shape `(n_points,)`.
- `reference_state_concentrations` must have one row per concentration row.
  For reference-state types `0`, `1`, `3`, and `4`, use an empty row such as
  `[[]]` for one point. For type `2`, provide explicit compositions with shape
  `(n_points, n_components)`, for example `[[0.0, 1.0, 0.0]]`.
- The CLI uses compiled C++ defaults for switches not explicitly handled in
  `code/bindings_forCLI.cpp`. The Python API exposes more options and is the
  recommended interface for advanced workflows.

## Python API reference

### Call sequence

The pybind11 API is intentionally small and dictionary-based:

```python
import openCOSMORS

openCOSMORS.loadMolecules(options, parameters, components)
openCOSMORS.loadCalculations(calculations)
calculations = openCOSMORS.calculate(parameters, calculations)
```

The design keeps the C++ model as the single source of truth while letting
Python users work with dictionaries, lists, and NumPy arrays. This avoids
duplicating the data model as Python classes, keeps pybind11 conversion code
thin, and makes it easy to assemble calculations in scripts, notebooks, or
optimization loops.

Lifecycle:

1. `initialize()` resets module state. It is called automatically before the
   first molecule load, but can be called manually when starting over.
2. `loadMolecules(options, parameters, components)` reads COSMO files and
   initializes molecule-level data. Call this once for a component set.
3. `loadCalculations(calculations, reload=False)` registers calculation
   dictionaries and prepares C++/Eigen data structures.
4. `calculate(parameters, calculations, reloadConcentrations=False,
   reloadReferenceConcentrations=False)` runs calculations and writes results
   back into the same dictionaries.

See [`bindings/run_example.py`](bindings/run_example.py) for a complete working
example.

### Calculation dictionary

Each entry in `calculations` is one calculation dictionary. The bundled example
loads three components:

```python
components = [
    "1112tetrachloroethane.orcacosmo",
    "methanol.orcacosmo",
    "water.orcacosmo",
]
```

Minimal one-point calculation:

```python
calculation = {
    "concentrations": np.array([[0.2, 0.5, 0.3]], dtype=np.float64),
    "temperatures": np.array([298.15], dtype=np.float64),
    "components": components,
    "component_indices": [0, 1, 2],
    "reference_state_types": np.array([3]),
    "reference_state_concentrations": np.array([[]], dtype=np.float64),
}
```

Dimensions:

- `n_points`: number of concentration/temperature rows.
- `n_components`: number of selected components.
- `n_interaction_matrices`: `1 + len(options["sw_SR_partialInteractionMatrices"])`.
  The leading matrix is `A_int`.

Required input fields:

| Field | Shape | Notes |
| --- | --- | --- |
| `component_indices` | `(n_components,)` | Integer indices into molecules loaded by `loadMolecules`. |
| `components` | usually `(n_components,)` | Used for Python-side bookkeeping; molecule selection uses `component_indices`. |
| `concentrations` | `(n_points, n_components)` | `float64`; each row must sum to one. |
| `temperatures` | `(n_points,)` | `float64`; row `i` belongs to `concentrations[i, :]`. |
| `reference_state_types` | `(n_points,)` | Integer reference-state code per row. |
| `reference_state_concentrations` | `(n_points, n_components)` or empty second dimension | Explicit only for reference mixture; first dimension must equal `n_points`. |

### Reference-state types

| Type | Meaning | `reference_state_concentrations` row |
| ---: | --- | --- |
| `0` | Pure components. References are generated internally as pure-component compositions. | Must be empty/zero-sum, for example `np.array([[]])` for one point. |
| `1` | Pure components only for neutral molecules. Charged components get no reference calculation. This path is marked untested in the C++ binding. | Must be empty/zero-sum. |
| `2` | Reference mixture. The same explicit reference composition is used for all components in that row. Also use this for reference-mixture/infinite-dilution setups. | Must have shape `(n_components,)` for that row and sum to one, for example `[0.0, 1.0, 0.0]`. |
| `3` | COSMO reference state. No composition reference is generated. | Must be empty/zero-sum. |
| `4` | COSMO reference state for solvation-energy calculations. | Must be empty/zero-sum. |

Examples:

```python
# COSMO reference state.
calculation["reference_state_types"] = np.array([3])
calculation["reference_state_concentrations"] = np.array([[]])

# Pure-component reference states.
calculation["reference_state_types"] = np.array([0])
calculation["reference_state_concentrations"] = np.array([[]])

# Explicit reference mixture: pure methanol in a three-component system.
calculation["reference_state_types"] = np.array([2])
calculation["reference_state_concentrations"] = np.array([[0.0, 1.0, 0.0]])
```

### Result arrays

The Python API expects result arrays to be allocated before
`loadCalculations()`. The binding maps these arrays directly into C++/Eigen data
structures, so shapes must match exactly:

| Field | Required when | Shape |
| --- | --- | --- |
| `index` | always | scalar integer |
| `ln_gamma_x_SR_combinatorial_calc` | always | `(n_points, n_components)` |
| `ln_gamma_x_SR_residual_calc` | always | `(n_points, n_components)` |
| `ln_gamma_x_SR_calc` | always | `(n_points, n_components)` |
| `dGsolv` | optional solvation output | `(n_points, 1)` |
| `contact_statistics` | `sw_SR_calculateContactStatisticsAndAdditionalProperties > 0` | `(n_points, n_components, n_components)` |
| `average_surface_energies` | `sw_SR_calculateContactStatisticsAndAdditionalProperties > 0` | `(n_points, n_interaction_matrices, n_components, n_components)` |
| `partial_molar_energies` | `sw_SR_calculateContactStatisticsAndAdditionalProperties == 2` | `(n_points, n_interaction_matrices, n_components)` |

The helper `fill_missing_calculation_structures()` in
[`bindings/run_example.py`](bindings/run_example.py) shows the recommended
allocation pattern.

### Direct NumPy binding caveats

`loadCalculations()` does not copy result arrays into independent C++ storage.
pybind11 converts Python objects to `py::array_t<double>`, and the C++ code
creates Eigen `Map`/`TensorMap` views over the arrays' memory.

Mapped arrays:

- `ln_gamma_x_SR_combinatorial_calc`
- `ln_gamma_x_SR_residual_calc`
- `ln_gamma_x_SR_calc`
- `dGsolv`, when present
- `contact_statistics`, when enabled
- `average_surface_energies`, when enabled
- `partial_molar_energies`, when enabled

Caveats:

- Allocate mapped arrays with `dtype=np.float64`.
- Use C-contiguous NumPy arrays. Avoid sliced, transposed, Fortran-order,
  object-dtype, or non-contiguous arrays.
- Do not replace, resize, reshape, or reassign mapped result arrays after
  calling `loadCalculations()`.
- Keep the `calculations` list and dictionaries alive until `calculate()`
  finishes.
- To change the number of points or components, allocate a new calculation
  dictionary and call `loadCalculations(..., reload=True)`.
- `calculate()` writes results in-place into the same arrays and returns the
  same Python list for convenience.
- Use `reloadConcentrations=True` or `reloadReferenceConcentrations=True` only
  after updating same-shape input arrays in-place.

## MATLAB usage through Python

MATLAB users should call the Python module instead, matching the recommendation in
[issue #2](https://github.com/TUHH-TVT/openCOSMO-RS_cpp/issues/2).
This keeps the maintained public interface in one place: C++ core plus pybind11.

MathWorks documents two relevant MATLAB-to-Python entry points:

- Use the `py.` prefix to import and call Python modules directly.
- Use `pyrun` or `pyrunfile` to execute Python statements or scripts from
  MATLAB. `pyrunfile` is available from R2021b.

See MathWorks' documentation on
[calling Python from MATLAB](https://www.mathworks.com/help/matlab/call-python-libraries.html),
[directly calling Python functionality](https://www.mathworks.com/help/matlab/matlab_external/ways-to-call-python-from-matlab.html),
and [passing data between MATLAB and Python](https://www.mathworks.com/help/matlab/matlab_external/passing-data-to-python.html).

### MATLAB smoke test

```matlab
repo = "C:\path\to\openCOSMO-RS_cpp";
pythonExe = "C:\path\to\python.exe";

pyenv("Version", pythonExe, "ExecutionMode", "OutOfProcess");

buildDir = fullfile(repo, "build", "Release");  % Windows multi-config build
bindingsDir = fullfile(repo, "bindings");

py.sys.path.insert(int32(0), buildDir);
py.sys.path.insert(int32(0), bindingsDir);

pyrunfile(fullfile(bindingsDir, "run_example.py"));
```

On macOS or Linux, use `build` instead of `build/Release` if the extension is
placed directly in the build directory:

```matlab
buildDir = fullfile(repo, "build");
```

### Recommended MATLAB integration pattern

For production MATLAB usage, write a small Python wrapper module instead of
constructing nested calculation dictionaries directly in MATLAB. The wrapper
can follow `bindings/run_example.py`, allocate NumPy arrays with the correct
shapes and dtypes, call `openCOSMORS`, and return only the arrays needed by
MATLAB.

```python
# file: matlab_openCOSMORS_bridge.py
def run_case(component_paths, concentrations, temperature):
    import numpy as np
    import openCOSMORS

    # Build options, parameters, and calculation dictionaries here.
    # Use bindings/run_example.py as the template.
    return result
```

Call the wrapper from MATLAB:

```matlab
bridge = py.importlib.import_module("matlab_openCOSMORS_bridge");
result = bridge.run_case(componentPaths, concentrations, 298.15);
```

### MATLAB array-layout caveat

MATLAB stores numeric arrays in column-major order. C/C++ row-major arrays and
default NumPy arrays use the opposite memory layout. This matters because
`openCOSMORS.loadCalculations()` maps NumPy output arrays directly as row-major
Eigen views.

Safest pattern: convert MATLAB-originated arrays inside the Python wrapper.

```python
import numpy as np

def as_openCOSMORS_array(value, shape=None):
    arr = np.ascontiguousarray(value, dtype=np.float64)
    if shape is not None:
        arr = arr.reshape(shape, order="C")
    return arr

concentrations = as_openCOSMORS_array(concentrations_from_matlab, (1, 3))
temperatures = as_openCOSMORS_array([temperature_from_matlab], (1,))

calculation = {
    "concentrations": concentrations,
    "temperatures": temperatures,
    "reference_state_concentrations": as_openCOSMORS_array(np.empty((1, 0)), (1, 0)),
    "ln_gamma_x_SR_combinatorial_calc": np.zeros((1, 3), dtype=np.float64, order="C"),
    "ln_gamma_x_SR_residual_calc": np.zeros((1, 3), dtype=np.float64, order="C"),
    "ln_gamma_x_SR_calc": np.zeros((1, 3), dtype=np.float64, order="C"),
}
```

For higher-dimensional result arrays:

```python
contact_statistics = np.zeros((n_points, n_components, n_components),
                              dtype=np.float64, order="C")
average_surface_energies = np.zeros((n_points, n_interaction_matrices,
                                     n_components, n_components),
                                    dtype=np.float64, order="C")
```

Additional MATLAB notes:

- Call `pyenv` before the first `py.` command.
- If Python is already loaded, change environments only after `terminate(pyenv)`
  for out-of-process Python, or after restarting MATLAB for in-process Python.
- The Python extension must match the Python executable, operating system, CPU
  architecture, and Python ABI that MATLAB loads.
- MATLAB can convert numeric MATLAB arrays to NumPy arrays when NumPy is
  available, but conversion behavior differs by MATLAB release. A Python wrapper
  keeps this repo's direct NumPy binding requirements explicit.
- If Python returns a NumPy array, MATLAB can display it as a Python object and
  convert it with numeric conversion functions such as `double(...)`.

## Related projects

- [openCOSMO-RS_py](https://github.com/TUHH-TVT/openCOSMO-RS_py)
- [LVPP sigma profile database](https://github.com/lvpp/sigma)
- [Benchmark COSMO-SAC implementation](https://github.com/usnistgov/COSMOSAC)
- [Pysac](https://github.com/lvpp/pysac)
