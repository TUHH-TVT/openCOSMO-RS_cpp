
# openCOSMO-RS

This is an open source implementation of the COSMO-RS model that uses multiple descriptor sigma profiles. The corresponding publication to this repository is:

#### Clone it 

After cloning this repository to build you will need to get the respective submodules running the following commands inside the cloned repository folder:
- git submodule init
- git submodule update
***
#### Build for python
> Build on **Windows** with Visual Studio (tested with vs community 2019).
> - Please check if you have set the correct directories for the openCOSMORS project. Right now it assumes standard Anaconda installation on Windows.
>   - Python include directory
>     - Properties > Configuration Properties > C/C++ > Additional Include Directories
>   - Python library directory
>     - Properties > Configuration Properties > Linker > Additional Library Directories
>   - [More information here](https://docs.microsoft.com/en-us/visualstudio/python/working-with-c-cpp-python-in-visual-studio?view=vs-2019)
> - Please check whether the specified instruction set matches your target hardware
>     - Properties > Configuration Properties > C/C++ > Code Generation > Enable Enhanced Instruction Set
>     - Supported are AVX, AVX2, if you select any other, SSE3 will be used.
> - A small python project is also included showing how to use it from within visual studio.

> Build on **Linux** with gcc (tested with gcc version 9.2.0 and python 3.6):
> - Execute the following command within the bindings folder after adding the corrects paths. (_pybind11_include_folder_, _eigen_include_folder_, _python_include_folder_)
> - Specify the _parallelization_flag_ [-msse3, -mavx, -mfma]:
> - g++ -fopenmp _parallelization_flag_ -O3 -Wall -shared -std=c++14 -fPIC \`python3 -m pybind11 --includes\` bindings_forPython.cpp -o openCOSMORS\`python3-config --extension-suffix\` -I _pybind11_include_folder_ -I _eigen_include_folder_ -I _python_include_folder_
>- Example:
>> g++ -fopenmp -mavx -O3 -Wall -shared -std=c++14 -fPIC \`python3 -m pybind11 --includes\` ../code/bindings_forPython.cpp -o openCOSMORS\` python3-config --extension-suffix\` -I ../pybind11/include -I ../eigen -I /usr/include/python3.6m

***
#### Build for MATLAB
> Build on **Windows** and **Linux** (tested only on windows with MSVC and MATLAB 2019a):
> - Ensure your matlab include folder is set correctly in the file _code/bindings_forMATLAB.cpp_:
>    - #include "_matlabroot_\extern\include\mex.hpp"
>    - #include "_matlabroot_\\extern\include\mexAdapter.hpp"
> - _matlabroot_ can be found by executing said command in MATLAB
> - Execute the following file within the bindings folder: _compile_mex.m_
***
#### Running
> An exemplary file to run the model on python and matlab is included in the bindings folder.
