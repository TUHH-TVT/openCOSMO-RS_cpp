%{
    c++ implementation of openCOSMO-RS including multiple segment descriptors
    @authors: Simon Mueller & Andrés González de Castilla, 2022
%}

clear mex mh;

output_dir = '../bindings';
cd('../code');

srcFile = 'bindings_forMATLAB.cpp';
mexfilename = 'openCOSMORS';

%% compiler switch options
compiler_switches = { '-D__AVX__'}; % vectorization level {'-D__SSE3__', '-D__AVX__', '-D__FMA__'}
compiler_switches(end + 1) = {['-I' fullfile(pwd,'../eigen')]}; % eigen library source code files
% compiler_switches(end + 1) = {'COMPFLAGS="$COMPFLAGS /openmp"'};% use openmp multi-threaded calculation
% compiler_switches(end + 1) = {'-g'}; % compile with debug flag
% compiler_switches(end + 1) = {'-DMEASURE_TIME'}; % other debug options {'-DMEASURE_TIME', '-DDEBUG_INFO'}

mex(compiler_switches{:}, '-output',  mexfilename, srcFile);

if exist(output_dir, 'dir') ~= 7
    mkdir(output_dir)
end

% optional renaming
mexfiles = dir([mexfilename '.mex*']);
for i = 1:length(mexfiles)
    movefile(mexfiles(i).name, fullfile(pwd, output_dir, mexfiles(i).name));
end

cd(output_dir)

%% To debug MEX from VS
%  - Compile MEX with compiler flag '-g'
%  - Hit Ctrl + Alt + P in VS
%  - Attach to MATLAB mex host process
%  - Execute MEX from within MATLAB

%  as an example you can do the following:
%  - execute "mh = mexhost";
%  - attach VS to the mexhost process as described above
%  - execute "run_example(mh);"