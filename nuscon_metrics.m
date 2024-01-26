% 1. install python and conda (or miniconda);
% 2. create conda environment: conda create --name nuscon python=3.9
% 3. activate this environment: conda activate nuscon
% 4. install corresponding python packages: pip install numpy tqdm networkx
% 5. set pyenv in Matlab, for example: pyenv(Version="C:\Users\14291\miniconda3\envs\nuscon\python.exe")
% more information about pyenv setting, please reference to "https://www.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html"

clear; clc; close all;
% addpath('data\nuscon\Ny_Nz')
% addpath('data\nuscon\quinine')
% addpath('data\nuscon\azithromycin')
addpath('data\nuscon\quinine')

% rec_result = load('Recsoft_Ny_Nz.mat'); 
% rec_result = load('Rechard_Ny_Nz.mat'); 
% rec_result = load('Rechard_azithromycin.mat'); 
% rec_result = load('Recsoft_azithromycin.mat');
rec_result = load('Recsympoisson_quinine.mat');
% rec_screen = rec_result.Rec_SCREEN;
rec_screenes = rec_result.Rec_SCREENES;
ideal_spec = rec_result.Spec2D_Ideal;

% peaks_position = load('psoft_Ny_Nz.mat');
% peaks_position = load('phard_Ny_Nz.mat');
% peaks_position = load('phard_azithromycin.mat');
% peaks_position = load('psoft_azithromycin.mat');
peaks_position = load('psympoisson_quinine.mat');
% position_screen = peaks_position.locSCREEN;
position_screenes = peaks_position.locSCREENES;
position_ideal = peaks_position.locIdeal;

mod = py.importlib.import_module('utils');
py.importlib.reload(mod);

% The order of the result array: freq_accuracy, intensity_linearity, true_positive_rate, false_positive_rate 
% results_screenes = double(py.utils.nuscon_metrics(py.numpy.array(ideal_spec), py.numpy.array(rec_screenes), ...
%     py.numpy.array(position_screenes), py.numpy.array(position_ideal), 1));
% results_screen = double(py.utils.nuscon_metrics(py.numpy.array(ideal_spec), py.numpy.array(rec_screen), ...
%     py.numpy.array(position_screen), py.numpy.array(position_ideal), 1));

%% quinine
results_screenes = double(py.utils.nuscon_metrics(py.numpy.array(ideal_spec), py.numpy.array(rec_screenes), ...
    py.numpy.array(position_screenes), py.numpy.array(position_ideal), 1));
