%--------------------------------------------------------------------%
% SYNCOM - A Nonlinear Synthetic Rope Numerical Computation Software
%
% Copyright 2020 Jessica Nguyen <nvnguyen@umass.edu>
%
% Module 2: Inversed constitutive model. Computes applying stress from
%           the provided strain time history.
%
% Matlab script: template for using the dynamic link library SYNCOM.dll
%                to evaluate stress/strain constitutive relation.
%
% Input: total strain time history
% Output: applying stress, nonlinear visco-elastic strain, 
%         and visco-plastic strain.
%--------------------------------------------------------------------%

clear; clc;

%% Setup input (stress/strain);
module = 1;
input_file = 'Setting.xml';

%% Read strain input from text file;
fileID = fopen('inputData/SYNCOM_Input_epsilon.txt', 'r');
fgetl(fileID);
eps = fscanf(fileID, '%f');
fclose(fileID);

%% Creates data respository;
sigma_cal = zeros(1,length(eps));
eps_ve = eps; eps_vp = eps;

%% Load *dll, solve, and extract results;
% Load library
loadlibrary('SYNCOM_DLL','SYNCOM_API.h');

% Initialization.
calllib('SYNCOM_DLL','initialize', module, input_file);

% Call solver
for i = 2:length(eps)
    calllib('SYNCOM_DLL','SYNCOM', eps(i));
    eps_ve(i) = calllib('SYNCOM_DLL','extract_eps_ve');
    eps_vp(i) = calllib('SYNCOM_DLL','extract_eps_vp');
    sigma_cal(i) = calllib('SYNCOM_DLL','extract_sigma');
end

%% Unload library.
unloadlibrary('SYNCOM_DLL');





