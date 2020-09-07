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
fileID = fopen('inputData/SYNCOM_Input_mod2hs.txt', 'r');
fgetl(fileID);
inputdata = fscanf(fileID, '%f', [2 inf]);
time = inputdata(1,:);
eps = inputdata(2,:);
fclose(fileID);

% Create dt vector from the time history;
dt = zeros(1,length(time));
for i = 2:length(dt)
    dt(i) = time(i) - time(i-1);
end

%% Creates data respository;
sigma_cal = zeros(1,length(eps));
eps_ve = eps; eps_vp = eps;

%% Load *dll, solve, and extract results;
% Load library
loadlibrary('SynCOM_API','SynCOM_API.h');

% Initialization.
calllib('SynCOM_API','initializeSC', module, input_file);
err = 0; 
% Call solver
tic
for i = 2:length(eps)
    if (err == 0)
        err = calllib('SynCOM_API','SynCOM', eps(i), dt(i));
        eps_ve(i) = calllib('SynCOM_API','extract_eps_ve');
        eps_vp(i) = calllib('SynCOM_API','extract_eps_vp');
        sigma_cal(i) = calllib('SynCOM_API','extract_sigma');
    else
        unloadlibrary('SynCOM_API');
        err('Error detected. Check SynCOM_Log.txt for details.!');
    end
end
toc

%% Unload library.
calllib('SynCOM_API','close_app');
unloadlibrary('SynCOM_API');

%% Plots results;
plot(time, sigma_cal);



