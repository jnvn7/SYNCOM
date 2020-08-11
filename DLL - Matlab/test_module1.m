%--------------------------------------------------------------------%
% SYNCOM - A Nonlinear Synthetic Rope Numerical Computation Software
%
% Copyright 2020 Jessica Nguyen <nvnguyen@umass.edu>
%
% Module 1: Compute nonlinear strain (deformation) 
%           with provided stress history.
%
% Matlab script: template for using the dynamic link library SYNCOM.dll
%                to evaluate stress/strain constitutive relation.
%
% Input: stress time history
% Output: nonlinear total strain, visco-elastic strain, 
%         and visco-plastic strain.
%--------------------------------------------------------------------%

clear; clc;

%% Setup input (stress/strain);
module = 0;
input_file = 'Setting.xml';

%% Sinusoidal wave example. 
% %% Option 1 - Input from Matlab function;
% dt = 0.005;
% time = dt:dt:10;
% A = 0.05; T = 1; w = 2*pi/T;
% sigma = [0 0.05 0.1 -A*cos(w*time)+0.2];
% n = length(time);
% for i = 1:n
%     if (sigma(i) < 0)
%         sigma(i) = 0;
%     else
%     end
% end    

%% Option 2 - Read from text file;
fileID = fopen('inputData/SYNCOM_Input_mod1.txt', 'r');
fgetl(fileID);
sigma = fscanf(fileID, '%f');
fclose(fileID);

%% Creates data respository;
eps = zeros(1,length(sigma));
eps_ve = eps; eps_vp = eps;

%% Load *dll, solve, and extract results;
% Load library
loadlibrary('SYNCOM_API','SYNCOM_API.h');

% Initialization.
calllib('SYNCOM_API','initialize', module, input_file);

% Call solver
for i = 2:length(sigma)
    calllib('SYNCOM_API','SYNCOM', sigma(i));
    eps(i) = calllib('SYNCOM_API','extract_eps');
    eps_ve(i) = calllib('SYNCOM_API','extract_eps_ve');
    eps_vp(i) = calllib('SYNCOM_API','extract_eps_vp');
end

%% Unload library.
unloadlibrary('SYNCOM_API');



