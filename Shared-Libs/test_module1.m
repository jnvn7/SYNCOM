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
fileID = fopen('inputData/SynCOM_Input_mod1hs.txt', 'r');
fgetl(fileID);
inputdata = fscanf(fileID, '%f', [2 inf]);
time = inputdata(1,:);
sigma = inputdata(2,:);
fclose(fileID);

% Create dt vector from the time history;
dt = zeros(1,length(time));
for i = 2:length(dt)
    dt(i) = time(i) - time(i-1);
end

%% Creates data respository;
eps = zeros(1,length(sigma));
eps_ve = eps; eps_vp = eps;

%% Load *dll, solve, and extract results;
% Load library
loadlibrary('SynCOM_API','SynCOM_API.h');

% Initialization.
calllib('SynCOM_API','initializeSC', module, input_file);
err = 0; 
% Call solver
tic
for i = 2:length(sigma)
    if (err == 0)
        err = calllib('SynCOM_API','SynCOM', sigma(i), dt(i));
        eps(i) = calllib('SynCOM_API','extract_eps');
        eps_ve(i) = calllib('SynCOM_API','extract_eps_ve');
        eps_vp(i) = calllib('SynCOM_API','extract_eps_vp');
    else
        unloadlibrary('SynCOM_API');
        err('Error detected. Check SynCOM_Log.txt for details.!');
    end
end
toc

%% Close the application and unload library.
calllib('SynCOM_API','close_app');
unloadlibrary('SynCOM_API');

%% Plots results;
plot(time, eps);



