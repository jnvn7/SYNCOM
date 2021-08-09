clear; clc;

%% Setup;
X = zeros(6,1); % platform position
XD = zeros(6,1); % platform velocity
TEnd = 8000; % simulation time in seconds
dt = 0.005; % coupling time step size (time between MoorDyn calls)
N = TEnd/dt;
Ts = 0:dt:(N-1)*dt; % time step array
FairTens1 = zeros(N,1); % array for storing fairlead 1 tension time series
FLines_temp = zeros(1,6); % going to make a pointer so LinesCalc can modify FLines
FLines_p = libpointer('doublePtr', FLines_temp); % access returned value with FLines_p.value

% Nodal position pointers;
numLine = 1;
numNode = 20;  
nodePOS = zeros(numNode, 3);
nodePOS_ptr = libpointer('doublePtr', nodePOS(1, :));

% Pointer to store fairleads' strains;
epsFL = zeros(N, 1);
epsFL_ptr = libpointer('doublePtr', epsFL(1, :));

d = 0.193;
A = d*d*pi/4;
MBL = 4.1647e+08*A;

%% Initialization
%tic
loadlibrary('MoorDynSC','MoorDynSC.h'); % load MoorDyn DLL
calllib('MoorDynSC','LinesInit',X,XD) % initialize MoorDyn
FairTens1(1) = calllib('MoorDynSC','GetFairTen',1);
 
%% Simulation
ti = 1; sigma = [0.1 0.35 0];
XD(1) = 3.5; % Fairlead speed for ramping up or down to a certain stress level;

% Specification for the case of sigma = 0. If this is set, fairlead oscillates
% around the current stress level;
Amp2 = 1; T = 10; freq = 2*pi/T;
cycles = 20; 

for j = 1:length(sigma) 
    %Type 1: ramp up to a specified stress level;
    if (j == 1 || (sigma(j) >= sigma(j-1)))
        while (FairTens1(ti)/MBL < sigma(j))
            X = X + XD*dt;
            calllib('MoorDynSC', 'LinesCalc', X, XD, FLines_p, Ts(ti), dt);
            FairTens1(ti+1) = calllib('MoorDynSC','GetFairTen',1);
            
            calllib('MoorDynSC','SC_GetEpsFL', numLine, epsFL_ptr);
            epsFL(ti+1,:) = epsFL_ptr.value;
            
            ti = ti + 1;
        end
    %Type 2: decrease to a specified stress level;
    elseif (sigma(j) < sigma(j-1) && sigma(j)~= 0)
        while (FairTens1(ti)/MBL > sigma(j))
            X = X - XD*dt;
            calllib('MoorDynSC', 'LinesCalc', X, XD, FLines_p, Ts(ti), dt);
            FairTens1(ti+1) = calllib('MoorDynSC','GetFairTen',1);
            
            calllib('MoorDynSC','SC_GetEpsFL', numLine, epsFL_ptr);
            epsFL(ti+1,:) = epsFL_ptr.value;
            
            ti = ti + 1;
        end
    elseif(sigma(j) == 0)
        %Type 3: oscillate around stress levels for a specified number of cycles;
        t2 = Ts(ti);
        for i = 1:ceil(T*cycles/dt)
            XD(1) = -Amp2*freq*cos(freq*(Ts(ti)- t2)); X = X + XD*dt; 
            calllib('MoorDynSC', 'LinesCalc', X, XD, FLines_p, Ts(ti), dt);
            FairTens1(ti+1) = calllib('MoorDynSC','GetFairTen',1); %

            calllib('MoorDynSC','SC_GetEpsFL', numLine, epsFL_ptr);
            epsFL(ti+1,:) = epsFL_ptr.value;
            ti = ti + 1;
        end
    else
        disp('No scenario selected');
    end
end
%toc

%% Plotting results;
plot(Ts(1:ti), FairTens1(1:ti)/MBL);
xlabel('Time (s)');
ylabel('Fairlead stress (%MBL)');

calllib('MoorDynSC','LinesClose'); % close MoorDyn
unloadlibrary MoorDynSC; % unload library (never forget to do this!)








