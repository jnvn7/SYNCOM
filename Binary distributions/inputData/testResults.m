%% Read data from text file;
fileID = fopen('SYNCOM_Output_mod1.csv', 'r');
fgetl(fileID);
inputdata = fscanf(fileID, '%f', [5 inf]);
time1 = inputdata(1,:);
sigma = inputdata(2,:);
eps_cal = inputdata(3,:);
fclose(fileID);

%% Read data from text file;
fileID = fopen('SYNCOM_Output_mod2.csv', 'r');
fgetl(fileID);
inputdata = fscanf(fileID, '%f', [5 inf]);
time2 = inputdata(1,:);
sigma_cal = inputdata(2,:);
eps = inputdata(3,:);
fclose(fileID);

%% Plot comparison;
figure(1)
plot(time1, sigma, time2, sigma_cal,'*');

figure(2)
plot(time1, eps_cal, time2, eps, '*');
