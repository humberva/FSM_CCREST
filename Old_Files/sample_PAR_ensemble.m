%Sampling parameter estimates from basins in CONUS
clear;
clc;
load('USbasins_HydroParameters.mat');
load('USbasins_MorphParameters.mat');

%Filter out "bad" basins: Affected by Water Pixels - Need to recompute
%basin averaged parameter values
%Add contraints based on selected parameters
good_idx = find(uscrestpar(:,2) > 0 & usbasinlags(:,2) > 4 & usbasinlags(:,2) < 14 & uscrestpar(:,1) < 800);

uscrestpar = uscrestpar(good_idx,:);
usbasinmorph = usbasinmorph(good_idx,:);
usbasinlags = usbasinlags(good_idx,:);
uskwpars = uskwpars(good_idx,:);

% PKE - Multiplier to convert PET
% rundata.Parameters(1) = 1;
% PIM - The impervious area ratio. This ratio is stored as a percentage between 0 and 100.
% rundata.Parameters(2) = uscrestpar(idx,4).*100;
% PWM - The maximum soil water capacity (depth integrated pore space) of the soil layer. (mm)
% rundata.Parameters(3) = uscrestpar(idx,1);
% PFC - The soil saturated hydraulic conductivity (Ksat ,mm/hr).
% rundata.Parameters(4) = 1./uscrestpar(idx,2); %(mm/hr)^-1
% LEAKO - The overland reservoir discharge multiplier.
% rundata.Parameters(5) = min(1,exp(-1.0582*log(usbasinlags(idx,2)*10.8448)+0.8384));
% LEAKI - The interflow reservoir discharge multiplier.
% rundata.Parameters(6) = (uscrestpar(idx,2)*3/100)/(usbasinmorph(idx,6)*usbasinmorph(idx,5));
% PB - The exponent of the variable infiltration curve.
% rundata.Parameters(7) = uscrestpar(idx,3);
% Nqo - NASH Routing Parameters - Number of quick linear reservoir tanks
% rundata.Parameters(8) = 3;

%Random Sample
shuffled_sets = randperm(numel(good_idx));

sample_idxs = shuffled_sets(1:1000);

figure;
usmapgeo(usbasinmorph(:,3),usbasinmorph(:,4), 'k',20,[],'quan', 'true');
hold all;
scatterm(usbasinmorph(sample_idxs,3),usbasinmorph(sample_idxs,4), 'r', 'filled');

figure;
subplot(2,3,1);
[counts,bins] = hist(uscrestpar(sample_idxs,4).*100,20);
bar(bins,counts, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'BarWidth', 1);
set(gca, 'FontSize', 12);
xlabel('$\alpha_1$', 'interpreter', 'latex', 'fontsize',20);
ylabel('Number of Occurrences', 'FontSize', 14);

subplot(2,3,2);
[counts,bins] = hist(uscrestpar(sample_idxs,1),20);
bar(bins,counts, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'BarWidth', 1);
set(gca, 'FontSize', 12);
xlabel('$\alpha_2$', 'interpreter', 'latex', 'fontsize',20);
ylabel('Number of Occurrences', 'FontSize', 14);

subplot(2,3,3);
[counts,bins] = hist(1./uscrestpar(sample_idxs,2),20);
bar(bins,counts, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'BarWidth', 1);
set(gca, 'FontSize', 12);
xlabel('$\alpha_3$', 'interpreter', 'latex', 'fontsize',20);
ylabel('Number of Occurrences', 'FontSize', 14);

subplot(2,3,4);
leako = exp(-1.0582*log(usbasinlags(sample_idxs,2)*10.8448)+0.8384);
leako(leako>1) = 1;
[counts,bins] = hist(leako,20);
bar(bins,counts, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'BarWidth', 1);
set(gca, 'FontSize', 12);
xlabel('$\alpha_4$', 'interpreter', 'latex', 'fontsize',20);
ylabel('Number of Occurrences', 'FontSize', 14);

subplot(2,3,5);
[counts,bins] = hist((uscrestpar(sample_idxs,2)*3/100)./(usbasinmorph(sample_idxs,6).*usbasinmorph(sample_idxs,5)),20);
bar(bins,counts, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'BarWidth', 1);
set(gca, 'FontSize', 12);
xlabel('$\alpha_5$', 'interpreter', 'latex', 'fontsize',20);
ylabel('Number of Occurrences', 'FontSize', 14);

subplot(2,3,6);
[counts,bins] = hist(uscrestpar(sample_idxs,3),20);
bar(bins,counts, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'BarWidth', 1);
set(gca, 'FontSize', 12);
xlabel('$\alpha_6$', 'interpreter', 'latex', 'fontsize',20);
ylabel('Number of Occurrences', 'FontSize', 14);

% save('Random_Selection_of_Parameter_Estimates_CONUS.mat', 'sample_idxs', 'uscrestpar', 'usbasinmorph', 'usbasinlags', 'uskwpars');

%% PCA Analysis of Parameters' Ensembles

%% Run ensemble test
% load('Random_Selection_of_Parameter_Estimates_CONUS.mat');
load('Behavioral_Constraint_Random_Selection_of_Parameter_Estimates_CONUS.mat');
ensPIM = uscrestpar(sample_idxs,4).*100;
ensPWM = uscrestpar(sample_idxs,1);
ensPFC = 1./uscrestpar(sample_idxs,2);
%Original model
ensLEAKO = exp(-1.0582*log(usbasinlags(sample_idxs,2)*10.8448)+0.8384);
%Modified subjectively
% ensLEAKO = exp(-1.0582*log(usbasinlags(sample_idxs,2)*100.8448)+0.8384);
% ensLEAKO(ensLEAKO>1) = 1;
ensLEAKI = (uscrestpar(sample_idxs,2)*3/100)./(usbasinmorph(sample_idxs,6).*usbasinmorph(sample_idxs,5));
ensPB = uscrestpar(sample_idxs,3);

%Estimated peak of rain storm based on Analytical Basin Lag Time
rainPeakdate = datenum('19-Jun-2006 05:56:05');

cont = 1;
hr = 0; %:12:444

%Input Data: Tar River basin
dd = dlmread('/Users/humbertovergara/Documents/Trabajo/Parameter_and_QPEs_Error_Characterization_for_DA/TS_Data/MPE_tarboro_basin_avg_time_series.csv', ',', 1,1);
rundata.PP = dd(:,2); %load('~/Documents/MS_Work/Thesis_Work/HyMOD/Rain_MPE_200206-200912.txt'); 
rundata.PP_period = datenum('01-Jan-2002 01:00:00'):1/24:datenum('31-Dec-2009 23:00:00');
rundata.PET = dlmread('~/Documents/MS_Work/Thesis_Work/HyMOD/PET_mmhr_Tarboro_2002-2009.txt', ' ', 1,0);  
rundata.PET_period = datenum('01-Jan-2002 00:00:00'):1/24:datenum('31-Dec-2009 23:00:00');

%Real Observations
rundata.obsQ  = dlmread('~/Documents/MS_Work/Thesis_Work/HyMOD/02083500_ObsQ_2002-2010_hrly.csv', ',', 1,1); %5709.28 km2 Tarboro drainage area  
rundata.obsQ_period = datenum('01-Jan-2002 00:00:00'):1/24:datenum('30-Jun-2010 21:00:00');
rundata.barea = 5709;  %Basin Area in km^2

%Simulation Data and Settings
rundata.tstep = 1; %Time Step in hours
% rundata.tperiod = st_date:rundata.tstep/24:datenum(2006,6,19,17,0,0)+20;
rundata.tperiod = datenum(2006,6,6,0,0,0):rundata.tstep/24:datenum(2006,6,25,23,0,0);
% rundata.tperiod = datenum(2002,1,1,1,0,0):rundata.tstep/24:datenum(2009,12,31,23,0,0);
% rundata.tperiod = datenum(2003,1,1,1,0,0):rundata.tstep/24:datenum(2006,10,1,0,0,0);


%Model Settings
rundata.runoff_model = 'crest';
rundata.class = 'Lumped';
rundata.mode = 'Deterministic';
caldata = load('Calibration_started_20130827214217_finished_20130901005456.mat');

figure;

tic;
for idx = 1:numel(sample_idxs)
    st_date = datenum(2006,6,6,0+hr,0,0);

rundata.Parameters(1) = 1;
rundata.Parameters(2) = ensPIM(idx);
rundata.Parameters(3) = ensPWM(idx);
rundata.Parameters(4) = ensPFC(idx);
rundata.Parameters(5) = ensLEAKO(idx);
rundata.Parameters(6) = ensLEAKI(idx);
rundata.Parameters(7) = ensPB(idx);
rundata.Parameters(8) = 3;

% fprintf('Area = %f\nLeako = %f\nLeaki = %f\n\n', usbasinmorph(idx,2), rundata.Parameters(5:6));

%Approximation Parameters
alpha.alpha1 = 10^2;
beta.beta1 = 10^1;
alpha.alpha2 = alpha.alpha1;
beta.beta2 = beta.beta1;
alpha.alpha3 = alpha.alpha1;
beta.beta3 = beta.beta1;
alpha.alpha5 = alpha.alpha1;
beta.beta5 = beta.beta1;
alpha.alpha6 = alpha.alpha1;
beta.beta6 = beta.beta1;
alpha.alpha7 = alpha.alpha1; %sensitive when 10^0
beta.beta7 = beta.beta1;
alpha.alpha8 = alpha.alpha1; %sensitive
beta.beta8 = beta.beta1;
alpha.alpha9 = alpha.alpha1; %sensitive
beta.beta9 = beta.beta1;
alpha.alpha10 = alpha.alpha1;
beta.beta10 = beta.beta1;
alpha.alpha11 = alpha.alpha1;
beta.beta11 = beta.beta1;
alpha.alpha12 = alpha.alpha1;
beta.beta12 = beta.beta1;

rundata.alpha = alpha;
rundata.beta = beta;

%CREST States
States_0.IWU = 0;
States_0.ISO = [0 0 0 0 0];
States_0.ISU = 0;
rundata.States_0 = States_0;

%Subset data to period of interest
[new_rundata] = subset(rundata,0);
rundata = new_rundata;

%Run CCREST
States_0.IWU = 10; %Stored_States.SM(hr+1);
States_0.ISU = 0.1; %Stored_States.ISU(hr+1);
States_0.ISO = [0.5 0.2 0.2]; %Stored_States.ISO(:,hr+1);
rundata.States = States_0;
% ref_prediction = HPRO(rundata);
prediction = ccrest(rundata);
plot(prediction.Q, 'Color', [0.7 0.7 0.7]);
hold all;
end
plot(rundata.obsQ, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
toc;