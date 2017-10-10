clear; 
clc;

addpath('./FSM_Functions/');
addpath('~/Documents/Trabajo/HPro_Tool/Code/hpro/OldVersions_of_Functions/');

% load('H_Isabel_FSE_Analysis_Stored_States.mat'); %24-Sep-2003 12:00:00
% load('USbasins_HydroParameters.mat');
% load('USbasins_MorphParameters.mat');
% tarbasin_idx = 1319; %1320;
% load('TS_Alberto_FSE_Analysis_Stored_States.mat');
% 
% %Pre-selected basins
% bidx = find(uscrestpar(:,2) > 0); % & uscrestpar(:,2) < 5.5);
% uscrestpar = uscrestpar(bidx,:);
% usbasinlags = usbasinlags(bidx,:);
% usbasinmorph = usbasinmorph(bidx,:);

%% Randomly selected basins
load('Behavioral_Constraint_Random_Selection_of_Parameter_Estimates_CONUS.mat');
ensPIM = uscrestpar(sample_idxs,4).*100;
ensPWM = uscrestpar(sample_idxs,1);
ensPFC = 1./uscrestpar(sample_idxs,2);
ensLEAKO = exp(-1.0582*log(usbasinlags(sample_idxs,2)*10.8448)+0.8384);
ensLEAKO(ensLEAKO>1) = 1;
ensLEAKI = (uscrestpar(sample_idxs,2)*3/100)./(usbasinmorph(sample_idxs,6).*usbasinmorph(sample_idxs,5));
ensPB = uscrestpar(sample_idxs,3);


%Estimated peak of rain storm based on Analytical Basin Lag Time
rainPeakdate = datenum('19-Jun-2006 05:56:05');
%Calibrated leako based on TS Alberto and peak timing
% rundata.Parameters(1) = 1;
% rundata.Parameters(2) = uscrestpar(tarbasin_idx,4).*100;
% rundata.Parameters(3) = uscrestpar(tarbasin_idx,1);
% rundata.Parameters(4) = uscrestpar(tarbasin_idx,2)/10; %cm/hr
% rundata.Parameters(5) = 0.0157; %leako_sample;
% rundata.Parameters(6) = (uscrestpar(tarbasin_idx,2)*3/100)/(usbasinmorph(tarbasin_idx,6)*usbasinmorph(tarbasin_idx,5));
% rundata.Parameters(7) = uscrestpar(tarbasin_idx,3);
% rundata.Parameters(8) = 3;

Optimal_leako = 0.0157;
leako_sample = 0.005:0.005:0.1;
cont = 1;
% for 
hr = 0; %:12:444
randomized_par_sample = randperm(size(usbasinlags,1));
% for idx = tarbasin_idx %randomized_par_sample(1:3)
% for idx = 1:20
% for leako_sample = 1 %0.005:0.005:0.2 %1
    st_date = datenum(2006,6,6,0+hr,0,0);
%Input Data: Tar River basin
dd = dlmread('/Users/humbertovergara/Documents/Trabajo/Parameter_and_QPEs_Error_Characterization_for_DA/TS_Data/MPE_tarboro_basin_avg_time_series.csv', ',', 1,1);
rundata.PP = dd(:,2); %load('~/Documents/MS_Work/Thesis_Work/HyMOD/Rain_MPE_200206-200912.txt'); 
rundata.PP_period = datenum('01-Jan-2002 01:00:00'):1/24:datenum('31-Dec-2009 23:00:00');
rundata.PET = dlmread('~/Documents/MS_Work/Thesis_Work/HyMOD/PET_mmhr_Tarboro_2002-2009.txt', ' ', 1,0);  
rundata.PET_period = datenum('01-Jan-2002 00:00:00'):1/24:datenum('31-Dec-2009 23:00:00');

%Real Observations
rundata.obsQ  = dlmread('~/Documents/MS_Work/Thesis_Work/HyMOD/02083500_ObsQ_2002-2010_hrly.csv', ',', 1,1); %5709.28 km2 Tarboro drainage area  
rundata.obsQ_period = datenum('01-Jan-2002 00:00:00'):1/24:datenum('30-Jun-2010 21:00:00');

%Synthetic Observations
% load('Synthetic_Observations_CCREST_Simulation_2002-2009.mat');
% rundata.obsQ = synth_obs.Q';
% rundata.obsQ_period = synth_obs.period;
rundata.barea = 5709;  %Basin Area in km^2

%Simulation Data and Settings
rundata.tstep = 1; %Time Step in hours
% rundata.tperiod = st_date:rundata.tstep/24:datenum(2006,6,19,17,0,0)+20;
rundata.tperiod = datenum(2006,6,6,0,0,0):rundata.tstep/24:datenum(2006,6,25,23,0,0);
% rundata.tperiod = datenum(2002,1,1,1,0,0):rundata.tstep/24:datenum(2009,12,31,23,0,0);
% rundata.tperiod = datenum(2003,1,1,1,0,0):rundata.tstep/24:datenum(2006,10,1,0,0,0);

%Steady-State Analysis
steady_state.activate = 0;
steady_state.mode = 'unit';
if (steady_state.activate == 1)
    rundata.tperiod = st_date:rundata.tstep/24:datenum(2006,12,20,0,0,0)+10;
    switch steady_state.mode
        case 'dying'
            %T.S. Alberto
            rundata.PP(rundata.PP_period >= rundata.tperiod(225)) = 0;
        case 'flat'
            %T.S. Alberto
            rundata.PP(rundata.PP_period >= rundata.tperiod(210)) = max(rundata.PP(rundata.PP_period >= rundata.tperiod(1) & rundata.PP_period <= rundata.tperiod(end)));
        case 'unit'
            %T.S. Alberto
            rundata.PP(rundata.PP_period >= rundata.tperiod(190)) = 1;
    end
end

%Model Settings
rundata.runoff_model = 'crest';
rundata.class = 'Lumped';
rundata.mode = 'Deterministic';
caldata = load('Calibration_started_20130827214217_finished_20130901005456.mat');
% allparsets_order = randperm(100000);
% Parameters.PKE = 1; %Multiplier to convert PET
% Parameters.PIM = 30; %The impervious area ratio. This ratio is stored as a percentage between 0 and 100.
% Parameters.PWM = 100; %The maximum soil water capacity (depth integrated pore space) of the soil layer. (mm)
% Parameters.PFC = 0.5; %The soil saturated hydraulic conductivity (Ksat ,mm/hr).
% Parameters.LEAKO = 0.025; %The overland reservoir discharge multiplier.
% Parameters.LEAKI = 0.0001; %The interflow reservoir discharge multiplier.
% Parameters.PB = 100; %The exponent of the variable infiltration curve.
% 
% %NASH Routing Parameters
% Parameters.Nqo = 3;

%Optimal Set
rundata.Parameters = caldata.prediction.BestPars;

% %Random Set
% valid_set = caldata.prediction.EvalPars(caldata.prediction.EvalPars(:,10) > prctile(caldata.prediction.EvalPars(:,10),10) & caldata.prediction.EvalPars(:,3) > 200 & caldata.prediction.EvalPars(:,7) < 100,:);
% % rundata.Parameters = valid_set(end,1:8);
% allparsets_order = randperm(size(valid_set,1));
% rundata.Parameters = valid_set(allparsets_order(1),1:8);
% % rundata.Parameters = caldata.prediction.EvalPars(25000,1:8);

%Manual perturbations
% rundata.Parameters(1) = 1;
% rundata.Parameters(2) = uscrestpar(tarbasin_idx,4).*100;
% rundata.Parameters(3) = uscrestpar(tarbasin_idx,1);
% rundata.Parameters(4) = 1./uscrestpar(tarbasin_idx,2); %(mm/hr)^-1
% rundata.Parameters(5) = exp(-1.0582*log(usbasinlags(tarbasin_idx,2)*10.8448)+0.8384);
% rundata.Parameters(6) = (uscrestpar(tarbasin_idx,2)*3/100)/(usbasinmorph(tarbasin_idx,6)*usbasinmorph(tarbasin_idx,5));
% rundata.Parameters(7) = uscrestpar(tarbasin_idx,3);
% rundata.Parameters(8) = 3;

% for idx = 1:20
% 
% rundata.Parameters(1) = 1;
% rundata.Parameters(2) = uscrestpar(idx,4).*100;
% rundata.Parameters(3) = uscrestpar(idx,1);
% rundata.Parameters(4) = 1./uscrestpar(idx,2); %(mm/hr)^-1
% rundata.Parameters(5) = min(1,exp(-1.0582*log(usbasinlags(idx,2)*10.8448)+0.8384));
% rundata.Parameters(6) = (uscrestpar(idx,2)*3/100)/(usbasinmorph(idx,6)*usbasinmorph(idx,5));
% rundata.Parameters(7) = uscrestpar(idx,3);
% rundata.Parameters(8) = 3;

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

for idx = 1:20

rundata.Parameters(1) = 1;
rundata.Parameters(2) = uscrestpar(idx,4).*100;
rundata.Parameters(3) = uscrestpar(idx,1);
rundata.Parameters(4) = 1./uscrestpar(idx,2); %(mm/hr)^-1
rundata.Parameters(5) = min(1,exp(-1.0582*log(usbasinlags(idx,2)*10.8448)+0.8384));
rundata.Parameters(6) = (uscrestpar(idx,2)*3/100)/(usbasinmorph(idx,6)*usbasinmorph(idx,5));
rundata.Parameters(7) = uscrestpar(idx,3);
rundata.Parameters(8) = 3;

prediction = ccrest(rundata);

% Run Forward Sensitivity Method (FSM)
[correction, V, U, A, H, Dh, condnum,rankcoef] = fsmccrest(rundata,prediction,alpha,beta);

%Iterate with FSM
doitagain = 0;
while (doitagain == 1)
    close all; clc;
% correction = zeros(5,1);

%Make Correction and run CCREST
%New CREST States
States_0.IWU = States_0.IWU+(correction(1)*100/rundata.Parameters(3));
States_0.ISU = States_0.ISU+correction(2);
States_0.ISO = States_0.ISO+correction(3:5)';
new_Parameters = rundata.Parameters+[correction(6:end); 0]';

% %Poor man's Boundary Constraints Enforcement
% States_0.IWU = max(0,States_0.IWU+(correction(1)*100/rundata.Parameters(3)));
% States_0.ISU = max(0,States_0.ISU+correction(2));
% States_0.ISO = States_0.ISO+correction(3:5)';
% States_0.ISO(1) = max(0,States_0.ISO(1));
% States_0.ISO(2) = max(0,States_0.ISO(2));
% States_0.ISO(3) = max(0,States_0.ISO(3));
% 
% new_Parameters = rundata.Parameters+[correction(6:end); 0]';
% new_Parameters(1) = min(max(0,new_Parameters(1)),1);
% new_Parameters(2) = min(max(0,new_Parameters(2)),100);
% new_Parameters(3) = max(0,new_Parameters(3));
% new_Parameters(4) = max(0,new_Parameters(4));
% new_Parameters(5) = min(max(0,new_Parameters(5)),1);
% new_Parameters(6) = min(max(0,new_Parameters(6)),1);
% new_Parameters(7) = max(0,new_Parameters(7));

rundata.Parameters = new_Parameters;
rundata.States = States_0;
newprediction = ccrest(rundata);

% %Plot Results
% %Plot and Stats
st_idx = find(rundata.obsQ_period == rundata.tperiod(1));
en_idx = find(rundata.obsQ_period == rundata.tperiod(end));
rundata.obsQ = rundata.obsQ(st_idx:en_idx);
% [bias(1),rmse(1),nsce(1)] = nanhydrostat(rundata.obsQ',prediction.Q);
% [bias(2),rmse(2),nsce(2)] = nanhydrostat(rundata.obsQ',newprediction.Q);
% [bias(3),rmse(3),nsce(3)] = nanhydrostat(ref_prediction.Q,newprediction.Q);
% 
obsq_h = area(rundata.tperiod,rundata.obsQ, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7]); hold all;
% crest_h = plot(rundata.tperiod,ref_prediction.Q, 'DisplayName', 'CREST', 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
nofsm_h = plot(rundata.tperiod,prediction.Q, 'DisplayName', 'CCREST Before FSM', 'Color', 'b');
fsm_h = plot(rundata.tperiod,newprediction.Q, 'DisplayName', 'CCREST After FSM', 'Color', 'r');
set(gca, 'FontSize', 12);
xlabel('Time (Hours)', 'FontSize', 14); ylabel('Streamflow (m^3/s)', 'FontSize', 14);
datetick('x', 'dd/mm', 'keepticks', 'keeplimits');
% legend([obsq_h crest_h nofsm_h fsm_h], 'Observed Q', 'CREST', 'CCREST Before FSM', 'CCREST After FSM');

doitagain = input('Again?: ');
prediction = newprediction;
end

% %Plot Results
% %Plot and Stats
st_idx = find(rundata.obsQ_period == rundata.tperiod(1));
en_idx = find(rundata.obsQ_period == rundata.tperiod(end));
rundata.obsQ = rundata.obsQ(st_idx:en_idx);
% [bias(1),rmse(1),nsce(1)] = nanhydrostat(rundata.obsQ',prediction.Q);
% [bias(2),rmse(2),nsce(2)] = nanhydrostat(rundata.obsQ',newprediction.Q);
% [bias(3),rmse(3),nsce(3)] = nanhydrostat(ref_prediction.Q,newprediction.Q);
% 
figure(1)
% obsq_h = area(rundata.tperiod,rundata.obsQ, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7]); 
hold all;
% crest_h = plot(rundata.tperiod,ref_prediction.Q, 'DisplayName', 'CREST', 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');
nofsm_h = plot(rundata.tperiod,prediction.Q, 'DisplayName', 'CCREST'); %, 'Color', 'b'); %[0.7 0.7 0.7]); 
set(gca, 'FontSize', 12, 'Xlim', [rundata.tperiod(1) rundata.tperiod(end)]);
xlabel('Time (Hours)', 'FontSize', 14); ylabel('Streamflow (m^3/s)', 'FontSize', 14);
datetick('x', 'dd/mm', 'keepticks', 'keeplimits');
% close;

%Plot Sensitivity
figure(2);
cont = 1;
for i = 1:size(V,1)
    for j = 1:size(V,2)
        subplot(size(V,1),size(V,2),cont);
%         [ybs] = get(gca, 'Ylim');
%         plot(datenum(2006,6,25,23,0,0)*ones(1,length(ybs(1):0.05:ybs(2))), ybs(1):0.05:ybs(2), 'Color', 'k', 'LineStyle', '--')
        reV = reshape(V(i,j,:),1,size(V,3));
        plot(rundata.tperiod,reV(1:end-1)); %, 'Color', 'k', 'LineStyle', '-'); hold all;
        hold all;
        title(['V_k x_', num2str(i), '/x_', num2str(j), '_(_0_)'], 'FontSize', 12);
        set(gca, 'Xlim', [rundata.tperiod(1) rundata.tperiod(end)], 'Xtick', []);
%         set(gca, 'Xlim', [datenum(2006,6,6,0,0,0) datenum(2006,6,25,23,0,0)], 'Xtick', []);
%         set(gca, 'Xlim', [datenum(2006,6,6,0,0,0) datenum(2006,7,9,17,0,0)], 'Xtick', []);
%         set(gca, 'Xlim', [datenum(2006,6,6,0,0,0) datenum(2006,7,9,17,0,0)], 'Xtick', []);
%         xlabel('Time Period (hourly)', 'FontSize', 10);
        cont = cont + 1;
    end
end
% 
figure(3);
cont = 1;
% CCmat.pearson = zeros(size(U,1),size(U,2));
% CCmat.spearman = zeros(size(U,1),size(U,2));
% CCmat.kendall = zeros(size(U,1),size(U,2));
for i = 1:size(U,1)
    for j = 1:size(U,2)
        subplot(size(U,1),size(U,2),cont);
        reU = reshape(U(i,j,:),1,size(U,3));
        plot(rundata.tperiod,reU(1:end-1)); hold all;
%         if (i == 1) 
%             scatter(reU(1:24:end),prediction.SM(1:24:end), '.', 'k');
%             predVAR = prediction.SM(1:24:end);
%         else if (i == 2)
%                 scatter(reU(1:24:end),prediction.ISU(1:24:end), '.', 'k');
%                 predVAR = prediction.ISU(1:24:end);
%             else
%                 scatter(reU(1:24:end),prediction.ISO(i-2,1:24:end), '.', 'k');
%                 predVAR = prediction.ISO(i-2,1:24:end);
%             end
%         end
%         set(gca, 'Xlim', [min([min(reU(1:24:end)), -0.00001]) max([max(reU(1:24:end)),0.00001])], 'Ylim', [min(predVAR) max(predVAR)]);
%         CCmat.pearson(i,j) = corr(prediction.Q',reU(1:end-1)', 'Type', 'Pearson');
%         CCmat.spearman(i,j) = corr(prediction.Q',reU(1:end-1)', 'Type', 'Spearman');
%         CCmat.kendall(i,j) = corr(prediction.Q',reU(1:end-1)', 'Type', 'Kendall');
        title(['x_', num2str(i), '/a_', num2str(j)], 'FontSize', 10);
        set(gca, 'Xlim', [rundata.tperiod(1) rundata.tperiod(end)], 'Xtick', []);
        
%         set(gca, 'Xlim', [datenum(2006,6,7,0,0,0) datenum(2006,6,26,0,0,0)], 'Xtick', []);
%         xlabel('Time Period (hourly)', 'FontSize', 10);
        cont = cont + 1;
    end
end

obs_peak_time = rundata.tperiod(rundata.obsQ == max(rundata.obsQ));
sim_peak_time = rundata.tperiod(prediction.Q == max(prediction.Q));
diff_ptime_hrs(cont) = (sim_peak_time - obs_peak_time(1))*24;
analytical_lag_peaktime_hrs(cont) = (sim_peak_time - rainPeakdate)*24;
lag_peaktime_hrs(cont) = (sim_peak_time - rundata.tperiod(rundata.PP == max(rundata.PP)))*24;
cont = cont + 1;
end

% obsq_h = plot(rundata.tperiod,rundata.obsQ, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);

% lx = log(0.005:0.005:0.2);
% ly = log(lag_peaktime_hrs);
% 
% pfit = polyfit(ly,lx,1);
% figure(1)
% plot(lx,pfit(1)*lx+pfit(2))
% 
% %% Animation Sensitivity in Time
% for t = 1:size(V,3)
%     Vmat = reshape(V(:,:,t),size(V,1),size(V,2));
% end

%% Condition Number Analysis
% Ls = [1 6 12 24 24*3 24*5 24*7 24*7*2 24*30 24*30*3 24*30*6 24*30*12 24*30*12*2 24*30*12*3 24*30*12*4 24*30*12*5 24*30*12*6 24*30*12*7];
% CN_surf = zeros(length(Ls),12);
% for i = 1:18
%     for j = 1:12
%         CN_surf(i,j) = cond(H(1:Ls(i),1:j)'*H(1:Ls(i),1:j));
%     end
% end
% 
% %Moving Window (Period) Analysis
% cont = 1;
% for i = 1:90*24:61367  
%     CN_series(cont) = cond(H(i:i+8760,:)'*H(i:i+8760,:));
%     cont = cont + 1;
% end