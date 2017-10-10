clear; clc;

%Input Data: Tar River basin
dd = dlmread('/Users/humbertovergara/Documents/Trabajo/Parameter_and_QPEs_Error_Characterization_for_DA/TS_Data/MPE_tarboro_basin_avg_time_series.csv', ',', 1,1);
rundata.PP = dd(:,2); %load('~/Documents/MS_Work/Thesis_Work/HyMOD/Rain_MPE_200206-200912.txt'); 
rundata.PP_period = datenum('01-Jan-2002 10:00:00'):1/24:datenum('31-Dec-2009 23:00:00');
rundata.PET = dlmread('~/Documents/MS_Work/Thesis_Work/HyMOD/PET_mmhr_Tarboro_2002-2009.txt', ' ', 1,0);  
rundata.PET_period = datenum('01-Jan-2002 00:00:00'):1/24:datenum('31-Dec-2009 23:00:00');
rundata.obsQ  = dlmread('~/Documents/MS_Work/Thesis_Work/HyMOD/02083500_ObsQ_2002-2010_hrly.csv', ',', 1,1); %5709.28 km2 Tarboro drainage area  
rundata.obsQ_period = datenum('01-Jan-2002 00:00:00'):1/24:datenum('30-Jun-2010 21:00:00');
barea = 5709;  %Basin Area in km^2
runset.barea = barea;

%Simulation Data and Settings
runset.tstep = 1; %Time Step in hours
runset.tperiod = datenum(2002,2,1,1,0,0):1/24:datenum(2009,12,31,23,0,0);
rundata.runset = runset;

%Subset data to period of interest
[new_rundata] = subset(rundata,0);
rundata = new_rundata;

%CREST Parameters
refset = load('~/Documents/Trabajo/NESSF/NESSF2013/Sat_Ensemble_DA_Demo/Tar2006_hrly_DREAM_bestPar.mat');
Parameters.PKE = refset.bestPar(1); %Multiplier to convert PET
Parameters.PIM = refset.bestPar(2); %The impervious area ratio. This ratio is stored as a percentage between 0 and 100.
Parameters.PWM = refset.bestPar(3); %The maximum soil water capacity (depth integrated pore space) of the soil layer. (mm)
Parameters.PFC = refset.bestPar(4); %The soil saturated hydraulic conductivity (Ksat ,mm/hr).
Parameters.LEAKO = refset.bestPar(5); %The overland reservoir discharge multiplier.
Parameters.LEAKI = refset.bestPar(6); %The interflow reservoir discharge multiplier.
Parameters.PB = refset.bestPar(7); %The exponent of the variable infiltration curve.

%NASH Routing Parameters
Parameters.Nqo = 3;

%Approximation Parameters
alpha.alpha1 = 10^20;
beta.beta1 = 10^3;
alpha.alpha2 = 10^20;
beta.beta2 = 10^3;
alpha.alpha3 = 10^20;
beta.beta3 = 10^3;
alpha.alpha5 = 10^20;
beta.beta5 = 10^3;
alpha.alpha6 = 10^20;
beta.beta6 = 10^3;
alpha.alpha7 = 10^20; %sensitive when 10^0
beta.beta7 = 10^3;
alpha.alpha8 = 10^20; %sensitive
beta.beta8 = 10^3;
alpha.alpha9 = 10^20; %sensitive
beta.beta9 = 10^3;
alpha.alpha10 = 10^20;
beta.beta10 = 10^3;
alpha.alpha11 = 10^20;
beta.beta11 = 10^3;
alpha.alpha12 = 10^20;
beta.beta12 = 10^3;

Parameters.alpha = alpha;
Parameters.beta = beta;

%CREST States
States_0.IWU = 0;
States_0.ISO = [0 0 0];
States_0.ISU = 0;
rundata.States = States_0;

%Run CCREST
prediction = crest(rundata,refset.bestPar);

%Run Forward Sensitivity Method (FSM)
% [correction, V, U, A, H, Dh] = fsmccrest(rundata,Parameters,prediction,alpha,beta);
correction = [0 0 0 0 0]';

%Make Correction and run CCREST
%New CREST States
States_0.IWU = States_0.IWU+(correction(1)*100/Parameters.PWM);
States_0.ISU = States_0.ISU+correction(2);
States_0.ISO = States_0.ISO+correction(3:5)';
rundata.States = States_0;
newprediction = ccrest(rundata,Parameters);

% %% Plot Results
% %Plot and Stats
% subplot(1,2,1);
% maxv = max([prediction.I(:)' newprediction.M2(:)']);
% minv = min([prediction.I(:)' newprediction.M2(:)']);
% plot(minv(1):maxv(1),minv(1):maxv(1), 'LineStyle', '-', 'Color', 'k'); hold all;
% scatter(prediction.I,newprediction.M2, '.', 'k');
% set(gca, 'FontSize', 12, 'Xtick', 0:2:10, 'Xlim', [minv(1) maxv(1)], 'Ylim', [minv(1) maxv(1)], 'PlotBoxAspectRatio', [1 1 1]);
% ylabel('CCREST I_t', 'FontSize', 16); xlabel('CREST I_t', 'FontSize', 16);
% 
% subplot(1,2,2);
% maxv = max([prediction.ER(:)' newprediction.M3(:)']);
% minv = min([prediction.ER(:)' newprediction.M3(:)']);
% plot(minv(1):maxv(1),minv(1):maxv(1), 'LineStyle', '-', 'Color', 'k'); hold all;
% scatter(prediction.ER,newprediction.M3, '.', 'k');
% set(gca, 'FontSize', 12, 'Xlim', [minv(1) maxv(1)], 'Ylim', [minv(1) maxv(1)], 'PlotBoxAspectRatio', [1 1 1]);
% ylabel('CCREST ER_t', 'FontSize', 16); xlabel('CREST ER_t', 'FontSize', 16);
% 
% % saveas(gcf, 'CCREST_vs_CREST_ScatterPlots_It_ERt_Tar_2002-2009.fig');
% % saveas(gcf, 'CCREST_vs_CREST_ScatterPlots_It_ERt_Tar_2002-2009.png');
% % close;
% 
% figure;
% subplot(1,2,1);
% maxv = max([prediction.ERI(:)' newprediction.M4(:)']);
% minv = min([prediction.ERI(:)' newprediction.M4(:)']);
% plot(minv(1):maxv(1),minv(1):maxv(1), 'LineStyle', '-', 'Color', 'k'); hold all;
% scatter(prediction.ERI,newprediction.M4, '.', 'k');
% set(gca, 'FontSize', 12, 'Xlim', [minv(1) maxv(1)], 'Ylim', [minv(1) maxv(1)], 'PlotBoxAspectRatio', [1 1 1]);
% ylabel('CCREST ER_I_,_t', 'FontSize', 16); xlabel('CREST ER_I_,_t', 'FontSize', 16);
% 
% subplot(1,2,2);
% maxv = max([prediction.temX(:)' newprediction.M5(:)']);
% minv = min([prediction.temX(:)' newprediction.M5(:)']);
% plot(minv(1):maxv(1),minv(1):maxv(1), 'LineStyle', '-', 'Color', 'k'); hold all;
% scatter(prediction.temX,newprediction.M5, '.', 'k');
% set(gca, 'FontSize', 12, 'Xlim', [minv(1) maxv(1)], 'Ylim', [minv(1) maxv(1)], 'PlotBoxAspectRatio', [1 1 1]);
% ylabel('CCREST temX_t', 'FontSize', 16); xlabel('CREST temX_t', 'FontSize', 16);
% 
% % saveas(gcf, 'CCREST_vs_CREST_ScatterPlots_ERIt_temXt_Tar_2002-2009.fig');
% % saveas(gcf, 'CCREST_vs_CREST_ScatterPlots_ERIt_temXt_Tar_2002-2009.png');
% % close;
% 
% % subplot(1,2,1);
% % maxv = max([prediction.W0(:)' newprediction.SM(:)']);
% % minv = min([prediction.W0(:)' newprediction.SM(:)']);
% % plot(minv(1):maxv(1),minv(1):maxv(1), 'LineStyle', '-', 'Color', 'k'); hold all;
% % scatter(prediction.W0(1:end-1),newprediction.SM(2:end-1), '.', 'k');
% % set(gca, 'FontSize', 12, 'Xlim', [minv(1) maxv(1)], 'Ylim', [minv(1) maxv(1)], 'PlotBoxAspectRatio', [1 1 1]);
% % ylabel('CCREST SM_t', 'FontSize', 16); xlabel('CREST SM_t', 'FontSize', 16);
% % 
% % 
% % subplot(1,2,2);
% % maxv = max([prediction.ERO(:)' newprediction.M7(:)']);
% % minv = min([prediction.ERO(:)' newprediction.M7(:)']);
% % plot(minv(1):maxv(1),minv(1):maxv(1), 'LineStyle', '-', 'Color', 'k'); hold all;
% % scatter(prediction.ERO,newprediction.M7, '.', 'k');
% % set(gca, 'FontSize', 12, 'Xlim', [minv(1) maxv(1)], 'Ylim', [minv(1) maxv(1)], 'PlotBoxAspectRatio', [1 1 1]);
% % ylabel('CCREST ER_O_,_t', 'FontSize', 16); xlabel('CREST ER_O_,_t', 'FontSize', 16);
% % 
% % saveas(gcf, 'CCREST_vs_CREST_ScatterPlots_EROt_SMt_Tar_2002-2009.fig');
% % saveas(gcf, 'CCREST_vs_CREST_ScatterPlots_EROt_SMt_Tar_2002-2009.png');
% % close;
% 
% figure;
% st_idx = find(rundata.obsQ_period == runset.tperiod(1));
% en_idx = find(rundata.obsQ_period == runset.tperiod(end));
% rundata.obsQ = rundata.obsQ(st_idx:en_idx);
% [bias(1),rmse(1),nsce(1)] = nanhydrostat(rundata.obsQ',prediction.Q);
% [bias(2),rmse(2),nsce(2)] = nanhydrostat(rundata.obsQ',newprediction.Q);
% [bias(3),rmse(3),nsce(3)] = nanhydrostat(prediction.Q,newprediction.Q);
% 
% % area(rundata.obsQ); hold all;
% area(runset.tperiod,prediction.Q); hold all;
% plot(runset.tperiod,newprediction.Q, 'DisplayName', 'After FSM', 'Color', 'r');
% set(gca, 'FontSize', 12);
% xlabel('Time (Hours)', 'FontSize', 14); ylabel('Streamflow (m^3/s)', 'FontSize', 14);
% 
% % close; plot(newprediction.test);
% % x = newprediction.test;
% % x(x>0 & x<1) = NaN;
% % y = isnan(x);
% % length(find(y == 1))
% 
% %% Plot Sensitivity
% % figure;
% % cont = 1;
% % for i = 1:size(V,1)
% %     for j = 1:size(V,2)
% %         subplot(size(V,1),size(V,2),cont);
% %         plot(reshape(V(i,j,:),1,size(V,3)));
% %         title(['V_k x_', num2str(i), '/x_', num2str(j), '_(_0_)'], 'FontSize', 14);
% %         %xlabel('Time (hours)', 'FontSize', 10);
% %         cont = cont + 1;
% %     end
% % end