clear;
clc;

% load('SERVER_Results/Ensemble_Simulation_Initial_Conditions.mat');
load('Random_Initial_Conditions_100members_TS_Alberto_02083500.mat');
load('SERVER_Results/Ensemble_Simulation_Initial_Conditions_Perturbation_Analysis.mat');

%% Sample of states
figure;
subplot(2,3,1);
[counts, bin_centers] = hist(random_states.SM);
bar(bin_centers,(counts./sum(counts)).*100);
xlabel('Soil Saturation (mm)', 'FontSize', 14);
ylabel('Relative count (%)', 'FontSize', 14);
set(gca, 'Xlim', [min(bin_centers)-round(mean(abs(diff(bin_centers)))) max(bin_centers)+round(mean(abs(diff(bin_centers))))], 'FontSize', 12);

subplot(2,3,2);
[counts, bin_centers] = hist(random_states.ISU);
bar(bin_centers,(counts./sum(counts)).*100);
xlabel('Slow reservoir content (mm)', 'FontSize', 14);
ylabel('Relative count (%)', 'FontSize', 14);
set(gca, 'Xlim', [min(bin_centers)-round(mean(abs(diff(bin_centers)))) max(bin_centers)+round(mean(abs(diff(bin_centers))))], 'FontSize', 12);

subplot(2,3,3);
[counts, bin_centers] = hist(random_states.ISO(1,:));
bar(bin_centers,(counts./sum(counts)).*100);
xlabel('Quick reservoir 1 content (mm)', 'FontSize', 14);
ylabel('Relative count (%)', 'FontSize', 14);
set(gca, 'Xlim', [min(bin_centers)-round(mean(abs(diff(bin_centers)))) max(bin_centers)+round(mean(abs(diff(bin_centers))))], 'FontSize', 12);

subplot(2,3,4);
[counts, bin_centers] = hist(random_states.ISO(2,:));
bar(bin_centers,(counts./sum(counts)).*100);
xlabel('Quick reservoir 2 content (mm)', 'FontSize', 14);
ylabel('Relative count (%)', 'FontSize', 14);
set(gca, 'Xlim', [min(bin_centers)-round(mean(abs(diff(bin_centers)))) max(bin_centers)+round(mean(abs(diff(bin_centers))))], 'FontSize', 12);

subplot(2,3,5);
[counts, bin_centers] = hist(random_states.ISO(3,:));
bar(bin_centers,(counts./sum(counts)).*100);
xlabel('Quick reservoir 3 content (mm)', 'FontSize', 14);
ylabel('Relative count (%)', 'FontSize', 14);
set(gca, 'Xlim', [min(bin_centers)-round(mean(abs(diff(bin_centers)))) max(bin_centers)+round(mean(abs(diff(bin_centers))))], 'FontSize', 12);

%% Model Solution
figure;
subplot(2,3,1);
plot(streamflow', 'Color', [0.7 0.7 0.7]);
xlabel('Simulation Time (hours)', 'FontSize', 14);
ylabel('Streamflow (m^3/s)', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 12);

subplot(2,3,2);
plot(soil_moisture', 'Color', [0.7 0.7 0.7]);
xlabel('Simulation Time (hours)', 'FontSize', 14);
ylabel('Soil Saturation (mm/hr)', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 12);

subplot(2,3,3);
plot(slow_reservoir', 'Color', [0.7 0.7 0.7]);
xlabel('Simulation Time (hours)', 'FontSize', 14);
ylabel('Slow reservoir flow (mm/hr)', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 12);

subplot(2,3,4);
plot(quick_reservoirs(:,:,1)', 'Color', [0.7 0.7 0.7]);
xlabel('Simulation Time (hours)', 'FontSize', 14);
ylabel('Quick reservoir flow 1 (mm/hr)', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 12);

subplot(2,3,5);
plot(quick_reservoirs(:,:,2)', 'Color', [0.7 0.7 0.7]);
xlabel('Simulation Time (hours)', 'FontSize', 14);
ylabel('Quick reservoir flow 2 (mm/hr)', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 12);

subplot(2,3,6);
plot(quick_reservoirs(:,:,3)', 'Color', [0.7 0.7 0.7]);
xlabel('Simulation Time (hours)', 'FontSize', 14);
ylabel('Quick reservoir flow 3 (mm/hr)', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 12);

%% FSM functions

cont = 0;
n = size(ensU,1);
p = size(ensU,2);
figure;
for i = 1:n
   for j = 1:p
       cont = cont + 1;
       subplot(n,p,cont);
       plot(reshape(ensU(i,j,:,:),481,100), 'Color', [0.7 0.7 0.7]);
       set(gca, 'XTick', [], 'YTick', []);
   end
end

figure;
cont = 0;
for i = 1:n
   for j = 1:n
       cont = cont + 1;
       subplot(n,n,cont);
       plot(reshape(ensV(i,j,:,:),481,100), 'Color', [0.7 0.7 0.7]);
       set(gca, 'XTick', [], 'YTick', []);
   end
end

%% Covariance analysis

for t = 1:480
    Cov_Q(t) = nanvar(streamflow(:,t));
    Cov_SM(t) = nanvar(soil_moisture(:,t));
    
    for i = 1:n
       for j = 1:n
           covTEMP = nancov(streamflow(:,t),reshape(ensV(i,j,t,:),1,100));
           Cov_Q_V(i,j,t) = covTEMP(2);
           
           covTEMP = nancov(soil_moisture(:,t),reshape(ensV(i,j,t,:),1,100));
           Cov_SM_V(i,j,t) = covTEMP(2);
           
           covTEMP = nancov(slow_reservoir(:,t),reshape(ensV(i,j,t,:),1,100));
           Cov_SQ_V(i,j,t) = covTEMP(2);
           
           covTEMP = nancov(reshape(quick_reservoirs(:,t,1),1,100),reshape(ensV(i,j,t,:),1,100));
           Cov_QQ1_V(i,j,t) = covTEMP(2);
           
           covTEMP = nancov(reshape(quick_reservoirs(:,t,2),1,100),reshape(ensV(i,j,t,:),1,100));
           Cov_QQ2_V(i,j,t) = covTEMP(2);
           
           covTEMP = nancov(reshape(quick_reservoirs(:,t,3),1,100),reshape(ensV(i,j,t,:),1,100));
           Cov_QQ3_V(i,j,t) = covTEMP(2);
       end
    end
    
%     for i = 1:n
%        for j = 1:p
%            Cov_U(i,j,t) = nanvar(reshape(ensU(i,j,t,:),1,100));
%        end
%     end
end

%% Perturbation analysis
figure;
cont = 0;
Covar_PertXv = nan(5,5,481);
for i = 1:n
    for t = 1:481
       cont = cont + 1;
%        subplot(3,2,cont);
%        plot(reshape(pertXv(i,:,:),481,100), 'Color', [0.7 0.7 0.7]);
       Covar_PertXv(:,:,t) = pertXv(i,:,t);
%        set(gca, 'XTick', [], 'YTick', []);
    end
end