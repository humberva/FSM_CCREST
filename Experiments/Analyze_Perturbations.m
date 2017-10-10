%% Brief description of script
%This script compares the evolution of the ensemble of model solutions
%given a sample of initial conditions to the evolution of the perturbations
%through FSM.

%% Load data
clear;
clc;

% load('Ensemble_Simulation_Initial_Conditions_Perturbation_Analysis_n1000ens_MaxVariance.mat');
% load('Ensemble_Simulation_Initial_Conditions_Perturbation_Analysis_n1000ens_BalancedVariance.mat');
load('Ensemble_Simulation_Initial_Conditions_Perturbation_Analysis_n1000ens_BalancedVariance_No_alpha_beta_correction.mat');
% load('Ensemble_Simulation_Initial_Conditions_Perturbation_Analysis_n1000ens_MinVariance.mat');

load('../Selected_Ensembles_of_States_n1000.mat')

%% Plot histograms of model states samples
% figure;
% set(gcf, 'Position', [440,-152,1507,950]);
% for x = 1:5
%     subplot(2,3,x);
%     [count, bins] = hist(maxVar_States(x,:),20);
%     bar(bins,count, 'EdgeColor', 'k', 'FaceColor', [0.7 0.7 0.7], 'BarWidth', 1);
%     ylabel('Number of Occurrences', 'FontSize', 14);
%     xlabel(['$x_',num2str(x),'$'], 'interpreter', 'latex', 'FontSize', 18);
%     set(gca, 'Xlim', [min(bins) max(bins)], 'FontSize', 12);
% end


%% Pre-allocation of variables
% Variables to hold traces of covariance matrices
totVar_pertVk = nan(1,size(pertXv,3));
totVar_ens_x = nan(1,size(pertXv,3));
totVar_ens_hat_x = nan(1,size(pertXv,3));

%Variables for model solution approximation through FSM
mean_xi = nan(5,size(pertXv,3));

x1_hat_fs = nan(1000,size(pertXv,3));
x2_hat_fs = nan(1000,size(pertXv,3));
x3_hat_fs = nan(1000,size(pertXv,3));
x4_hat_fs = nan(1000,size(pertXv,3));
x5_hat_fs = nan(1000,size(pertXv,3));

for t = 1:size(pertXv,3)-1
    totVar_pertVk(t) = trace(nancov(reshape(pertXv(:,:,t),1000,5)));
    totVar_ens_x(t) = trace(nancov([soil_moisture(:,t), slow_reservoir(:,t),reshape(quick_reservoirs(:,t,:),1000,3)]));
    mean_xi(:,t) = nanmean([soil_moisture(:,t), slow_reservoir(:,t),reshape(quick_reservoirs(:,t,:),1000,3)])';
    
    %Model solution approximation through FSM
    %Option 1: Mean of the ensemble of model solutions
    x1_hat_fs(:,t) = mean_xi(1,t)+reshape(pertXv(1,:,t),1000,1);
    x2_hat_fs(:,t) = mean_xi(2,t)+reshape(pertXv(2,:,t),1000,1);
    x3_hat_fs(:,t) = mean_xi(3,t)+reshape(pertXv(3,:,t),1000,1);
    x4_hat_fs(:,t) = mean_xi(4,t)+reshape(pertXv(4,:,t),1000,1);
    x5_hat_fs(:,t) = mean_xi(5,t)+reshape(pertXv(5,:,t),1000,1);
    
%     %Option 2: Mean of the ensemble of initial states to run model forward
%     x1_hat_fs(:,t) = mean_prediction.SM(t)+reshape(pertXv(1,:,t),1000,1);
%     x2_hat_fs(:,t) = mean_prediction.ISU(t)+reshape(pertXv(2,:,t),1000,1);
%     x3_hat_fs(:,t) = mean_prediction.ISO(1,t)+reshape(pertXv(3,:,t),1000,1);
%     x4_hat_fs(:,t) = mean_prediction.ISO(2,t)+reshape(pertXv(4,:,t),1000,1);
%     x5_hat_fs(:,t) = mean_prediction.ISO(3,t)+reshape(pertXv(5,:,t),1000,1);
    
    totVar_ens_hat_x(t) = trace(nancov([x1_hat_fs(:,t), x2_hat_fs(:,t), x3_hat_fs(:,t), x4_hat_fs(:,t), x5_hat_fs(:,t)]));
end

%Variance of each ensemble model component solution
var_ens_x(1,:) = nanvar(soil_moisture);
var_ens_x(2,:) = nanvar(slow_reservoir);
var_ens_x(3,:) = nanvar(reshape(quick_reservoirs(:,:,1),1000,480));
var_ens_x(4,:) = nanvar(reshape(quick_reservoirs(:,:,2),1000,480));
var_ens_x(5,:) = nanvar(reshape(quick_reservoirs(:,:,3),1000,480));

%Variance of each ensemble FS approximation solution
var_x_hat_fs(1,:) = nanvar(x1_hat_fs);
var_x_hat_fs(2,:) = nanvar(x2_hat_fs);
var_x_hat_fs(3,:) = nanvar(x3_hat_fs);
var_x_hat_fs(4,:) = nanvar(x4_hat_fs);
var_x_hat_fs(5,:) = nanvar(x5_hat_fs);

%Compute quantiled time series for easy plot overlaying
quant_pdf = 99.99;
quant_pdf_low = 100 - quant_pdf;

%Ensemble model solution
quant_up_ens_x(1,:) = prctile(soil_moisture,quant_pdf);
quant_up_ens_x(2,:) = prctile(slow_reservoir,quant_pdf);
quant_up_ens_x(3,:) = prctile(reshape(quick_reservoirs(:,:,1),1000,480),quant_pdf);
quant_up_ens_x(4,:) = prctile(reshape(quick_reservoirs(:,:,2),1000,480),quant_pdf);
quant_up_ens_x(5,:) = prctile(reshape(quick_reservoirs(:,:,3),1000,480),quant_pdf);

quant_lo_ens_x(1,:) = prctile(soil_moisture,quant_pdf_low);
quant_lo_ens_x(2,:) = prctile(slow_reservoir,quant_pdf_low);
quant_lo_ens_x(3,:) = prctile(reshape(quick_reservoirs(:,:,1),1000,480),quant_pdf_low);
quant_lo_ens_x(4,:) = prctile(reshape(quick_reservoirs(:,:,2),1000,480),quant_pdf_low);
quant_lo_ens_x(5,:) = prctile(reshape(quick_reservoirs(:,:,3),1000,480),quant_pdf_low);

%Ensemble FS approximation
quant_up_ens_hat_x(1,:) = prctile(x1_hat_fs,quant_pdf);
quant_up_ens_hat_x(2,:) = prctile(x2_hat_fs,quant_pdf);
quant_up_ens_hat_x(3,:) = prctile(x3_hat_fs,quant_pdf);
quant_up_ens_hat_x(4,:) = prctile(x4_hat_fs,quant_pdf);
quant_up_ens_hat_x(5,:) = prctile(x5_hat_fs,quant_pdf);

quant_lo_ens_hat_x(1,:) = prctile(x1_hat_fs,quant_pdf_low);
quant_lo_ens_hat_x(2,:) = prctile(x2_hat_fs,quant_pdf_low);
quant_lo_ens_hat_x(3,:) = prctile(x3_hat_fs,quant_pdf_low);
quant_lo_ens_hat_x(4,:) = prctile(x4_hat_fs,quant_pdf_low);
quant_lo_ens_hat_x(5,:) = prctile(x5_hat_fs,quant_pdf_low);

%% Comparison between mean of ensemble solution and mean initial conditions
figure;
for x = 1:5
    subplot(2,3,x); %+5)
    plot(mean_xi(x,:)', 'Color', 'b');
    hold all;
    ylabstr = ['$x_', num2str(x),'$'];
    ylabel(ylabstr, 'interpreter', 'latex', 'fontsize',20)
    set(gca, 'Xlim', [0 480], 'FontSize', 12);
    xlabel('Time Steps', 'FontSize', 14);
end

subplot(2,3,1);
plot(mean_prediction.SM, 'Color', 'r');
hleg = legend({'$\bar x_k$', '$x_k = M(\bar x_0)$'}, 'Interpreter', 'latex');
set(hleg, 'fontsize',16);

subplot(2,3,2);
plot(mean_prediction.ISU, 'Color', 'r');

for x = 1:3
    subplot(2,3,x+2)
    plot(mean_prediction.ISO(x,:), 'Color', 'r');
    hold all;
    ylabstr = ['$x_', num2str(x+2),'$'];
    ylabel(ylabstr, 'interpreter', 'latex', 'fontsize',20)
    set(gca, 'Xlim', [0 480], 'FontSize', 12);
    xlabel('Time Steps', 'FontSize', 14);
end

%% Ensemble Simulations
figure;
set(gcf, 'Position', [156,-152,2336,950]);
subplot(2,5,1)
plot(soil_moisture', 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('$x_1$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Soil Moisture', 'FontSize', 16);

subplot(2,5,2)
plot(slow_reservoir', 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('$x_2$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Slow Reservoir', 'FontSize', 16);

subplot(2,5,3)
plot(quick_reservoirs(:,:,1)', 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('$x_3$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Quick Reservoir 1', 'FontSize', 16);

subplot(2,5,4)
plot(quick_reservoirs(:,:,2)', 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('$x_4$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Quick Reservoir 2', 'FontSize', 16);

subplot(2,5,5)
plot(quick_reservoirs(:,:,3)', 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('$x_5$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Quick Reservoir 3', 'FontSize', 16);

for x = 1:5
    subplot(2,5,x+5)
    eval(['xk_hat_fs = x', num2str(x), '_hat_fs;']);
    plot(xk_hat_fs', 'Color', [0.7 0.7 0.7]);
    hold all;
    plot(0:480,zeros(1,481), 'Color', 'k');
    ylabstr = ['$\hat{x}_', num2str(x),'$'];
    ylabel(ylabstr, 'interpreter', 'latex', 'fontsize',20)
    set(gca, 'Xlim', [0 480], 'FontSize', 12);
    xlabel('Time Steps', 'FontSize', 14);
end

%% Extra figure to overlay plots using quantiles of ensembles
figure;
for x = 1:5
    subplot(2,3,x);
    h_xq = area(quant_up_ens_x(x,:), 'EdgeColor', [0.7 0.7 0.7], 'FaceColor', [0.7 0.7 0.7]);
    hold all;
    h_xq2 = area(quant_lo_ens_x(x,:), 'EdgeColor', 'w', 'FaceColor', 'w');
    h_hxq = plot(quant_up_ens_hat_x(x,:), 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
    h_hxq2 = plot(quant_lo_ens_hat_x(x,:), 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
    legend([h_xq, h_hxq], {['$x_', num2str(x), '$'], ['$\hat{x}_', num2str(x), '$']}, 'Interpreter', 'latex', 'FontSize', 12);
    set(gca, 'Xlim', [0 480], 'FontSize', 12);
    ylabel(['$x_', num2str(x), '$'], 'interpreter', 'latex', 'fontsize',20);
    xlabel('Time Steps', 'FontSize', 14);
    title([num2str(quant_pdf_low), ' - ', num2str(quant_pdf), ' Quantiles'], 'FontSize', 14);
end

%% Non-linearity Analysis
figure;
set(gcf, 'Position', [440,-152,1507,950]);

subplot(3,2,1);
plot(var_ens_x(1,:), 'LineStyle', '-', 'Color', 'k');
hold all;
plot(var_x_hat_fs(1,:),'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
text(380, (min(var_ens_x(1,:)) + 0.6*(max(var_ens_x(1,:))-min(var_ens_x(1,:)))), ['Linear Corr = ', num2str(round(corr(var_ens_x(1,:)',var_x_hat_fs(1,1:480)')*100)/100, '%.2f')], 'FontSize', 14);
hleg = legend({'$x_1$', '$\hat{x}_1$'}, 'Interpreter', 'latex');
set(hleg, 'fontsize',16);
set(gca, 'FontSize', 12);
ylabel('Total Variance', 'fontsize',18)
xlabel('Time steps', 'FontSize', 14);
title('Soil Moisture', 'FontSize', 16);

subplot(3,2,2);
plot(var_ens_x(2,:), 'LineStyle', '-', 'Color', 'k');
hold all;
plot(var_x_hat_fs(2,:),'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
text(380, (min(var_ens_x(2,:)) + 0.6*(max(var_ens_x(2,:))-min(var_ens_x(2,:)))), ['Linear Corr = ', num2str(round(corr(var_ens_x(2,:)',var_x_hat_fs(2,1:480)')*100)/100, '%.2f')], 'FontSize', 14);
hleg = legend({'$x_2$', '$\hat{x}_2$'}, 'Interpreter', 'latex');
set(hleg, 'fontsize',16);
set(gca, 'FontSize', 12);
ylabel('Total Variance', 'fontsize',18)
xlabel('Time steps', 'FontSize', 14);
title('Slow Reservoir', 'FontSize', 16);

subplot(3,2,3);
plot(var_ens_x(3,:), 'LineStyle', '-', 'Color', 'k');
hold all;
plot(var_x_hat_fs(3,:),'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
text(380, (min(var_ens_x(3,:)) + 0.6*(max(var_ens_x(3,:))-min(var_ens_x(3,:)))), ['Linear Corr = ', num2str(round(corr(var_ens_x(3,:)',var_x_hat_fs(3,1:480)')*100)/100, '%.2f')], 'FontSize', 14);
hleg = legend({'$x_3$', '$\hat{x}_3$'}, 'Interpreter', 'latex');
set(hleg, 'fontsize',16);
set(gca, 'FontSize', 12);
ylabel('Total Variance', 'fontsize',18);
xlabel('Time steps', 'FontSize', 14);
title('Quick Reservoir 1', 'FontSize', 16);

subplot(3,2,4);
plot(var_ens_x(4,:), 'LineStyle', '-', 'Color', 'k');
hold all;
plot(var_x_hat_fs(4,:),'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
text(380, (min(var_ens_x(4,:)) + 0.6*(max(var_ens_x(4,:))-min(var_ens_x(4,:)))), ['Linear Corr = ', num2str(round(corr(var_ens_x(4,:)',var_x_hat_fs(4,1:480)')*100)/100, '%.2f')], 'FontSize', 14);
hleg = legend({'$x_4$', '$\hat{x}_4$'}, 'Interpreter', 'latex');
set(hleg, 'fontsize',16);
set(gca, 'FontSize', 12);
ylabel('Total Variance', 'fontsize',18);
xlabel('Time steps', 'FontSize', 14);
title('Quick Reservoir 2', 'FontSize', 16);

subplot(3,2,5);
plot(var_ens_x(5,:), 'LineStyle', '-', 'Color', 'k');
hold all;
plot(var_x_hat_fs(5,:),'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
text(380, (min(var_ens_x(5,:)) + 0.6*(max(var_ens_x(5,:))-min(var_ens_x(5,:)))), ['Linear Corr = ', num2str(round(corr(var_ens_x(5,:)',var_x_hat_fs(5,1:480)')*100)/100, '%.2f')], 'FontSize', 14);
hleg = legend({'$x_5$', '$\hat{x}_5$'}, 'Interpreter', 'latex');
set(hleg, 'fontsize',16);
set(gca, 'FontSize', 12);
ylabel('Total Variance', 'fontsize',18);
xlabel('Time steps', 'FontSize', 14);
title('Quick Reservoir 3', 'FontSize', 16);

subplot(3,2,6);
plot(totVar_ens_x, 'LineStyle', '-', 'Color', 'k');
hold all;
plot(totVar_ens_hat_x,'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
text(380, (min(totVar_ens_x) + 0.6*(max(totVar_ens_x)-min(totVar_ens_x))), ['Linear Corr = ', num2str(round(corr(totVar_ens_x(1:480)',totVar_ens_hat_x(1:480)')*100)/100, '%.2f')], 'FontSize', 14);
hleg = legend({'$x$', '$\hat{x}$'}, 'Interpreter', 'latex');
set(hleg, 'fontsize',16);
set(gca, 'FontSize', 12);
ylabel('Total Variance', 'fontsize',18);
xlabel('Time steps', 'FontSize', 14);
title('Trace of Covariances', 'FontSize', 16);



% %% 4-D View of variance dynamics
% for t = 1 %:10:480
%     subplot(1,2,1);
%     scatter3(soil_moisture(:,t),slow_reservoir(:,t),quick_reservoirs(:,t,3),10,'k', 'filled');
%     xlabel('Soil Moisture'); ylabel('Slow Reservoir'); zlabel('Quick Reservoir 1');
%     set(gca, 'Xlim', [60 210], 'Ylim', [15 40], 'Zlim', [0 75]);
%     
%     subplot(1,2,2);
%     scatter3(reshape(pertXv(1,:,t),1000,1)',reshape(pertXv(2,:,t),1000,1)',reshape(pertXv(3,:,t),1000,1)',(reshape(pertXv(4,:,t),1000,1)'+10).^2,reshape(pertXv(5,:,t),1000,1)', 'filled');
%     xlabel('Soil Moisture'); ylabel('Slow Reservoir'); zlabel('Quick Reservoir 1');
%     set(gca, 'Xlim', [-55 0], 'Ylim', [-12 8], 'Zlim', [-15 10]);
%     title(['t = ', num2str(t)]);
% %     colorbar; caxis([-5 4])
%     getframe();
% end