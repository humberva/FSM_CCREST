%% Brief description of script
%This script compares the evolution of the ensemble of model solutions
%given a sample of parameter estimates to the evolution of the perturbations
%through FSM.

%% Load data
clear;
% clc;
Par = 'LEAKI';

fid = fopen([Par, '_Sensitivity_Tests_Variance_Ratio_Summary_Table.csv'], 'a+');
fid2 = fopen([Par, '_Sensitivity_Tests_Variance_Difference_Summary_Table.csv'], 'a+');

% load('Ensemble_Simulation_Parameters_Perturbation_Analysis_nsize(pertXu,2)ens.mat');
% benchmark = load('Nominal_Sensitivity_Parameters/Ensemble_Simulation_Nominal_Parameters_Perturbation_Analysis_baseline.mat');
% load('Nominal_Sensitivity_Parameters/Ensemble_Simulation_Nominal_Parameters_Perturbation_Analysis_PIM_40.001_0.001_40.01.mat');
% load('Nominal_Sensitivity_Parameters/ccrest_Ensemble_Simulation_Nominal_Parameters_Perturbation_Analysis_PIM_40.001_0.001_40.01.mat');
% load('Nominal_Sensitivity_Parameters/Ensemble_Simulation_Nominal_Parameters_Perturbation_Analysis_PIM_10_20_30_40.mat');
% load('Nominal_Sensitivity_Parameters/ccrest_Ensemble_Simulation_Nominal_Parameters_Perturbation_Analysis_PIM_2.2_0.001_2.21.mat');
load('Nominal_Sensitivity_Parameters/Ensemble_Simulation_Nominal_Parameters_Perturbation_Analysis_LEAKI_0.01_to_0.03_in_0.001s.mat');
%% Pre-allocation of variables
% Variables to hold traces of covariance matrices
totVar_pertUk = nan(1,size(pertXu,3));
totVar_ens_x = nan(1,size(pertXu,3));
totVar_ens_hat_x = nan(1,size(pertXu,3));

%Variables for model solution approximation through FSM
mean_xi = nan(5,size(pertXu,3));

x1_hat_fs = nan(size(pertXu,2),size(pertXu,3));
x2_hat_fs = nan(size(pertXu,2),size(pertXu,3));
x3_hat_fs = nan(size(pertXu,2),size(pertXu,3));
x4_hat_fs = nan(size(pertXu,2),size(pertXu,3));
x5_hat_fs = nan(size(pertXu,2),size(pertXu,3));

for t = 1:size(pertXu,3)-1
    totVar_pertUk(t) = trace(nancov(reshape(pertXu(:,:,t),size(pertXu,2),5)));
    totVar_ens_x(t) = trace(nancov([soil_moisture(:,t), slow_reservoir(:,t),reshape(quick_reservoirs(:,t,:),size(pertXu,2),3)]));
    mean_xi(:,t) = nanmean([soil_moisture(:,t), slow_reservoir(:,t),reshape(quick_reservoirs(:,t,:),size(pertXu,2),3)])';

    %Model solution approximation through FSM
    %Option 1: Mean of the ensemble of model solutions
    x1_hat_fs(:,t) = mean_xi(1,t)+reshape(pertXu(1,:,t),size(pertXu,2),1);
    x2_hat_fs(:,t) = mean_xi(2,t)+reshape(pertXu(2,:,t),size(pertXu,2),1);
    x3_hat_fs(:,t) = mean_xi(3,t)+reshape(pertXu(3,:,t),size(pertXu,2),1);
    x4_hat_fs(:,t) = mean_xi(4,t)+reshape(pertXu(4,:,t),size(pertXu,2),1);
    x5_hat_fs(:,t) = mean_xi(5,t)+reshape(pertXu(5,:,t),size(pertXu,2),1);
    
%     %Option 2: Mean of the ensemble of initial states to run model forward
%     x1_hat_fs(:,t) = benchmark.soil_moisture(t)+reshape(pertXu(1,:,t),size(pertXu,2),1);
%     x2_hat_fs(:,t) = benchmark.slow_reservoir(t)+reshape(pertXu(2,:,t),size(pertXu,2),1);
%     x3_hat_fs(:,t) = benchmark.quick_reservoirs(:,t,1)+reshape(pertXu(3,:,t),size(pertXu,2),1);
%     x4_hat_fs(:,t) = benchmark.quick_reservoirs(:,t,2)+reshape(pertXu(4,:,t),size(pertXu,2),1);
%     x5_hat_fs(:,t) = benchmark.quick_reservoirs(:,t,3)+reshape(pertXu(5,:,t),size(pertXu,2),1);
    
%     %Show SM as percentage (%)
%     soil_moisture(:,t) = soil_moisture(:,t)./Parameter_sets(2,:)';
%     x1_hat_fs(:,t) = x1_hat_fs(:,t)./Parameter_sets(2,:)';

    totVar_ens_hat_x(t) = trace(nancov([x1_hat_fs(:,t), x2_hat_fs(:,t), x3_hat_fs(:,t), x4_hat_fs(:,t), x5_hat_fs(:,t)]));
end

%Variance of each ensemble model component solution
var_ens_x(1,:) = nanvar(soil_moisture);
var_ens_x(2,:) = nanvar(slow_reservoir);
var_ens_x(3,:) = nanvar(reshape(quick_reservoirs(:,:,1),size(pertXu,2),480));
var_ens_x(4,:) = nanvar(reshape(quick_reservoirs(:,:,2),size(pertXu,2),480));
var_ens_x(5,:) = nanvar(reshape(quick_reservoirs(:,:,3),size(pertXu,2),480));

%Variance of each ensemble FS approximation solution
var_x_hat_fs(1,:) = nanvar(x1_hat_fs);
var_x_hat_fs(2,:) = nanvar(x2_hat_fs);
var_x_hat_fs(3,:) = nanvar(x3_hat_fs);
var_x_hat_fs(4,:) = nanvar(x4_hat_fs);
var_x_hat_fs(5,:) = nanvar(x5_hat_fs);


%% Ensemble Simulations
figure;
set(gcf, 'Position', [156,-152,2336,950]);
subplot(2,5,1)
plot(soil_moisture'); %, 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('$x_1$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Soil Moisture', 'FontSize', 16);

c_par_i = 5;
for leg_i = 1:size(Parameter_sets,2)
    legend_str{leg_i} = num2str(Parameter_sets(c_par_i,leg_i));
end

legend(legend_str);

subplot(2,5,2)
plot(slow_reservoir'); %, 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('$x_2$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Slow Reservoir', 'FontSize', 16);

subplot(2,5,3)
plot(quick_reservoirs(:,:,1)'); %, 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('$x_3$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Quick Reservoir 1', 'FontSize', 16);

subplot(2,5,4)
plot(quick_reservoirs(:,:,2)'); %, 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('$x_4$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Quick Reservoir 2', 'FontSize', 16);

subplot(2,5,5)
plot(quick_reservoirs(:,:,3)'); %, 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('$x_5$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Quick Reservoir 3', 'FontSize', 16);

for x = 1:5
    subplot(2,5,x+5)
    eval(['xk_hat_fs = x', num2str(x), '_hat_fs;']);
    plot(xk_hat_fs'); %, 'Color', [0.7 0.7 0.7]);
    hold all;
    plot(0:480,zeros(1,481), 'Color', 'k');
    ylabstr = ['$\hat{x}_', num2str(x),'$'];
    ylabel(ylabstr, 'interpreter', 'latex', 'fontsize',20)
    set(gca, 'Xlim', [0 480], 'FontSize', 12);
    xlabel('Time Steps', 'FontSize', 14);
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

fprintf(fid, '%.4f - %.4f - %.4f, %f, ', min(Parameter_sets(c_par_i,:)), mean(diff(Parameter_sets(c_par_i,:))),  ...
    max(Parameter_sets(c_par_i,:)), nanmedian(var_x_hat_fs(1,1:480)./var_ens_x(1,:)));
fprintf(fid, '%f, ', nanmedian(var_x_hat_fs(2,1:480)./var_ens_x(2,:)));
fprintf(fid, '%f, ', nanmedian(var_x_hat_fs(3,1:480)./var_ens_x(3,:)));
fprintf(fid, '%f, ', nanmedian(var_x_hat_fs(4,1:480)./var_ens_x(4,:)));
fprintf(fid, '%f, ', nanmedian(var_x_hat_fs(5,1:480)./var_ens_x(5,:)));
fprintf(fid, '%f\n', nanmedian(totVar_ens_hat_x(1:480)./totVar_ens_x(1:480)));

fclose(fid);

fprintf(fid2, '%.4f - %.4f - %.4f, %f, ', min(Parameter_sets(c_par_i,:)), mean(diff(Parameter_sets(c_par_i,:))),  ...
    max(Parameter_sets(c_par_i,:)), nanmedian(var_x_hat_fs(1,1:480)-var_ens_x(1,:)));
fprintf(fid2, '%f, ', nanmedian(var_x_hat_fs(2,1:480)-var_ens_x(2,:)));
fprintf(fid2, '%f, ', nanmedian(var_x_hat_fs(3,1:480)-var_ens_x(3,:)));
fprintf(fid2, '%f, ', nanmedian(var_x_hat_fs(4,1:480)-var_ens_x(4,:)));
fprintf(fid2, '%f, ', nanmedian(var_x_hat_fs(5,1:480)-var_ens_x(5,:)));
fprintf(fid2, '%f\n', nanmedian(totVar_ens_hat_x(1:480)-totVar_ens_x(1:480)));

fclose(fid2);

% %% Non-linearity Analysis
% figure;
% set(gcf, 'Position', [440,-152,1507,950]);
% 
% subplot(3,2,1);
% plotyy(1:480, var_ens_x(1,:), var_x_hat_fs(1,:), 'LineStyle', '-', 'Color', 'k');
% hold all;
% plot(var_x_hat_fs(1,:),'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
% text(380, (min(var_ens_x(1,:)) + 0.6*(max(var_ens_x(1,:))-min(var_ens_x(1,:)))), ['Linear Corr = ', num2str(round(corr(var_ens_x(1,:)',var_x_hat_fs(1,1:480)')*100)/100, '%.2f')], 'FontSize', 14);
% hleg = legend({'$x_1$', '$\hat{x}_1$'}, 'Interpreter', 'latex');
% set(hleg, 'fontsize',16);
% set(gca, 'FontSize', 12);
% ylabel('Total Variance', 'fontsize',18)
% xlabel('Time steps', 'FontSize', 14);
% title('Soil Moisture', 'FontSize', 16);
% 
% subplot(3,2,2);
% plot(var_ens_x(2,:), 'LineStyle', '-', 'Color', 'k');
% hold all;
% plot(var_x_hat_fs(2,:),'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
% text(380, (min(var_ens_x(2,:)) + 0.6*(max(var_ens_x(2,:))-min(var_ens_x(2,:)))), ['Linear Corr = ', num2str(round(corr(var_ens_x(2,:)',var_x_hat_fs(2,1:480)')*100)/100, '%.2f')], 'FontSize', 14);
% hleg = legend({'$x_2$', '$\hat{x}_2$'}, 'Interpreter', 'latex');
% set(hleg, 'fontsize',16);
% set(gca, 'FontSize', 12);
% ylabel('Total Variance', 'fontsize',18)
% xlabel('Time steps', 'FontSize', 14);
% title('Slow Reservoir', 'FontSize', 16);
% 
% subplot(3,2,3);
% plot(var_ens_x(3,:), 'LineStyle', '-', 'Color', 'k');
% hold all;
% plot(var_x_hat_fs(3,:),'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
% text(380, (min(var_ens_x(3,:)) + 0.6*(max(var_ens_x(3,:))-min(var_ens_x(3,:)))), ['Linear Corr = ', num2str(round(corr(var_ens_x(3,:)',var_x_hat_fs(3,1:480)')*100)/100, '%.2f')], 'FontSize', 14);
% hleg = legend({'$x_3$', '$\hat{x}_3$'}, 'Interpreter', 'latex');
% set(hleg, 'fontsize',16);
% set(gca, 'FontSize', 12);
% ylabel('Total Variance', 'fontsize',18);
% xlabel('Time steps', 'FontSize', 14);
% title('Quick Reservoir 1', 'FontSize', 16);
% 
% subplot(3,2,4);
% plot(var_ens_x(4,:), 'LineStyle', '-', 'Color', 'k');
% hold all;
% plot(var_x_hat_fs(4,:),'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
% text(380, (min(var_ens_x(4,:)) + 0.6*(max(var_ens_x(4,:))-min(var_ens_x(4,:)))), ['Linear Corr = ', num2str(round(corr(var_ens_x(4,:)',var_x_hat_fs(4,1:480)')*100)/100, '%.2f')], 'FontSize', 14);
% hleg = legend({'$x_4$', '$\hat{x}_4$'}, 'Interpreter', 'latex');
% set(hleg, 'fontsize',16);
% set(gca, 'FontSize', 12);
% ylabel('Total Variance', 'fontsize',18);
% xlabel('Time steps', 'FontSize', 14);
% title('Quick Reservoir 2', 'FontSize', 16);
% 
% subplot(3,2,5);
% plot(var_ens_x(5,:), 'LineStyle', '-', 'Color', 'k');
% hold all;
% plot(var_x_hat_fs(5,:),'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
% text(380, (min(var_ens_x(5,:)) + 0.6*(max(var_ens_x(5,:))-min(var_ens_x(5,:)))), ['Linear Corr = ', num2str(round(corr(var_ens_x(5,:)',var_x_hat_fs(5,1:480)')*100)/100, '%.2f')], 'FontSize', 14);
% hleg = legend({'$x_5$', '$\hat{x}_5$'}, 'Interpreter', 'latex');
% set(hleg, 'fontsize',16);
% set(gca, 'FontSize', 12);
% ylabel('Total Variance', 'fontsize',18);
% xlabel('Time steps', 'FontSize', 14);
% title('Quick Reservoir 3', 'FontSize', 16);
% 
% subplot(3,2,6);
% plot(totVar_ens_x, 'LineStyle', '-', 'Color', 'k');
% hold all;
% plot(totVar_ens_hat_x,'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
% text(380, (min(totVar_ens_x) + 0.6*(max(totVar_ens_x)-min(totVar_ens_x))), ['Linear Corr = ', num2str(round(corr(totVar_ens_x(1:480)',totVar_ens_hat_x(1:480)')*100)/100, '%.2f')], 'FontSize', 14);
% hleg = legend({'$x$', '$\hat{x}$'}, 'Interpreter', 'latex');
% set(hleg, 'fontsize',16);
% set(gca, 'FontSize', 12);
% ylabel('Total Variance', 'fontsize',18);
% xlabel('Time steps', 'FontSize', 14);
% title('Trace of Covariances', 'FontSize', 16);
% 
% % %% 4-D View of variance dynamics
% % for t = 1:1:480
% % %     subplot(1,2,1);
% %     scatter3(soil_moisture(:,t),quick_reservoirs(:,t,1),slow_reservoir(:,t),10,'k', 'filled');
% %     xlabel('Soil Moisture'); zlabel('Slow Reservoir'); ylabel('Quick Reservoir 1');
% %     set(gca, 'Xlim', [60 210], 'Zlim', [15 40], 'Ylim', [0 35]);
% %     
% % %     subplot(1,2,2);
% % %     scatter3(reshape(pertXu(1,:,t),size(pertXu,2),1)',reshape(pertXu(2,:,t),size(pertXu,2),1)',reshape(pertXu(3,:,t),size(pertXu,2),1)',(reshape(pertXu(4,:,t),size(pertXu,2),1)'+10).^2,reshape(pertXu(5,:,t),size(pertXu,2),1)', 'filled');
% % %     xlabel('Soil Moisture'); ylabel('Slow Reservoir'); zlabel('Quick Reservoir 1');
% % %     set(gca, 'Xlim', [-55 0], 'Ylim', [-12 8], 'Zlim', [-15 10]);
% % %     title(['t = ', num2str(t)]);
% % %     colorbar; caxis([-5 4])
% %     getframe();
% % end