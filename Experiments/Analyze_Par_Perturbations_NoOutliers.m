clear;
clc;

load('Ensemble_Simulation_Parameters_Perturbation_Analysis_n1000ens.mat');

% load('../Selected_Ensembles_of_States_n1000.mat')

%% 
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


%% 
totVar_pertUk = nan(1,size(pertXu,3));
totVar_xk = nan(1,size(pertXu,3));
mean_xi = nan(5,size(pertXu,3));
x1_perts = nan(1000,size(pertXu,3));
x2_perts = nan(1000,size(pertXu,3));
x3_perts = nan(1000,size(pertXu,3));
x4_perts = nan(1000,size(pertXu,3));
x5_perts = nan(1000,size(pertXu,3));

%Normalize by stddev
std_pertXu = pertXu;
nostd_pertXu = pertXu;

%Find bounds based on streamflow to ignore "non-behavioral" solutions
low_b = prctile(streamflow,40);
up_b = prctile(streamflow,60);
bounded_members = zeros(1000,size(pertXu,3));
for t = 1:size(pertXu,3)-1
    good_idxs = find(streamflow(:,t) > low_b(t) & streamflow(:,t) < up_b(t));
    bounded_members(good_idxs,t) = 1;
end

for t = 1:size(pertXu,3)-1
    good_idxs = find(bounded_members(:,t)==1);
    n_good_idxs = numel(good_idxs);
    for x = 1:5
        std_pertXu(x,good_idxs,t) = reshape(pertXu(x,good_idxs,t),n_good_idxs,1)./nanstd(reshape(pertXu(x,good_idxs,t),n_good_idxs,1));
    end
    totVar_pertUk(t) = trace(nancov(reshape(pertXu(:,good_idxs,t),n_good_idxs,5)));
    totVar_xk(t) = trace(nancov([soil_moisture(good_idxs,t), slow_reservoir(good_idxs,t),reshape(quick_reservoirs(good_idxs,t,:),n_good_idxs,3)]));
    mean_xi(:,t) = nanmean([soil_moisture(good_idxs,t), slow_reservoir(good_idxs,t),reshape(quick_reservoirs(good_idxs,t,:),n_good_idxs,3)])';
    x1_perts(good_idxs,t) = (soil_moisture(good_idxs,t) - mean_xi(1,t))./nanstd(soil_moisture(good_idxs,t));
    x2_perts(good_idxs,t) = (slow_reservoir(good_idxs,t) - mean_xi(2,t))./nanstd(slow_reservoir(good_idxs,t));
    x3_perts(good_idxs,t) = (quick_reservoirs(good_idxs,t,1) - mean_xi(3,t))./nanstd(quick_reservoirs(good_idxs,t,1));
    x4_perts(good_idxs,t) = (quick_reservoirs(good_idxs,t,2) - mean_xi(4,t))./nanstd(quick_reservoirs(good_idxs,t,2));
    x5_perts(good_idxs,t) = (quick_reservoirs(good_idxs,t,3) - mean_xi(5,t))./nanstd(quick_reservoirs(good_idxs,t,3));
    
    var_pertVk_x1(t) = nanvar(reshape(pertXu(1,good_idxs,t),n_good_idxs,1));
    var_pertVk_x2(t) = nanvar(reshape(pertXu(2,good_idxs,t),n_good_idxs,1));
    var_pertVk_x3(t) = nanvar(reshape(pertXu(3,good_idxs,t),n_good_idxs,1));
    var_pertVk_x4(t) = nanvar(reshape(pertXu(4,good_idxs,t),n_good_idxs,1));
    var_pertVk_x5(t) = nanvar(reshape(pertXu(5,good_idxs,t),n_good_idxs,1));
end

pertXu = std_pertXu;

%% Ensemble Simulations
figure;
set(gcf, 'Position', [156,-152,2336,950]);
subplot(2,5,1)
% plot(soil_moisture', 'Color', [0.7 0.7 0.7]);
plot(x1_perts', 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('pert. $x_1$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Soil Moisture', 'FontSize', 16);

subplot(2,5,2)
% plot(slow_reservoir', 'Color', [0.7 0.7 0.7]);
plot(x2_perts', 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('pert. $x_2$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Slow Reservoir', 'FontSize', 16);

subplot(2,5,3)
% plot(quick_reservoirs(:,:,1)', 'Color', [0.7 0.7 0.7]);
plot(x3_perts', 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('pert. $x_3$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Quick Reservoir 1', 'FontSize', 16);

subplot(2,5,4)
% plot(quick_reservoirs(:,:,2)', 'Color', [0.7 0.7 0.7]);
plot(x4_perts', 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('pert. $x_4$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Quick Reservoir 2', 'FontSize', 16);

subplot(2,5,5)
% plot(quick_reservoirs(:,:,3)', 'Color', [0.7 0.7 0.7]);
plot(x5_perts', 'Color', [0.7 0.7 0.7]);
hold all;
plot(0:480,zeros(1,481), 'Color', 'k');
ylabel('pert. $x_5$', 'interpreter', 'latex', 'fontsize',20);
xlabel('Time Steps', 'FontSize', 14);
set(gca, 'Xlim', [0 480], 'FontSize', 14);
title('Quick Reservoir 3', 'FontSize', 16);

for x = 1:5
    subplot(2,5,x+5)
    plot(reshape(pertXu(x,:,:),1000,481)', 'Color', [0.7 0.7 0.7]);
    hold all;
    plot(0:480,zeros(1,481), 'Color', 'k');
    ylabstr = ['$\delta x_', num2str(x),'$'];
    ylabel(ylabstr, 'interpreter', 'latex', 'fontsize',20)
    set(gca, 'Xlim', [0 480], 'FontSize', 12);
    xlabel('Time Steps', 'FontSize', 14);
end

%% Non-linearity Analysis
figure;
set(gcf, 'Position', [440,-152,1507,950]);

subplot(3,2,1);
[AX,hp1,hp2] = plotyy(1:480,nanvar(soil_moisture),1:480,var_pertVk_x1(1:480));
text(380, (min(nanvar(soil_moisture)) + 0.6*(max(nanvar(soil_moisture))-min(nanvar(soil_moisture)))), ['Linear Corr = ', num2str(round(corr(nanvar(soil_moisture)',var_pertVk_x1(1:480)')*100)/100, '%.2f')], 'FontSize', 14);
set(hp1, 'LineStyle', '-', 'Color', 'k');
set(hp2, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
hleg = legend([hp1 hp2], 'x_1', '\delta x_1');
set(hleg, 'fontsize',16);
set(AX(1),  'YColor', 'k', 'FontSize', 12);
set(AX(2),  'YColor', 'k', 'FontSize', 12);
ylabRstr = '$\delta x_1$ Variance';
ylabLstr = '$x_1$ Variance';
ylabel(AX(1), ylabLstr, 'interpreter', 'latex', 'fontsize',18)
ylabel(AX(2), ylabRstr, 'interpreter', 'latex', 'FontSize', 18);
xlabel('Time steps', 'FontSize', 14);
title('Soil Moisture', 'FontSize', 16);

subplot(3,2,2);
[AX,hp1,hp2] = plotyy(1:480,nanvar(slow_reservoir),1:480,var_pertVk_x2(1:480));
text(380, (min(nanvar(slow_reservoir)) + 0.6*(max(nanvar(slow_reservoir))-min(nanvar(slow_reservoir)))), ['Linear Corr = ', num2str(round(corr(nanvar(slow_reservoir)',var_pertVk_x2(1:480)')*100)/100, '%.2f')], 'FontSize', 14);
set(hp1, 'LineStyle', '-', 'Color', 'k');
set(hp2, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
hleg = legend([hp1 hp2], 'x_2', '\delta x_2');
set(hleg, 'fontsize',16);
set(AX(1),  'YColor', 'k', 'FontSize', 12);
set(AX(2),  'YColor', 'k', 'FontSize', 12);
ylabRstr = '$\delta x_2$ Variance';
ylabLstr = '$x_2$ Variance';
ylabel(AX(1), ylabLstr, 'interpreter', 'latex', 'fontsize',18)
ylabel(AX(2), ylabRstr, 'interpreter', 'latex', 'FontSize', 18);
xlabel('Time steps', 'FontSize', 14);
title('Slow Reservoir', 'FontSize', 16);

subplot(3,2,3);
[AX,hp1,hp2] = plotyy(1:480,nanvar(quick_reservoirs(:,:,1)),1:480,var_pertVk_x3(1:480));
text(380, (min(nanvar(quick_reservoirs(:,:,1))) + 0.6*(max(nanvar(quick_reservoirs(:,:,1)))-min(nanvar(quick_reservoirs(:,:,1))))), ['Linear Corr = ', num2str(round(corr(nanvar(quick_reservoirs(:,:,1))',var_pertVk_x3(1:480)')*100)/100, '%.2f')], 'FontSize', 14);
set(hp1, 'LineStyle', '-', 'Color', 'k');
set(hp2, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
hleg = legend([hp1 hp2], 'x_3', '\delta x_3');
set(hleg, 'fontsize',16);
set(AX(1),  'YColor', 'k', 'FontSize', 12);
set(AX(2),  'YColor', 'k', 'FontSize', 12);
ylabRstr = '$\delta x_3$ Variance';
ylabLstr = '$x_3$ Variance';
ylabel(AX(1), ylabLstr, 'interpreter', 'latex', 'fontsize',18)
ylabel(AX(2), ylabRstr, 'interpreter', 'latex', 'FontSize', 18);
xlabel('Time steps', 'FontSize', 14);
title('Quick Reservoir 1', 'FontSize', 16);

subplot(3,2,4);
[AX,hp1,hp2] = plotyy(1:480,nanvar(quick_reservoirs(:,:,2)),1:480,var_pertVk_x4(1:480));
text(380, (min(nanvar(quick_reservoirs(:,:,2))) + 0.88*(max(nanvar(quick_reservoirs(:,:,2)))-min(nanvar(quick_reservoirs(:,:,2))))), ['Linear Corr = ', num2str(round(corr(nanvar(quick_reservoirs(:,:,2))',var_pertVk_x4(1:480)')*100)/100, '%.2f')], 'FontSize', 14);
set(hp1, 'LineStyle', '-', 'Color', 'k');
set(hp2, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
hleg = legend([hp1 hp2], 'x_4', '\delta x_4');
set(hleg, 'fontsize',16);
set(AX(1),  'YColor', 'k', 'FontSize', 12);
set(AX(2),  'YColor', 'k', 'FontSize', 12);
ylabRstr = '$\delta x_4$ Variance';
ylabLstr = '$x_4$ Variance';
ylabel(AX(1), ylabLstr, 'interpreter', 'latex', 'fontsize',18)
ylabel(AX(2), ylabRstr, 'interpreter', 'latex', 'FontSize', 18);
xlabel('Time steps', 'FontSize', 14);
title('Quick Reservoir 2', 'FontSize', 16);

subplot(3,2,5);
[AX,hp1,hp2] = plotyy(1:480,nanvar(quick_reservoirs(:,:,3)),1:480,var_pertVk_x5(1:480));
text(380, (min(nanvar(quick_reservoirs(:,:,3))) + 0.6*(max(nanvar(quick_reservoirs(:,:,3)))-min(nanvar(quick_reservoirs(:,:,3))))), ['Linear Corr = ', num2str(round(corr(nanvar(quick_reservoirs(:,:,3))',var_pertVk_x5(1:480)')*100)/100, '%.2f')], 'FontSize', 14);
set(hp1, 'LineStyle', '-', 'Color', 'k');
set(hp2, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
hleg = legend([hp1 hp2], 'x_5', '\delta x_5');
set(hleg, 'fontsize',16);
set(AX(1),  'YColor', 'k', 'FontSize', 12);
set(AX(2),  'YColor', 'k', 'FontSize', 12);
ylabRstr = '$\delta x_5$ Variance';
ylabLstr = '$x_5$ Variance';
ylabel(AX(1), ylabLstr, 'interpreter', 'latex', 'fontsize',18)
ylabel(AX(2), ylabRstr, 'interpreter', 'latex', 'FontSize', 18);
xlabel('Time steps', 'FontSize', 14);
title('Quick Reservoir 3', 'FontSize', 16);

subplot(3,2,6);
[AX,hp1,hp2] = plotyy(1:480,totVar_xk(1:end-1),1:480,totVar_pertUk(1:end-1));
text(380, (min((totVar_xk(1:end-1))) + 0.6*(max((totVar_xk(1:end-1)))-min((totVar_xk(1:end-1))))), ['Linear Corr = ', num2str(round(corr((totVar_xk(1:end-1))',totVar_pertUk(1:end-1)')*100)/100, '%.2f')], 'FontSize', 14);
set(hp1, 'LineStyle', '-', 'Color', 'k');
set(hp2, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2);
hleg = legend([hp1 hp2], 'x_5', '\delta x_5');
set(hleg, 'fontsize',16);
set(AX(1),  'YColor', 'k', 'FontSize', 12);
set(AX(2),  'YColor', 'k', 'FontSize', 12);
ylabRstr = '$\delta x$ Total Variance';
ylabLstr = '$x$ Total Variance';
ylabel(AX(1), ylabLstr, 'interpreter', 'latex', 'fontsize',18)
ylabel(AX(2), ylabRstr, 'interpreter', 'latex', 'FontSize', 18);
xlabel('Time steps', 'FontSize', 14);
title('Trace of Covariances', 'FontSize', 16);

% %% 4-D View of variance dynamics
% for t = 1:1:480
% %     subplot(1,2,1);
%     scatter3(soil_moisture(:,t),quick_reservoirs(:,t,1),slow_reservoir(:,t),10,'k', 'filled');
%     xlabel('Soil Moisture'); zlabel('Slow Reservoir'); ylabel('Quick Reservoir 1');
%     set(gca, 'Xlim', [60 210], 'Zlim', [15 40], 'Ylim', [0 35]);
%     
% %     subplot(1,2,2);
% %     scatter3(reshape(pertXu(1,:,t),1000,1)',reshape(pertXu(2,:,t),1000,1)',reshape(pertXu(3,:,t),1000,1)',(reshape(pertXu(4,:,t),1000,1)'+10).^2,reshape(pertXu(5,:,t),1000,1)', 'filled');
% %     xlabel('Soil Moisture'); ylabel('Slow Reservoir'); zlabel('Quick Reservoir 1');
% %     set(gca, 'Xlim', [-55 0], 'Ylim', [-12 8], 'Zlim', [-15 10]);
% %     title(['t = ', num2str(t)]);
% %     colorbar; caxis([-5 4])
%     getframe();
% end