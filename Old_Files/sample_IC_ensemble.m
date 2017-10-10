%Sample and Analyze States Covariance
clear;
clc;
load('Simulation_2002-2009.mat');

%Drop 1-year
prediction.SM = prediction.SM(8761:end);
prediction.ISU = prediction.ISU(8761:end);
prediction.ISO = prediction.ISO(:,8761:end);
prediction.Q = prediction.Q(8761:end);

sample_n = 1000;

States = prediction.SM;
States = [States; prediction.ISU];
States = [States; prediction.ISO];

meanStates = mean(States,2);

normStates = States;
normStates(1,:) = normStates(1,:)-meanStates(1);
normStates(2,:) = normStates(2,:)-meanStates(2);
normStates(3,:) = normStates(3,:)-meanStates(3);
normStates(4,:) = normStates(4,:)-meanStates(4);
normStates(5,:) = normStates(5,:)-meanStates(5);

subStates = nan(size(States,1),sample_n);

min_coeff_var = 2;
min_trace = 600;
max_trace = 300;

for trial_i = 1:5000
    randomsample = randperm(numel(prediction.ISU)-1);
    randomsample = randomsample(prediction.Q(randomsample) > 0);
    randomsample = randomsample(1:sample_n);
%     randomsample(1,:) = randperm(numel(prediction.ISU),sample_n);
%     randomsample(2,:) = randperm(numel(prediction.ISU),sample_n);
%     randomsample(3,:) = randperm(numel(prediction.ISU),sample_n);
%     randomsample(4,:) = randperm(numel(prediction.ISU),sample_n);
%     randomsample(5,:) = randperm(numel(prediction.ISU),sample_n);
    
%     subStates(1,:) = States(1,randomsample(1,:));
%     subStates(2,:) = States(2,randomsample(2,:));
% %     subStates(3:end,:) = States(3:end,randomsample(3,:));
%     subStates(3,:) = States(3,randomsample(3,:));
%     subStates(4,:) = States(4,randomsample(4,:));
%     subStates(5,:) = States(5,randomsample(5,:));
    
    subStates = States(:,randomsample);
    
    eigenvals = eig(corr(subStates(1:end,:)'));
    eigenvals_covar = eig(cov(subStates(1:end,:)'));
%     CovTrace(trial_i) = trace(cov(subStates'));
    figure(1)
    plot(sort(eigenvals, 'descend'), 'Color', [0.7 0.7 0.7]);
    hold all;
    
    figure(2)
    plot(sort(eigenvals_covar, 'descend'), 'Color', [0.7 0.7 0.7]);
    hold all;
    
    coeff_var(trial_i) = std(eigenvals)/mean(eigenvals);
    trace_cov(trial_i) = sum(eigenvals_covar);
    pct_var(trial_i,:) = cumsum(sort(eigenvals_covar, 'descend'))./trace_cov(trial_i);
    pct_norm_var(trial_i,:) = cumsum(sort(eigenvals, 'descend'))./5;
    
    if (trace_cov(trial_i) <= min_trace)
        min_trace = trace_cov(trial_i);
        minVar_States = subStates;
        minVar_eigenvals = eigenvals;
        minVar_eigenvals_cov = eigenvals_covar;
        minVar_randomsample = randomsample;
    end
    
    if (trace_cov(trial_i) >= max_trace)
        max_trace = trace_cov(trial_i);
        maxVar_States = subStates;
        maxVar_eigenvals = eigenvals;
        maxVar_eigenvals_cov = eigenvals_covar;
        maxVar_randomsample = randomsample;
    end
    
    if (coeff_var(trial_i) < min_coeff_var)
        min_coeff_var = coeff_var(trial_i);
        final_States = subStates;
        final_eigenvals = eigenvals;
        final_eigenvals_cov = eigenvals_covar;
        final_randomsample = randomsample;
    end
    
%    P = polyfit(1:5,log(eigenvals'),1);
%    parfit1(trial_i) = P(1);
%    parfit2(trial_i) = P(2);
end
% figure(1);
f1h1 = plot(sort(final_eigenvals, 'descend'), 'DisplayName', 'Least Coeff. Var', 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
f1h2 = plot(sort(maxVar_eigenvals, 'descend'),'DisplayName', 'Max Cov Trace', 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
f1h3 = plot(sort(minVar_eigenvals, 'descend'),'DisplayName', 'Min Cov Trace', 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
set(gca, 'XTick', 1:1:5, 'FontSize', 14)
ylabel('Correlation Eigenvalue', 'FontSize', 14);
grid on;
legend([f1h1 f1h2 f1h3]);
saveas(gcf, '5000_Samples_Corr_EigenValues_Max_Min_Trace_1k_members.fig');
close;

% figure(2);
f2h1 = plot(sort(final_eigenvals_cov, 'descend'), 'DisplayName', 'Least Coeff. Var', 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);
f2h2 = plot(sort(maxVar_eigenvals_cov, 'descend'),'DisplayName', 'Max Cov Trace', 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2);
f2h3 = plot(sort(minVar_eigenvals_cov, 'descend'),'DisplayName', 'Min Cov Trace', 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
set(gca, 'XTick', 1:1:5, 'FontSize', 14)
ylabel('Covariance Eigenvalue', 'FontSize', 14);
grid on;
legend([f2h1 f2h2 f2h3]);
saveas(gcf, '5000_Samples_Covar_EigenValues_Max_Min_Trace_1k_members.fig');
close;

% save('Selected_Ensembles_of_States.mat', 'max_trace', 'maxVar_States', 'maxVar_eigenvals', 'maxVar_eigenvals_cov', 'maxVar_randomsample', 'min_trace', 'minVar_States', 'minVar_eigenvals', 'minVar_eigenvals_cov', 'minVar_randomsample', 'min_coeff_var', 'final_States', 'final_eigenvals_cov', 'final_eigenvals', 'final_randomsample');