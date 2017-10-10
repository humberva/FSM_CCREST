clear; clc;
load('Calibration_started_20130827214217_finished_20130901005456.mat');
par_names = {'PKE', 'PIM', 'PWM', 'PFC', 'LEAKO', 'LEAKI', 'PB', 'Nq'};

%Sum of Squared Residuals
objF = abs(prediction.EvalPars(:,10));

% for i = 1:7
%     subplot(4,2,i);
%     %Sample based on objective function threshold
% %     sample = prediction.EvalPars(objF>= quantile(objF,0.25) & objF <= quantile(objF,0.75),i);
%     sample = prediction.EvalPars(objF <= quantile(objF,0.3),i);
%     %Fit PDFs
%     [gamma_p,gamma_ci] = gamfit(sample);
%     [muhat,sigmahat,muci,sigmaci] = normfit(sample);
%     [parmhat,parmci] = wblfit(sample);
%     
%     sample_min = min(sample);
%     sample_max = max(sample);
%     sample_range = sample_max - sample_min;
%     [n,xout] = hist(sample,20);
%     bar(xout,(n./length(sample)).*100, 'hist');
% %     hist(sample);
% %     set(gca, 'Xlim', [sample_min sample_max], 'Ylim', [0 100], 'FontSize', 12);
%     ylabel('Frequency (%)', 'FontSize', 12);
% %     xlabel('Parameter Values', 'FontSize', 12);
%     title(par_names{i}, 'FontSize', 14);
%     
%     hold all;
%     gamma_sample = gamrnd(gamma_p(1),gamma_p(2),1,length(sample));
%     [n_g,xout_g] = hist(gamma_sample,xout);
%     plot(xout,(n_g./length(sample)).*100, 'Color', 'r', 'LineWidth', 2);
%     
%     gauss_sample = normrnd(muhat,sigmahat,1,length(sample));
%     [n_n,xout_g] = hist(gauss_sample,xout);
%     plot(xout,(n_n./length(sample)).*100, 'Color', 'b', 'LineWidth', 2);
%     
% %     wbull_sample = wblrnd(parmhat(1),parmhat(2),1,length(sample));
% %     [n_w,xout_g] = hist(wbull_sample,xout);
% %     plot(xout,(n_w./length(sample)).*100, 'Color', 'g', 'LineWidth', 2);
% end

sample = prediction.EvalPars(objF <= quantile(objF,0.8),1:8);
ref_set = [0.5 0.5 200 2 0.3 0.1 10 3];

n = 100;
p = 8;

inflt = [1 1 1 1 1 1 1 1];

random_noise = randn(n,p);

covmat = cov(sample);

%Eigen-decomposition of covariance matrix
[COVeigenvec,COVeigenval] = eig(covmat);

%Loop through perturbations
ensemble_data.perturbed_sets = zeros(n,p);
for i = 1:n
    perturbations = zeros(1,p);
    pnoise = zeros(1,p);
    for j = 1:p
        pnoise(j) = inflt(p).*random_noise(i,j)*sqrt(COVeigenval(j,j));
        perturbations = perturbations + (pnoise(j).*COVeigenvec(:,j)');
    end
    
    %Output data
    ensemble_data.perturbations(i,:) = perturbations;
    ensemble_data.perturbed_sets(i,:) = ref_set+perturbations;
end

for i = 1:p
    subplot(4,2,i);
    hist(ensemble_data.perturbed_sets(:,i));
end

% %Check for feasibility of proposed members (if feasible bounds provided)
% if (isfield(rundata, 'feasible_ranges') == 1)
%     for i = 1:n
%         minvalues = ensemble_data.perturbed_sets(i,:) - rundata.feasible_ranges(1,:);
%         maxvalues = rundata.feasible_ranges(2,:) - ensemble_data.perturbed_sets(i,:);
%         if (min(minvalues) < 0 || min(maxvalues) < 0)
%             fprintf('HPro Tool Warning: At least one of the ensemble members contains values outside feasible bounds of parameters.\nConsider modifying ensemble specs.\n');
%             ensemble_data = [];
%             break;
%         end
%     end
% end